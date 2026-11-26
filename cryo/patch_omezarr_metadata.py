import argparse
import json
import shutil
import xml.etree.ElementTree as ET
from pathlib import Path


COLOR_TO_HEX = {
    "red": "FF0000",
    "green": "00FF00",
    "blue": "0000FF",
    "yellow": "FFFF00",
    "magenta": "FF00FF",
    "cyan": "00FFFF",
    "white": "FFFFFF",
    "gray": "FFFFFF",
    "grey": "FFFFFF",
}
def image_channels(metadata, channel_count=None):
    """Channel metadata in the order of the image's C axis.

    Every channel is kept, in its original index order, so that entry i describes
    image channel i. (A former filter kept only GFP/PE/ChloA/TL in that fixed order,
    which dropped e.g. AF568/Reflection and put names on the wrong channels.)
    channel_count is the C size of the written image; the list is cut or padded to it.
    """
    channels = sorted(
        (dict(channel) for channel in metadata.get("channels", [])),
        key=lambda channel: int(channel.get("index", 0)),
    )
    if channel_count is not None:
        channels = channels[:channel_count]
        for index in range(len(channels), channel_count):
            channels.append({"label": f"channel_{index}", "display": f"channel_{index}"})
    for index, channel in enumerate(channels):
        channel["index"] = index
    return channels


def find_zarr_root(path):
    root = Path(path)
    if is_zarr_root(root):
        return root

    for candidate in root.rglob("*"):
        if candidate.is_dir() and candidate.suffix in {".zarr", ".ome.zarr"}:
            if is_zarr_root(candidate):
                return candidate
    raise FileNotFoundError(f"Could not find OME-Zarr root below {root}")


def is_zarr_root(path):
    return (path / "zarr.json").exists() or (path / ".zattrs").exists() or (path / ".zgroup").exists()


def load_root_attributes(zarr_root):
    """Root group attributes of an OME-Zarr.

    zarr v3 (OME-Zarr 0.5): zarr.json -> attributes, OME metadata under attributes["ome"]
    zarr v2 (OME-Zarr 0.4): .zattrs, OME metadata (multiscales, omero) at the top level

    Returns (file to write, full document, custom attributes, OME attributes).
    """
    zarr_json = zarr_root / "zarr.json"
    if zarr_json.exists():
        document = json.loads(zarr_json.read_text(encoding="utf-8"))
        attributes = document.setdefault("attributes", {})
        return zarr_json, document, attributes, attributes.setdefault("ome", {})

    zattrs_path = zarr_root / ".zattrs"
    attrs = {}
    if zattrs_path.exists():
        attrs = json.loads(zattrs_path.read_text(encoding="utf-8"))
    return zattrs_path, attrs, attrs, attrs


def image_channel_count(zarr_root, ome_attrs):
    """C size of the full-resolution array (1 when eubi squeezed a single channel away)."""
    try:
        multiscale = ome_attrs["multiscales"][0]
        axes = [axis["name"] if isinstance(axis, dict) else axis for axis in multiscale["axes"]]
        array_dir = zarr_root / multiscale["datasets"][0]["path"]
    except (KeyError, IndexError, TypeError):
        return None
    for name in ("zarr.json", ".zarray"):
        array_meta = array_dir / name
        if array_meta.exists():
            shape = json.loads(array_meta.read_text(encoding="utf-8")).get("shape")
            break
    else:
        return None
    if not shape or len(shape) != len(axes):
        return None
    return int(shape[axes.index("c")]) if "c" in axes else 1


def parse_contrast_limits(value):
    if not value:
        return None, None
    value = str(value).strip().strip("()")
    parts = [part.strip() for part in value.split(",")]
    if len(parts) != 2:
        return None, None
    try:
        return float(parts[0]), float(parts[1])
    except ValueError:
        return None, None


def channel_color(channel):
    color = str(channel.get("color", "")).strip()
    return COLOR_TO_HEX.get(color.lower(), "FFFFFF")


def build_omero_channels(channels):
    omero_channels = []
    for index, channel in enumerate(channels):
        ## window start/end = display range (percentiles), min/max = data range
        start, end = channel.get("p_low"), channel.get("p_high")
        if start is None or end is None:
            start, end = parse_contrast_limits(channel.get("contrast_limits"))
        min_value = channel.get("min", start)
        max_value = channel.get("max", end)
        if start is None:
            start = min_value
        if end is None:
            end = max_value
        if start is None or end is None:
            start, end = 0.0, 1.0
        if min_value is None or max_value is None:
            min_value, max_value = start, end

        label = (
            channel.get("display")
            or channel.get("label")
            or f"channel_{channel.get('index', index)}"
        )
        omero_channels.append(
            {
                "active": True,
                "coefficient": 1,
                "color": channel_color(channel),
                "family": "linear",
                "inverted": False,
                "label": str(label),
                "window": {
                    "start": start,
                    "end": end,
                    "min": min_value,
                    "max": max_value,
                },
            }
        )
    return omero_channels


def compact_cryo_metadata(metadata, channels):
    compact = dict(metadata)
    if "channels" in compact:
        compact["channels"] = channels
        compact["size_c"] = len(channels)
    return compact


def xml_local_name(element):
    return str(element.tag).split("}", 1)[-1]


def format_ome_float(value):
    return f"{float(value):g}"


def patch_channel_attribute(element, attribute, value):
    if value is None:
        element.attrib.pop(attribute, None)
        return
    element.set(attribute, format_ome_float(value))


def patch_ome_xml_file(path, channels):
    try:
        tree = ET.parse(path)
    except ET.ParseError:
        return False

    root = tree.getroot()
    namespace = ""
    if root.tag.startswith("{"):
        namespace = root.tag.split("}", 1)[0][1:]
        ET.register_namespace("", namespace)

    if not channels:
        return False

    channel_elements = [
        element
        for element in root.iter()
        if xml_local_name(element) == "Channel"
    ]
    if not channel_elements:
        return False

    for index, element in enumerate(channel_elements[: len(channels)]):
        channel = channels[index]
        label = channel.get("display") or channel.get("label")
        fluor = channel.get("fluor")
        if label:
            element.set("Name", str(label))
        if fluor:
            element.set("Fluor", str(fluor))
        else:
            element.attrib.pop("Fluor", None)

        excitation = channel.get("excitation_wavelength_nm")
        emission = channel.get("emission_wavelength_nm")
        patch_channel_attribute(element, "ExcitationWavelength", excitation)
        patch_channel_attribute(element, "EmissionWavelength", emission)
        if excitation is None:
            element.attrib.pop("ExcitationWavelengthUnit", None)
        else:
            element.set("ExcitationWavelengthUnit", "nm")
        if emission is None:
            element.attrib.pop("EmissionWavelengthUnit", None)
        else:
            element.set("EmissionWavelengthUnit", "nm")

    tree.write(path, encoding="utf-8", xml_declaration=True)
    return True


def patch_ome_xml_files(zarr_root, channels):
    patched = []
    for path in sorted(zarr_root.rglob("*.ome.xml")):
        if patch_ome_xml_file(path, channels):
            patched.append(str(path))
    return patched


def copy_czi_metadata_xml(zarr_root, czi_metadata_xml):
    if not czi_metadata_xml:
        return None
    source = Path(czi_metadata_xml)
    if not source.exists() or source.stat().st_size == 0:
        return None

    metadata_dir = zarr_root / "czi_metadata"
    metadata_dir.mkdir(parents=True, exist_ok=True)
    target = metadata_dir / "czi_metadata.xml"
    shutil.copy2(source, target)
    return target


def main():
    parser = argparse.ArgumentParser(
        description="Patch extracted CRYO/CZI metadata into the root attributes of an OME-Zarr (zarr.json or .zattrs)."
    )
    parser.add_argument("--omezarr", required=True)
    parser.add_argument("--metadata-json", required=True)
    parser.add_argument("--czi-metadata-xml", default="")
    parser.add_argument("--log", required=True)
    args = parser.parse_args()

    zarr_root = find_zarr_root(args.omezarr)
    zattrs_path, document, attrs, ome_attrs = load_root_attributes(zarr_root)

    metadata = json.loads(Path(args.metadata_json).read_text(encoding="utf-8"))
    channels = image_channels(metadata, image_channel_count(zarr_root, ome_attrs))
    omero_channels = build_omero_channels(channels)
    if omero_channels:
        omero = ome_attrs.get("omero", {})
        omero["channels"] = omero_channels
        omero.setdefault("rdefs", {})["model"] = "color"
        ome_attrs["omero"] = omero

    attrs["cryo_metadata"] = compact_cryo_metadata(metadata, channels)
    zattrs_path.write_text(json.dumps(document, indent=2, sort_keys=True), encoding="utf-8")
    patched_ome_xmls = patch_ome_xml_files(zarr_root, channels)
    copied_czi_metadata = copy_czi_metadata_xml(zarr_root, args.czi_metadata_xml)

    with Path(args.log).open("w", encoding="utf-8") as handle:
        handle.write("omezarr\tmetadata_json\tpatched_zattrs\tchannels\tpatched_ome_xmls\tczi_metadata_xml\n")
        handle.write(
            f"{zarr_root}\t{args.metadata_json}\t{zattrs_path}\t"
            f"{len(omero_channels)}\t{len(patched_ome_xmls)}\t"
            f"{copied_czi_metadata or ''}\n"
        )


if __name__ == "__main__":
    main()
