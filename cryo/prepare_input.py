import argparse
import json
import os
import shutil
from pathlib import Path


def load_bioimage(path):
    try:
        from bioio import BioImage
        import numpy as np
    except ImportError as exc:
        raise RuntimeError(
            "bioio and numpy are required to convert CZI inputs to OME-TIFF"
        ) from exc

    image = BioImage(path)
    data = image.get_image_data("TCZYX")
    return np.asarray(data)


## display window percentiles: min/max is useless when a few pixels are saturated
## (common in the cryo CZIs), which made the images open almost black
CONTRAST_LOW_PERCENTILE = 0.1
CONTRAST_HIGH_PERCENTILE = 99.5


def channel_percentiles(channel, low, high):
    """(min, max, p_low, p_high) of one channel, exact via a histogram for integer data."""
    import numpy as np

    if np.issubdtype(channel.dtype, np.integer) and channel.size > 0:
        min_value = int(channel.min())
        max_value = int(channel.max())
        if min_value >= 0 and max_value < 2**24:
            cdf = np.cumsum(np.bincount(channel.ravel()))
            total = cdf[-1]
            p_low = int(np.searchsorted(cdf, total * low / 100.0))
            p_high = int(np.searchsorted(cdf, total * high / 100.0))
            return float(min_value), float(max_value), float(p_low), float(p_high)

    finite = channel[np.isfinite(channel)]
    if finite.size == 0:
        return None
    p_low, p_high = np.percentile(finite, [low, high])
    return float(finite.min()), float(finite.max()), float(p_low), float(p_high)


def channel_stats_from_tczyx(data):
    stats = []
    if data.ndim != 5:
        return stats

    for channel_index in range(data.shape[1]):
        values = channel_percentiles(
            data[:, channel_index, :, :, :],
            CONTRAST_LOW_PERCENTILE,
            CONTRAST_HIGH_PERCENTILE,
        )
        if values is None:
            continue
        min_value, max_value, p_low, p_high = values
        if p_high <= p_low:
            p_low, p_high = min_value, max_value
        stats.append(
            {
                "index": channel_index,
                "min": min_value,
                "max": max_value,
                "p_low": p_low,
                "p_high": p_high,
                "contrast_limits": f"({p_low:.6g},{p_high:.6g})",
            }
        )
    return stats


def ome_tiff_metadata(metadata):
    ome_metadata = {"axes": "TCZYX"}
    scale_keys = {
        "PhysicalSizeX": "x_scale_nm",
        "PhysicalSizeY": "y_scale_nm",
        "PhysicalSizeZ": "z_scale_nm",
    }
    for ome_key, metadata_key in scale_keys.items():
        value = metadata.get(metadata_key)
        if value is not None:
            ome_metadata[ome_key] = float(value)
            ome_metadata[f"{ome_key}Unit"] = "nm"

    channels = metadata.get("channels", [])
    if channels:
        channel_metadata = {
            "Name": [
                channel.get("display") or channel.get("label") or f"channel_{index}"
                for index, channel in enumerate(channels)
            ]
        }
        fluor = [channel.get("fluor") for channel in channels]
        excitation = [channel.get("excitation_wavelength_nm") for channel in channels]
        emission = [channel.get("emission_wavelength_nm") for channel in channels]
        if all(value is not None for value in fluor):
            channel_metadata["Fluor"] = fluor
        if all(value is not None for value in excitation):
            channel_metadata["ExcitationWavelength"] = excitation
            channel_metadata["ExcitationWavelengthUnit"] = ["nm"] * len(channels)
        if all(value is not None for value in emission):
            channel_metadata["EmissionWavelength"] = emission
            channel_metadata["EmissionWavelengthUnit"] = ["nm"] * len(channels)
        ome_metadata["Channel"] = channel_metadata
    return ome_metadata


def write_czi_as_ome_tiff(input_path, output_path, metadata):
    try:
        import tifffile
    except ImportError as exc:
        raise RuntimeError("tifffile is required to write intermediate OME-TIFF") from exc

    data = load_bioimage(input_path)
    channel_stats = channel_stats_from_tczyx(data)
    tifffile.imwrite(
        output_path,
        data,
        bigtiff=True,
        ome=True,
        photometric="minisblack",
        metadata=ome_tiff_metadata(metadata),
    )
    return channel_stats


def link_or_copy(source, target):
    try:
        os.symlink(source, target)
    except OSError:
        shutil.copy2(source, target)


def main():
    parser = argparse.ArgumentParser(
        description="Prepare CRYO raw inputs for EuBI-Bridge conversion."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--name", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--metadata-json", required=True)
    parser.add_argument("--prepared-log", required=True)
    args = parser.parse_args()

    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    if not input_path.exists():
        raise FileNotFoundError(f"Raw image does not exist: {input_path}")

    output_dir.mkdir(parents=True, exist_ok=True)
    suffix = input_path.suffix.lower()
    metadata = {}
    metadata_path = Path(args.metadata_json)
    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))

    if suffix == ".czi":
        prepared_path = output_dir / f"{args.name}.ome.tif"
        channel_stats = write_czi_as_ome_tiff(input_path, prepared_path, metadata)
        mode = "czi_to_ome_tiff"
    else:
        prepared_path = output_dir / input_path.name
        link_or_copy(input_path, prepared_path)
        channel_stats = []
        mode = "link_or_copy_original"

    metadata["prepared_input_path"] = str(prepared_path)
    metadata["prepared_input_mode"] = mode
    if channel_stats:
        metadata["channel_stats"] = channel_stats
        stats_by_index = {row["index"]: row for row in channel_stats}
        for channel in metadata.get("channels", []):
            stats = stats_by_index.get(channel.get("index"))
            if stats is not None:
                for key in ("min", "max", "p_low", "p_high", "contrast_limits"):
                    channel[key] = stats[key]
    metadata_path.write_text(
        json.dumps(metadata, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    with Path(args.prepared_log).open("w", encoding="utf-8") as handle:
        handle.write("name\traw_path\tprepared_path\tmode\n")
        handle.write(f"{args.name}\t{input_path}\t{prepared_path}\t{mode}\n")


if __name__ == "__main__":
    main()
