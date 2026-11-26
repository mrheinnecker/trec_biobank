"""Build a CRYO input table with the 10x overview of every grid from a file inventory.

The inventory is the output of (run on the parent folder of the site folders):

    find "$D" \\( -type f -o -type d \\) -printf '%y\\t%s\\t%TY-%Tm-%Td\\t%P\\n' | sort -t$'\\t' -k4

Layout: <year>-<site>/<sample>/.../<grid>_data/<file>.czi, e.g.
2022-Iceland/asIC01/asIC01a_data/asIC01a_10x-01.czi. The grid is taken from the folder,
because many files carry the name of another sample (the microscope kept the last name).

Per grid, the raw 10x acquisitions are the candidates. Left out: ZEN derivatives
(Orthogonal Projection, Split Scenes, Stitching, Shading Correction, Airyscan Processing),
manual _orth/_stitched exports, files < 15 MB (aborted acquisitions or the empty base
file of a multi-part series) and "-01" copies with the same size as their original.
Grids with exactly one candidate get convert=1, all other rows convert=0 with a note.
"""
import argparse
import csv
import re
from collections import defaultdict
from pathlib import PurePosixPath

DERIVED = re.compile(r"Orthogonal Projection|Split Scenes|Stitching|Shading Correction|Airyscan Processing")
EXPORT = re.compile(r"orth|stitched", re.IGNORECASE)
TEN_X = re.compile(r"(?<!\d)10x", re.IGNORECASE)
GRID_DIR = re.compile(r"^(?:as|asIC|WF)\d+[a-d]$")
MULTIPART = re.compile(r"\(\d+\)\.czi$", re.IGNORECASE)
MIN_SIZE = 15_000_000


def grid_of(parts):
    for folder in reversed(parts[:-1]):
        if "_data" in folder or GRID_DIR.match(folder):
            return re.sub(r"_data.*$", "", folder.replace("_airy", ""))
    return None


def sanitize(value):
    return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--inventory", required=True)
    parser.add_argument("--parent", required=True, help="absolute path of the inventoried folder on the cluster")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    files = []
    with open(args.inventory, encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) == 4 and row[0] == "f" and row[3].lower().endswith(".czi"):
                files.append((int(row[1]), row[3]))
    sizes = {path: size for size, path in files}

    candidates = defaultdict(list)
    for size, path in files:
        parts = PurePosixPath(path).parts
        name = parts[-1]
        if "tobediscarded" in path or not TEN_X.search(name):
            continue
        if DERIVED.search(name) or EXPORT.search(name) or size < MIN_SIZE:
            continue
        # "<x>-01.czi" next to "<x>.czi" with the same size: a ZEN copy
        original = re.sub(r"-\d{2}\.czi$", ".czi", path)
        if original != path and sizes.get(original) == size:
            continue
        grid = grid_of(parts)
        if grid is None:
            continue
        candidates[grid].append((size, path))

    rows = []
    for grid in sorted(candidates):
        options = sorted(candidates[grid], key=lambda item: item[1])
        for size, path in options:
            parts = PurePosixPath(path).parts
            year, _, site = parts[0].partition("-")
            stem = sanitize(PurePosixPath(path).stem)
            # dataset name always starts with the grid from the folder
            filename = stem if stem.startswith(grid) else f"{grid}_{stem}"
            notes = []
            blocked = len(options) > 1 or bool(MULTIPART.search(path))
            if len(options) > 1:
                notes.append(f"{len(options)} 10x candidates for this grid, choose")
            if MULTIPART.search(path):
                notes.append("multi-part CZI series, reading untested")
            if not PurePosixPath(path).name.startswith(grid):
                notes.append("file carries another sample name")
            rows.append(
                {
                    "FilePath": f"{args.parent.rstrip('/')}/{path}",
                    "filename": filename,
                    "site": site,
                    "year": year,
                    "sample": parts[1],
                    "grid_id": grid,
                    "acq_type": "confocal" if "confocal" in path.lower() else "widefield",
                    "size_mb": round(size / 1e6),
                    "convert": "0" if blocked else "1",
                    "note": "; ".join(notes),
                }
            )

    with open(args.output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    grids = {row["grid_id"] for row in rows}
    ready = [row for row in rows if row["convert"] == "1"]
    print(f"{len(rows)} rows, {len(grids)} grids, {len(ready)} with convert=1 "
          f"({sum(r['size_mb'] for r in ready) / 1000:.1f} GB)")


if __name__ == "__main__":
    main()
