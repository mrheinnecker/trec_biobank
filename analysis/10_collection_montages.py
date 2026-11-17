#!/usr/bin/env python3
"""Overview figures (montages) of the collection tables in <run>/tables/collections/, true to scale:
  <A>_vs_<B>.tsv (09_negative_pair_collections.R): groups <A>_only, <B>_only, both, none
  clusters_k<k>.tsv (06_organelle_clustering.R):   groups C1 .. Ck

For every table one JPEG <run>/montages_f<scale_factor>/<table>.jpg: the groups stacked from top to
bottom, each under a header band. All cells are shown at the same physical scale (big cells big,
small cells small).

For clusters_k<k> additionally <table>_tree.jpg: the clustering tree on top (branches coloured by
cluster; from tables/clust_tree_k<k>_*.tsv of step 06) and below one column of cells per cluster,
in tree order. The width of a cluster's column is proportional to its number of
cells (one strip of equal width per cell, as in the clustering heatmap), so the tree is not distorted
and columns end at different heights. The cell images are reduced (same factor for all cells, so still
true to scale) until the widest cell of every cluster fits its column; --tree_image_scale is the upper
limit of that factor. Between tree and cells, as in the clustering heatmap of step 06:
the metadata tracks (cluster, site, size_frac, sampling_time, life_status; legends on top) and the
presence/absence of every organelle used for the clustering (rows grouped, label with % of cells),
one strip per cell above the column with its image.

Scale: the overview PNGs of the TEM workflow all have the same pixel size (MRC pixel 1.766 nm x
preview factor 4 = 7.064 nm; the burned-in scale bar is 5 um). Every PNG is reduced by the same
--scale_factor (default 4 -> 28.3 nm per montage pixel; full resolution would exceed the
JPEG size limit), so the burned-in scale bars all have the
same length.

Layout: rectangle packing with the skyline bottom-left heuristic: the cells of a group are sorted
by height (tallest first) and each is put at the lowest position on the current skyline where it
fits (leftmost on ties), so small cells fill the gaps next to big ones. All groups share one width;
the width is chosen by bisection so that the stacked figure is about square.

Images: <png_dir>/<site>/<source_name without .zarr>.png, reduced once to grayscale JPEG thumbnails
cached in --thumb_dir/f<scale_factor>/. Cells without a PNG are listed in montages_f<scale_factor>/missing_pngs.tsv.

usage: python3 10_collection_montages.py --outdir <run folder> [--png_dir DIR] [--thumb_dir DIR]
       [--scale_factor 4] [--quality 85] [--workers 6] [--tables 'clusters_*']
"""

import argparse
import csv
import fnmatch
import re
import os
import sys
from multiprocessing import Pool

from PIL import Image, ImageDraw, ImageFont

Image.MAX_IMAGE_PIXELS = None  # the stitched TEM images are far above the default limit

PNG_PIXEL_NM = 1.766 * 4     # pixel size of the overview PNGs (MRC pixel x preview factor)
GAP = 8                      # px between tiles
BACKGROUND = 255             # white
FONT_PATHS = ["/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
              "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", required=True, help="run folder with tables/collections/")
    p.add_argument("--png_dir", default=os.environ.get("TEM_PNG_DIR", "/mnt/c/projects/tem_classification/pngs"))
    p.add_argument("--thumb_dir", default=os.environ.get("TEM_THUMB_DIR", "/mnt/c/projects/tem_classification/png_thumbs"))
    p.add_argument("--scale_factor", type=int, default=int(os.environ.get("TEM_SCALE_FACTOR", 4)),
                   help="PNG pixels per montage pixel (same for all cells -> true scale)")
    p.add_argument("--quality", type=int, default=85, help="JPEG quality of the montages")
    p.add_argument("--workers", type=int, default=6, help="parallel thumbnail jobs (each large PNG needs a few GB RAM)")
    p.add_argument("--tree_plot_scale", type=float, default=float(os.environ.get("TEM_TREE_PLOT_SCALE", 1.5)),
                   help="tree figure: size of tree, tracks, organelle rows, labels and legend (1 = as the montage headers)")
    p.add_argument("--tree_image_scale", type=float, default=float(os.environ.get("TEM_TREE_IMAGE_SCALE", 0.5)),
                   help="tree figure: largest size of the cell images relative to the montages (smaller if needed so that "
                        "the widest cell of every cluster fits its column; same for all cells -> true scale)")
    p.add_argument("--tables", default=os.environ.get("TEM_MONTAGE_TABLES", "*"),
                   help="only these collection tables (glob on the file name without .tsv, e.g. 'clusters_*')")
    # run_analysis.sh passes the options of the R steps too; they are not used here
    p.add_argument("--input")
    p.add_argument("--rounds")
    return p.parse_args()


def font(size, bold=False):
    try:
        return ImageFont.truetype(FONT_PATHS[1 if bold else 0], size)
    except OSError:
        return ImageFont.load_default()


# ---- thumbnails ----
def make_thumb(job):
    """PNG -> grayscale JPEG reduced by the common factor (cached); returns (name, thumb path or None, error)."""
    name, png, thumb, factor = job
    if os.path.exists(thumb) and os.path.getmtime(thumb) >= os.path.getmtime(png):
        return name, thumb, None
    try:
        with Image.open(png) as im:
            if im.mode in ("RGBA", "LA", "P"):
                im = im.convert("RGBA")
                bg = Image.new("RGBA", im.size, (255, 255, 255, 255))
                im = Image.alpha_composite(bg, im)
            im = im.convert("L")
            size = (max(1, round(im.width / factor)), max(1, round(im.height / factor)))
            if factor >= 4:
                im = im.reduce(factor // 2)            # fast integer step first
            im = im.resize(size, Image.LANCZOS)
            os.makedirs(os.path.dirname(thumb), exist_ok=True)
            im.save(thumb, "JPEG", quality=90)
        return name, thumb, None
    except Exception as e:  # noqa: BLE001 - report and continue with the other cells
        return name, None, str(e)


# ---- packing ----
def skyline_pack(sizes, width):
    """Skyline bottom-left packing of rectangles (w, h) into a strip of the given width.
    Returns (positions [(x, y)] in input order, used height). Rectangles are placed tallest first."""
    order = sorted(range(len(sizes)), key=lambda i: (-sizes[i][1], -sizes[i][0]))
    skyline = [[0, 0, width]]   # segments [x, y, w], left to right, covering [0, width)
    pos = [None] * len(sizes)
    for i in order:
        w, h = sizes[i][0] + GAP, sizes[i][1] + GAP
        best = None             # (top, y, x, segment index)
        for s in range(len(skyline)):
            x = skyline[s][0]
            if x + w > width:
                break
            y, covered, k = 0, 0, s
            while covered < w:  # highest segment under [x, x + w)
                y = max(y, skyline[k][1])
                covered += skyline[k][2] if k > s else skyline[k][2]
                k += 1
            if best is None or (y + h, x) < (best[0], best[2]):
                best = (y + h, y, x, s)
        if best is None:
            raise ValueError("rectangle wider than the strip")
        top, y, x, s = best
        pos[i] = (x, y)
        # new segment [x, x + w) at height top; cut the covered segments
        new = []
        for seg in skyline:
            sx, sy, sw = seg
            ex = sx + sw
            if ex <= x or sx >= x + w:
                new.append(seg)
                continue
            if sx < x:
                new.append([sx, sy, x - sx])
            if ex > x + w:
                new.append([x + w, sy, ex - (x + w)])
        new.append([x, top, w])
        new.sort(key=lambda g: g[0])
        merged = []             # merge neighbours of the same height
        for seg in new:
            if merged and merged[-1][1] == seg[1] and merged[-1][0] + merged[-1][2] == seg[0]:
                merged[-1][2] += seg[2]
            else:
                merged.append(seg)
        skyline = merged
    height = max((pos[i][1] + sizes[i][1] for i in range(len(sizes))), default=0)
    return pos, height


def figure_height(groups, width, header):
    return sum(header + skyline_pack(sizes, width)[1] + 2 * GAP for sizes in groups)


def choose_width(groups, header):
    """Common strip width of all groups so that the stacked figure is about square (bisection)."""
    widest = max((w for sizes in groups for w, _ in sizes), default=1) + GAP
    area = sum((w + GAP) * (h + GAP) for sizes in groups for w, h in sizes)
    lo, hi = widest, max(widest, int(area ** 0.5 * 2) + 1)
    for _ in range(25):
        if hi - lo <= 16:
            break
        mid = (lo + hi) // 2
        if figure_height(groups, mid, header) > mid:
            lo = mid
        else:
            hi = mid
    return hi


# ---- drawing ----
def labelled_tile(thumb, name, label_font, size=None):
    im = Image.open(thumb).convert("L")
    if size and size != im.size:
        im = im.resize(size, Image.LANCZOS)
    draw = ImageDraw.Draw(im)
    tb = draw.textbbox((0, 0), name, font=label_font)
    th = tb[3] - tb[1]
    draw.rectangle([0, im.height - th - 8, min(im.width, tb[2] - tb[0] + 8), im.height], fill=0)
    draw.text((4, im.height - th - 6 - tb[1]), name, fill=255, font=label_font)   # clipped to the tile
    return im


def render(groups, width, header, title, nm_per_px):
    """groups: list of (label, [(name, thumb, w, h)]) -> PIL image."""
    packed = [skyline_pack([(t[2], t[3]) for t in tiles], width) if tiles else ([], 0) for _, tiles in groups]
    total = header + sum(header + h + 2 * GAP for _, h in packed)
    canvas = Image.new("L", (width + 2 * GAP, total), BACKGROUND)
    draw = ImageDraw.Draw(canvas)

    # title + scale bar (10 um)
    draw.text((GAP, header * 0.15), title, fill=0, font=font(int(header * 0.55), bold=True))
    bar_px = round(10000 / nm_per_px)
    bx, by = canvas.width - GAP * 4 - bar_px, int(header * 0.35)
    draw.rectangle([bx, by, bx + bar_px, by + int(header * 0.2)], fill=0)
    draw.text((bx - int(header * 0.1), by), "10 µm", fill=0, font=font(int(header * 0.3)), anchor="ra")

    label_font = font(22)
    y = header
    for (label, tiles), (positions, h) in zip(groups, packed):
        draw.rectangle([0, y, canvas.width, y + header - GAP], fill=40)
        draw.text((GAP * 3, y + header * 0.2), f"{label}  ({len(tiles)} cells)", fill=255,
                  font=font(int(header * 0.45), bold=True))
        for (name, thumb, _, _), (px, py) in zip(tiles, positions):
            canvas.paste(labelled_tile(thumb, name, label_font), (GAP + px, y + header + py))
        y += header + h + 2 * GAP
    return canvas


# ---- tree + cluster columns ----
def hex_rgb(h):
    h = h.lstrip("#")
    return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))


TRACKS = ["cluster", "site", "size_frac", "sampling_time", "life_status"]
TOP_TRACKS = ["cluster", "site", "size_frac", "sampling_time"]   # between tree and organelles
BOTTOM_TRACKS = ["life_status"]                                  # below the organelles, above taxo_class


def render_tree(clusters, tiles, leaves, merges, colors, title, nm_per_px, header, annot=None,
                plot_scale=1.0, image_scale=1.0):
    """clusters: C1..Ck in tree order; tiles[c] = [(name, thumb, w, h)]; leaves: name -> (pos, cluster, number);
    annot: (tracks {name: row}, rows [(feature, group, prevalence)]) or None;
    plot_scale: size of the plot part (tree, tracks, rows, labels); image_scale: size of the cell images"""
    if image_scale != 1:
        tiles = {c: [(n, t, max(1, round(w * image_scale)), max(1, round(h * image_scale))) for n, t, w, h in ts]
                 for c, ts in tiles.items()}
        nm_per_px = nm_per_px / image_scale
    header = int(header * plot_scale)
    gap = GAP * 6
    track_h, row_h = int(header * 1.2), int(header * 0.9)
    label_font, group_font = font(int(row_h * 0.75)), font(int(row_h * 0.75), bold=True)
    legend_font, legend_bold = font(int(header * 0.75)), font(int(header * 0.75), bold=True)
    tracks, rows = annot if annot else ({}, [])
    groups_in_rows = list(dict.fromkeys(r[1] for r in rows))
    row_labels = [f"{f}  ({100 * prev:.0f}%)" for f, _, prev in rows]
    # vertical space between title and tree: legend band, tracks, organelle rows (group gaps hold the group name)
    legend_h = 0   # the legend is on the right side
    tracks_h = (len(TOP_TRACKS) * track_h + GAP * 4) if annot else 0
    heat_h = (len(rows) * row_h + len(groups_in_rows) * int(row_h * 1.2) + GAP * 4) if annot else 0
    # taxonomic class (where annotated): extra track below the organelles + coloured frame around the image
    taxo = {n: (r.get("taxo_class", ""), r.get("taxo_class_color", "")) for n, r in tracks.items()
            if r.get("taxo_class", "")}
    bottom_h = (GAP * 8 + (len(BOTTOM_TRACKS) + (1 if taxo else 0)) * track_h) if annot else 0
    heat_h += bottom_h
    extra_h = legend_h + tracks_h + heat_h
    probe = ImageDraw.Draw(Image.new("L", (1, 1)))
    texts = TRACKS + row_labels + (["taxo_class"] if taxo else []) if annot else []
    # right side: row labels next to the strips, then the legend (one block per track, one value per line)
    label_w = max((probe.textbbox((0, 0), t, font=label_font)[2] for t in texts), default=0)
    legend = []
    if annot:
        def track_legend(tr):
            vals = {}
            for row in tracks.values():
                if row[tr] != "":
                    vals[row[tr]] = row[tr + "_color"]
            return [(tr, None)] + sorted(vals.items(),
                                         key=lambda kv: [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", kv[0])])
        for tr in TOP_TRACKS:   # legend in the order of the strips
            legend += track_legend(tr)
        legend += [("organelle", None), ("present", "#2b2b2b"), ("absent", "#f2f2f2")]
        for tr in BOTTOM_TRACKS:
            legend += track_legend(tr)
        if taxo:
            legend.append(("taxo_class (frame)", None))
            legend += sorted(set(taxo.values()))
    swatch = int(header * 0.6)
    legend_w = max((probe.textbbox((0, 0), t, font=legend_bold if c is None else legend_font)[2]
                    + (0 if c is None else swatch + GAP * 2) for t, c in legend), default=0)
    left = GAP * 2
    right = (GAP * 6 + label_w + GAP * 12 + legend_w + GAP * 6) if annot else GAP * 2

    # column width of a cluster proportional to its number of cells (one strip of equal width per cell,
    # as in the clustering heatmap); the image scale is then set by the column whose widest cell needs
    # the most room, so that every cell fits its column (image_scale is the upper limit)
    n_cells = {c: sum(1 for _, cl, _ in leaves.values() if cl == c) for c in clusters}
    strip = max(8, int(header * 0.35))       # px per cell
    widths = {c: n_cells[c] * strip for c in clusters}
    width = sum(widths.values()) + gap * (len(clusters) - 1)
    fit = min(((widths[c] - 2 * GAP) / max(t[2] for t in tiles[c]) for c in clusters if tiles[c]), default=1.0)
    if fit < 1:
        tiles = {c: [(n, t, max(1, int(w * fit)), max(1, int(h * fit))) for n, t, w, h in ts] for c, ts in tiles.items()}
        nm_per_px = nm_per_px / fit
    print(f"  tree figure: {strip} px per cell, images at {nm_per_px:.1f} nm/px", flush=True)
    packs = {c: skyline_pack([(t[2], t[3]) for t in tiles[c]], widths[c]) if tiles[c] else ([], 0) for c in clusters}
    tree_h = int(width * 0.16)
    total_h = header + extra_h + tree_h + header + max(h for _, h in packs.values()) + 2 * GAP

    canvas = Image.new("RGB", (left + width + right, total_h), (255, 255, 255))
    label_x = left + width + GAP * 6
    legend_x = label_x + label_w + GAP * 12
    draw = ImageDraw.Draw(canvas)
    draw.text((GAP, header * 0.15), title, fill=(0, 0, 0), font=font(int(header * 0.55), bold=True))
    bar_px = round(10000 / nm_per_px)
    bx, by = canvas.width - GAP * 4 - bar_px, int(header * 0.35)
    draw.rectangle([bx, by, bx + bar_px, by + int(header * 0.2)], fill=(0, 0, 0))
    draw.text((bx - int(header * 0.1), by), "10 µm", fill=(0, 0, 0), font=font(int(header * 0.3)), anchor="ra")

    # legend on the right: track name in bold, then one colour box + value per line
    if annot:
        ly = header + GAP * 4
        line_h = int(header * 0.95)
        for text, col in legend:
            if col is None:
                ly += int(line_h * 0.4) if ly > header + GAP * 4 else 0
                draw.text((legend_x, ly), text, fill=(0, 0, 0), font=legend_bold)
            else:
                draw.rectangle([legend_x, ly + int(line_h * 0.1), legend_x + swatch, ly + int(line_h * 0.1) + swatch],
                               fill=hex_rgb(col), outline=(150, 150, 150))
                draw.text((legend_x + swatch + GAP * 2, ly), text, fill=(0, 0, 0), font=legend_font)
            ly += line_h

    # column positions
    starts, x = {}, left + GAP
    for c in clusters:
        starts[c] = x
        x += widths[c] + gap

    # leaf x: spread over the column of the cluster, in tree order
    node_x, node_h, node_cl = {}, {}, {}
    for c in clusters:
        members = sorted((pos, name) for name, (pos, cl, _) in leaves.items() if cl == c)
        for r, (pos, name) in enumerate(members):
            key = -pos  # leaves are keyed by their position
            node_x[key] = starts[c] + (r + 0.5) / len(members) * widths[c]
            node_h[key] = 0.0
            node_cl[key] = {c}
    pos_of_leaf = {i + 1: pos for i, pos in enumerate(leaves_by_number(leaves))}
    max_h = max(m[3] for m in merges) or 1
    top, bottom = header + legend_h + GAP, header + legend_h + tree_h
    lw = max(3, round(width / 2500 * plot_scale))

    def y_of(h):
        return bottom - (bottom - top) * h / max_h

    def ref(v):
        return -pos_of_leaf[-v] if v < 0 else v

    for step, child_a, child_b, height in merges:
        a, b = ref(child_a), ref(child_b)
        node_x[step] = (node_x[a] + node_x[b]) / 2
        node_h[step] = height
        node_cl[step] = node_cl[a] | node_cl[b]
        for child in (a, b):
            col = hex_rgb(colors[next(iter(node_cl[child]))]) if len(node_cl[child]) == 1 else (70, 70, 70)
            draw.line([(node_x[child], y_of(node_h[child])), (node_x[child], y_of(height))], fill=col, width=lw)
        col = hex_rgb(colors[next(iter(node_cl[step]))]) if len(node_cl[step]) == 1 else (70, 70, 70)
        draw.line([(node_x[a], y_of(height)), (node_x[b], y_of(height))], fill=col, width=lw)

    # metadata tracks and organelle presence/absence: one strip per cell, above its column
    if annot:
        spans = {}
        for c in clusters:
            members = sorted((pos, name) for name, (pos, cl, _) in leaves.items() if cl == c)
            for r, (_, name) in enumerate(members):
                spans[name] = (starts[c] + r * widths[c] / len(members), starts[c] + (r + 1) * widths[c] / len(members))
        def draw_track(tr, y):
            draw.text((label_x, y + (track_h - row_h) / 2), tr, fill=(0, 0, 0), font=label_font)
            for name, (x0, x1) in spans.items():
                if name in tracks:
                    draw.rectangle([x0, y, max(x0, x1 - 1), y + track_h - 4], fill=hex_rgb(tracks[name][tr + "_color"]))

        y = bottom + GAP * 4
        for tr in TOP_TRACKS:
            draw_track(tr, y)
            y += track_h
        y += GAP * 4
        current = None
        for (feature, group, _), label in zip(rows, row_labels):
            if group != current:
                draw.text((label_x, y + int(row_h * 0.1)), group, fill=(90, 90, 90), font=group_font)
                y += int(row_h * 1.2)
                current = group
            draw.text((label_x, y), label, fill=(0, 0, 0), font=label_font)
            for name, (x0, x1) in spans.items():
                if name in tracks:
                    present = tracks[name].get(feature) == "1"
                    draw.rectangle([x0, y, max(x0, x1 - 1), y + row_h - 3],
                                   fill=(43, 43, 43) if present else (242, 242, 242))
            y += row_h
        y += GAP * 8
        for tr in BOTTOM_TRACKS:
            draw_track(tr, y)
            y += track_h
        if taxo:
            draw.text((label_x, y + (track_h - row_h) / 2), "taxo_class", fill=(0, 0, 0), font=label_font)
            for name, (x0, x1) in spans.items():
                col = hex_rgb(taxo[name][1]) if name in taxo else (242, 242, 242)
                draw.rectangle([x0, y, max(x0, x1 - 1), y + track_h - 4], fill=col)

    # cluster columns
    label_font = font(22)
    y0 = header + legend_h + tree_h + tracks_h + heat_h
    for c in clusters:
        n_leaves = sum(1 for _, cl, _ in leaves.values() if cl == c)
        draw.rectangle([starts[c], y0, starts[c] + widths[c], y0 + header - GAP], fill=hex_rgb(colors[c]))
        draw.text((starts[c] + GAP * 2, y0 + header * 0.2), f"{c}  ({n_leaves} cells)", fill=(255, 255, 255),
                  font=font(int(header * 0.45), bold=True))
        positions, _ = packs[c]
        frame = max(6, int(header * 0.12))
        for (name, thumb, w, h), (px, py) in zip(tiles[c], positions):
            x_, y_ = starts[c] + px, y0 + header + py
            canvas.paste(labelled_tile(thumb, name, label_font, (w, h)), (x_, y_))
            if name in taxo:   # coloured frame (inside the tile) for cells with a taxonomic class
                draw.rectangle([x_, y_, x_ + w - 1, y_ + h - 1], outline=hex_rgb(taxo[name][1]), width=frame)
    return canvas


def leaves_by_number(leaves):
    """positions of the leaves in the original hclust label order (leaf number 1..n)"""
    return [pos for _, (pos, _, _) in sorted(leaves.items(), key=lambda kv: kv[1][2])]


def read_tree(table_dir, k):
    lf = os.path.join(table_dir, f"clust_tree_k{k}_leaves.tsv")
    mf = os.path.join(table_dir, f"clust_tree_k{k}_merges.tsv")
    if not (os.path.exists(lf) and os.path.exists(mf)):
        return None
    with open(lf, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    # name -> (position in tree order, cluster, leaf number in hclust label order)
    leaves = {r["name"]: (int(r["leaf_pos"]), r["cluster"], i + 1) for i, r in enumerate(rows)}
    colors = {r["cluster"]: r["color"] for r in rows}
    with open(mf, newline="") as fh:
        merges = [(int(r["step"]), int(r["left"]), int(r["right"]), float(r["height"]))
                  for r in csv.DictReader(fh, delimiter="\t")]
    tf = os.path.join(table_dir, f"clust_tree_k{k}_tracks.tsv")
    rf = os.path.join(table_dir, f"clust_tree_k{k}_rows.tsv")
    annot = None
    if os.path.exists(tf) and os.path.exists(rf):
        with open(tf, newline="") as fh:
            tracks = {r["name"]: r for r in csv.DictReader(fh, delimiter="\t")}
        with open(rf, newline="") as fh:
            rows = [(r["feature"], r["group"], float(r["prevalence"])) for r in csv.DictReader(fh, delimiter="\t")]
        annot = (tracks, rows)
    return leaves, merges, colors, annot


def main():
    args = parse_args()
    coll_dir = os.path.join(args.outdir, "tables", "collections")
    tables = sorted(f for f in os.listdir(coll_dir)
                    if f.endswith(".tsv") and f != "overview.tsv" and fnmatch.fnmatch(f[:-4], args.tables)) \
        if os.path.isdir(coll_dir) else []
    if not tables:
        sys.exit(f"No collection tables matching '{args.tables}' in {coll_dir} (run 06 or 04 + 09 first)")
    out_dir = os.path.join(args.outdir, f"montages_f{args.scale_factor}")   # one folder per scale factor
    os.makedirs(out_dir, exist_ok=True)

    collections = {}
    for f in tables:
        with open(os.path.join(coll_dir, f), newline="") as fh:
            collections[f[:-4]] = list(csv.DictReader(fh, delimiter="\t"))

    # thumbnails of all cells of all tables, once
    thumb_root = os.path.join(args.thumb_dir, f"f{args.scale_factor}")
    cells, missing = {}, []
    for rows in collections.values():
        for r in rows:
            if r["name"] in cells:
                continue
            stem = r["source_name"].removesuffix(".zarr")
            png = os.path.join(args.png_dir, r["site"], stem + ".png")
            if os.path.exists(png):
                cells[r["name"]] = (png, os.path.join(thumb_root, r["site"], stem + ".jpg"))
            else:
                missing.append((r["name"], r["site"], png))
                cells[r["name"]] = None
    jobs = [(n, p[0], p[1], args.scale_factor) for n, p in cells.items() if p]
    print(f"{len(jobs)} cells with a PNG, {len(missing)} without; thumbnails in {thumb_root}", flush=True)
    thumbs = {}
    with Pool(args.workers) as pool:
        for k, (name, thumb, err) in enumerate(pool.imap_unordered(make_thumb, jobs), 1):
            if err:
                missing.append((name, "", f"thumbnail failed: {err}"))
            else:
                thumbs[name] = thumb
            if k % 50 == 0 or k == len(jobs):
                print(f"  thumbnails {k}/{len(jobs)}", flush=True)
    dims = {}
    for name, thumb in thumbs.items():
        with Image.open(thumb) as im:
            dims[name] = im.size

    with open(os.path.join(out_dir, "missing_pngs.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["name", "site", "detail"])
        w.writerows(missing)

    nm_per_px = PNG_PIXEL_NM * args.scale_factor
    header = 120
    def natural(g):  # C2 before C10
        return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", g)]

    for pair, rows in collections.items():
        present = list(dict.fromkeys(r["group"] for r in rows))
        if "_vs_" in pair:   # negative pair: fixed order, empty groups shown too
            a, b = pair.split("_vs_", 1)
            order = [f"{a}_only", f"{b}_only", "both", "none"]
            order += [g for g in present if g not in order]
        else:                # clusters and others: groups in natural order
            order = sorted(present, key=natural)
        groups = [(g, [(r["name"], thumbs[r["name"]], *dims[r["name"]])
                       for r in rows if r["group"] == g and r["name"] in thumbs]) for g in order]
        width = choose_width([[(t[2], t[3]) for t in tiles] for _, tiles in groups], header)
        img = render(groups, width, header, pair.replace("_vs_", "  vs  "), nm_per_px)
        if max(img.size) > 65500:
            sys.exit(f"{pair}: montage {img.size} exceeds the JPEG limit; use a larger --scale_factor")
        out = os.path.join(out_dir, pair + ".jpg")
        img.save(out, "JPEG", quality=args.quality, optimize=True)
        print(f"{pair}: {img.width} x {img.height} px ({nm_per_px:.1f} nm/px), "
              f"{os.path.getsize(out) / 1e6:.1f} MB -> {out}", flush=True)

        # clusters: tree on top, one column of cells per cluster below
        m = re.fullmatch(r"clusters_k(\d+)", pair)
        tree = read_tree(os.path.join(args.outdir, "tables"), m.group(1)) if m else None
        if m and tree is None:
            print(f"{pair}: no tree files of step 06, tree figure skipped", flush=True)
        if tree:
            leaves, merges, colors, annot = tree
            clusters = [g for g, _ in groups]
            tiles = dict(groups)
            img = render_tree(clusters, tiles, leaves, merges, colors,
                              f"{pair}: clustering tree and cells per cluster", nm_per_px, header, annot,
                              args.tree_plot_scale, args.tree_image_scale)
            if max(img.size) > 65500:
                sys.exit(f"{pair}: tree figure {img.size} exceeds the JPEG limit; use a larger --scale_factor")
            out = os.path.join(out_dir, pair + "_tree.jpg")
            img.save(out, "JPEG", quality=args.quality, optimize=True)
            print(f"{pair}_tree: {img.width} x {img.height} px, {os.path.getsize(out) / 1e6:.1f} MB -> {out}", flush=True)


if __name__ == "__main__":
    main()
