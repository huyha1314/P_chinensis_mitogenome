import math
import os
import random
from dataclasses import dataclass

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO

random.seed(7)
np.random.seed(7)

SPECIES_NAME = "Parashorea chinensis"
INPUTS = [
    ("data/mt_VN1_contig1.gb", "1"),
    ("data/mt_VN1_contig2.gb", "2"),
]

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 12,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

COLOR_DICT = {
    "complex I (NADH dehydrogenase)": "#f2d64b",
    "complex II (succinate dehydrogenase)": "#4caf50",
    "complex III (ubichinol cytochrome c reductase)": "#cddc39",
    "complex IV (cytochrome c oxidase)": "#f48fb1",
    "ATP synthase": "#8bc34a",
    "cytochrome c biogenesis": "#8e24aa",
    "RNA polymerase": "#d32f2f",
    "ribosomal proteins (SSU)": "#d7ccc8",
    "ribosomal proteins (LSU)": "#8d6e63",
    "maturases": "#ff9800",
    "other genes": "#ce93d8",
    "ORFs": "#00bcd4",
    "transfer RNAs": "#3f51b5",
    "ribosomal RNAs": "#f44336",
    "origin of replication": "#ffffff",
    "polycistronic transcripts": "#eeeeee",
}

OUTER_FEATURE_BOTTOM = 0.92  
OUTER_FEATURE_TOP = 0.98     
INNER_FEATURE_BOTTOM = 0.85  
INNER_FEATURE_TOP = 0.91     
GC_BOTTOM = 0.38             
GC_TOP = 0.58                


@dataclass(frozen=True)
class MapFeature:
    name: str
    feature_type: str
    segments: tuple[tuple[int, int], ...]
    strand: int
    is_pseudo: bool

    @property
    def start(self) -> int:
        return min(start for start, _end in self.segments)

    @property
    def end(self) -> int:
        return max(end for _start, end in self.segments)

    @property
    def center(self) -> float:
        return (self.start + self.end) / 2.0


def feature_segments(feature) -> tuple[tuple[int, int], ...]:
    location = feature.location
    parts = getattr(location, "parts", [location])
    segments = []
    for part in parts:
        start = int(part.start)
        end = int(part.end)
        if end > start:
            segments.append((start, end))
    return tuple(segments)


def is_pseudogene(feature) -> bool:
    return (
        "pseudo" in feature.qualifiers
        or "pseudogene" in feature.qualifiers
        or "pseudogene" in " ".join(feature.qualifiers.get("note", [])).lower()
    )


def collect_features(record) -> list[MapFeature]:
    primary = []
    represented_genes = set()
    seen_segments = set()

    gene_name_map = {}
    for feature in record.features:
        if feature.type == "gene" and "gene" in feature.qualifiers:
            start = int(feature.location.start)
            end = int(feature.location.end)
            gene_name_map[(start, end)] = feature.qualifiers["gene"][0]

    for feature in record.features:
        if feature.type not in {"CDS", "tRNA", "rRNA"}:
            continue
        
        name = feature.qualifiers.get("gene", [""])[0]
        if not name:
            start = int(feature.location.start)
            end = int(feature.location.end)
            name = gene_name_map.get((start, end), "")
        if not name:
            name = feature.qualifiers.get("product", [""])[0]

        segments = feature_segments(feature)
        if not name or not segments:
            continue
            
        primary.append(MapFeature(name, feature.type, segments, feature.location.strand or 1, is_pseudogene(feature)))
        represented_genes.add(name.lower())
        seen_segments.add(segments)

    for feature in record.features:
        if feature.type != "gene":
            continue
            
        name = feature.qualifiers.get("gene", [""])[0]
        if not name:
            name = feature.qualifiers.get("product", [""])[0]
            
        segments = feature_segments(feature)
        
        if not name or not segments or name.lower() in represented_genes or segments in seen_segments:
            continue
            
        primary.append(MapFeature(name, feature.type, segments, feature.location.strand or 1, True))
        represented_genes.add(name.lower())
        seen_segments.add(segments)

    return sorted(primary, key=lambda item: (item.start, item.end, item.name))


def get_smart_labels(features: list[MapFeature]) -> list[MapFeature]:
    MAX_GAP = 5000
    
    split_features = []
    for f in features:
        segs = sorted(f.segments, key=lambda s: s[0])
        if not segs: continue
        
        clusters = []
        curr = [segs[0]]
        for s in segs[1:]:
            if s[0] - curr[-1][1] > MAX_GAP:
                clusters.append(curr)
                curr = [s]
            else:
                curr.append(s)
        clusters.append(curr)
        
        for c in clusters:
            split_features.append(MapFeature(f.name, f.feature_type, tuple(c), f.strand, f.is_pseudo))
            
    final_labels = []
    grouped = {}
    for f in split_features:
        grouped.setdefault(f.name.lower(), []).append(f)
        
    for name, feats in grouped.items():
        feats.sort(key=lambda x: x.start)
        cluster = [feats[0]]
        for f in feats[1:]:
            if f.start - cluster[-1].end <= MAX_GAP:
                cluster.append(f)
            else:
                merged_segments = tuple(s for feat in cluster for s in feat.segments)
                final_labels.append(MapFeature(
                    cluster[0].name, cluster[0].feature_type,
                    merged_segments, cluster[0].strand, cluster[0].is_pseudo
                ))
                cluster = [f]
        
        merged_segments = tuple(s for feat in cluster for s in feat.segments)
        final_labels.append(MapFeature(
            cluster[0].name, cluster[0].feature_type,
            merged_segments, cluster[0].strand, cluster[0].is_pseudo
        ))
        
    return sorted(final_labels, key=lambda f: f.start)


def functional_category(feature: MapFeature) -> tuple[str, str]:
    gene = feature.name.lower()
    if gene.startswith(("nad", "ndh")): category = "complex I (NADH dehydrogenase)"
    elif gene.startswith("sdh"): category = "complex II (succinate dehydrogenase)"
    elif gene.startswith("cob"): category = "complex III (ubichinol cytochrome c reductase)"
    elif gene.startswith("cox"): category = "complex IV (cytochrome c oxidase)"
    elif gene.startswith("atp"): category = "ATP synthase"
    elif gene.startswith("ccm"): category = "cytochrome c biogenesis"
    elif gene.startswith("rpo"): category = "RNA polymerase"
    elif gene.startswith("rps"): category = "ribosomal proteins (SSU)"
    elif gene.startswith("rpl"): category = "ribosomal proteins (LSU)"
    elif gene.startswith("trn") or feature.feature_type == "tRNA": category = "transfer RNAs"
    elif gene.startswith("rrn") or feature.feature_type == "rRNA": category = "ribosomal RNAs"
    elif gene.startswith("mat"): category = "maturases"
    elif gene.startswith("orf"): category = "ORFs"
    else: category = "other genes"
    return category, COLOR_DICT[category]


def theta(position: float, genome_size: int) -> float:
    return (position / genome_size) * 2.0 * math.pi


def repel_positions(positions: list[float], genome_size: int, min_gap: float) -> np.ndarray:
    if len(positions) < 2:
        return np.array(positions, dtype=float)

    adjusted = np.array(positions, dtype=float)
    n = len(adjusted)

    for _ in range(8000):
        moved = False
        for i in range(n):
            next_i = (i + 1) % n
            diff = (adjusted[next_i] - adjusted[i]) % genome_size
            if diff < min_gap:
                shift = min((min_gap - diff) / 2.0, min_gap * 0.05)
                adjusted[i] = (adjusted[i] - shift) % genome_size
                adjusted[next_i] = (adjusted[next_i] + shift) % genome_size
                moved = True
        if not moved:
            break
    return adjusted


def draw_arc_block(ax, genome_size: int, start: int, end: int, bottom: float, top: float, color: str):
    width = theta(end - start, genome_size)
    center = theta((start + end) / 2.0, genome_size)
    ax.bar(
        center, top - bottom, width=width, bottom=bottom,
        align="center", color=color, edgecolor="black", linewidth=0.45, zorder=4,
    )


def draw_gc_track(ax, record, window_size: int = 200, step_size: int = 50):
    sequence = str(record.seq).upper()
    genome_size = len(sequence)
    positions = []
    values = []

    for start in range(0, genome_size, step_size):
        window = sequence[start : min(start + window_size, genome_size)]
        called = window.count("A") + window.count("T") + window.count("G") + window.count("C")
        gc = ((window.count("G") + window.count("C")) / called) if called else 0.0
        positions.append(start + len(window) / 2.0)
        values.append(gc)

    avg_gc = float(np.mean(values))
    track_height = GC_TOP - GC_BOTTOM
    heights = np.array(values) * track_height
    avg_radius = GC_TOP - (avg_gc * track_height)
    angles = [theta(pos, genome_size) for pos in positions]

    angles.append(angles[0] + 2 * math.pi)
    heights = np.append(heights, heights[0])

    ax.bar(0, GC_TOP - GC_BOTTOM, width=2 * math.pi, bottom=GC_BOTTOM, color="#eeeeee", zorder=1)
    
    ax.fill_between(
        angles,
        GC_TOP - heights,
        GC_TOP,
        color="#bdbdbd",
        linewidth=0, 
        zorder=2,
    )
    ax.plot(angles, [avg_radius] * len(angles), color="#8f8f8f", linewidth=0.9, zorder=3)


def feature_track_limits(feature: MapFeature) -> tuple[float, float]:
    if feature.strand == -1: return INNER_FEATURE_BOTTOM, INNER_FEATURE_TOP
    return OUTER_FEATURE_BOTTOM, OUTER_FEATURE_TOP


def text_rotation_and_alignment(label_theta: float, align_outside: bool) -> tuple[float, str]:
    display_degrees = (90.0 - math.degrees(label_theta)) % 360.0
    is_bottom_half = 90.0 < display_degrees < 270.0
    
    if align_outside:
        ha = "right" if is_bottom_half else "left"
    else:
        ha = "left" if is_bottom_half else "right"
        
    rotation = display_degrees + 180.0 if is_bottom_half else display_degrees
    return rotation, ha


def draw_labels(ax, features: list[MapFeature], genome_size: int):
    outer_features = [f for f in features if f.strand == 1]
    inner_features = [f for f in features if f.strand == -1]

    def _draw_strand_labels(feats, leader_start, leader_end, label_radius, align_outside, fontsize=10.5):
        n = len(feats)
        if n == 0: return
        
        feats = sorted(feats, key=lambda f: f.center)
        positions = [f.center for f in feats]
        
        if align_outside:
            label_gap = genome_size * (3.8 / 360.0) 
        else:
            label_gap = genome_size * (4.8 / 360.0) 
            
        label_positions = repel_positions(positions, genome_size, label_gap)
        isolated_gap = genome_size * (6.0 / 360.0)

        for i, (feature, original_pos, label_pos) in enumerate(zip(feats, positions, label_positions)):
            label_theta = theta(label_pos, genome_size)
            rotation, ha = text_rotation_and_alignment(label_theta, align_outside)

            if n > 1:
                dist_prev = (original_pos - positions[(i - 1) % n]) % genome_size
                dist_next = (positions[(i + 1) % n] - original_pos) % genome_size
                is_crowded = (dist_prev < isolated_gap) or (dist_next < isolated_gap)
            else:
                is_crowded = False
            
            moved = min((label_pos - original_pos) % genome_size, (original_pos - label_pos) % genome_size)
            if moved > genome_size * 0.005: 
                is_crowded = True

            needs_leader = (len(feature.segments) > 1) or is_crowded

            if needs_leader:
                actual_label_radius = label_radius
                for s_start, s_end in feature.segments:
                    seg_center = (s_start + s_end) / 2.0
                    seg_theta = theta(seg_center, genome_size)

                    delta = label_theta - seg_theta
                    if delta > math.pi: delta -= 2 * math.pi
                    elif delta < -math.pi: delta += 2 * math.pi

                    t = np.linspace(0.0, 1.0, 36)
                    smooth = t * t * (3.0 - 2.0 * t)
                    curve_theta = seg_theta + delta * smooth
                    curve_radius = leader_start + (leader_end - leader_start) * t
                    
                    ax.plot(curve_theta, curve_radius, color="#777777", lw=0.55, alpha=0.78, zorder=5)
            else:
                actual_label_radius = leader_start + 0.015 if align_outside else leader_start - 0.015

            ax.text(
                label_theta, actual_label_radius, feature.name,
                rotation=rotation, ha=ha, va="center",
                rotation_mode="anchor", fontsize=fontsize, color="#111111", zorder=6,
            )

    # ĐÃ SỬA: Giảm độ dài râu nối xuống, làm chúng ngắn hơn đáng kể
    _draw_strand_labels(outer_features, OUTER_FEATURE_TOP, 1.030, 1.035, align_outside=True, fontsize=16)
    _draw_strand_labels(inner_features, INNER_FEATURE_BOTTOM, 0.820, 0.815, align_outside=False, fontsize=15)


def draw_scale_ticks(ax, genome_size: int):
    step = 50_000
    ticks = list(range(0, genome_size + 1, step))
    for value in ticks:
        angle = theta(value, genome_size)
        ax.plot([angle, angle], [0.915, 0.935], color="black", linewidth=0.7, zorder=5)
        label = "0" if value == 0 else f"{value // 1000} kb"
        degrees = math.degrees(angle) % 360
        rotation = degrees + 180 if 90 < degrees < 270 else degrees
        ax.text(
            angle, 0.955, label, rotation=rotation, ha="center",
            va="center", rotation_mode="anchor", fontsize=7.8, color="#222222", zorder=6,
        )


def draw_genome_map(ax, gb_file: str, chromosome: str):
    record = SeqIO.read(gb_file, "genbank")
    genome_size = len(record.seq)
    features = collect_features(record)

    ax.set_theta_direction(-1)
    ax.set_theta_offset(math.pi / 2.0)
    ax.set_ylim(0, 1.22)
    ax.grid(False)
    ax.set_axis_off()

    ax.bar(0, 0.008, width=2 * math.pi, bottom=0.912, color="black", zorder=3)
    draw_gc_track(ax, record)

    for feature in features:
        _category, color = functional_category(feature)
        feature_bottom, feature_top = feature_track_limits(feature)
        for start, end in feature.segments:
            draw_arc_block(ax, genome_size, start, end, feature_bottom, feature_top, color)

    label_features = get_smart_labels(features)
    draw_labels(ax, label_features, genome_size)

    species_parts = SPECIES_NAME.split(" ", 1)
    if len(species_parts) == 2:
        species = rf"$\mathit{{{species_parts[0]}\ {species_parts[1]}}}$"
    else:
        species = SPECIES_NAME
    ax.text(
        0, 0,
        f"{species}\nmitochondrion chromosome {chromosome}\ncomplete sequence\n{genome_size:,} bp",
        ha="center", va="center", fontsize=15, linespacing=1.25,
    )

    return {functional_category(feature)[0] for feature in features}


def main():
    os.makedirs("results", exist_ok=True)
    fig = plt.figure(figsize=(24, 13), facecolor="white") 
    axes = [fig.add_subplot(1, 2, i + 1, projection="polar") for i in range(2)]

    used_categories = set()
    for ax, (gb_file, chromosome) in zip(axes, INPUTS):
        used_categories |= draw_genome_map(ax, gb_file, chromosome)

    legend_handles = [
        mpatches.Patch(facecolor=COLOR_DICT[label], edgecolor="black", label=label)
        for label in COLOR_DICT
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center", ncol=4, fontsize=18, bbox_to_anchor=(0.5, 0.035),
        frameon=False, columnspacing=1.35, labelspacing=0.65, handleheight=1.05, handlelength=1.05,
    )
    fig.subplots_adjust(left=0.035, right=0.975, top=0.965, bottom=0.20, wspace=0.18)

    out_png = "results/Fig3_Circular_Genomic_Map.png"
    out_pdf = "results/Fig3_Circular_Genomic_Map.pdf"
    fig.savefig(out_png, dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(out_pdf, bbox_inches="tight", facecolor="white")
    print(f"Saved {out_png}")
    print(f"Saved {out_pdf}")


if __name__ == "__main__":
    main()