import os
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pycirclize import Circos
from pycirclize.parser import Genbank

# Exact OGDRAW-like color mapping based on user's image
COLOR_DICT = {
    'complex I (NADH dehydrogenase)': '#ffeb3b', # yellow
    'complex II (succinate dehydrogenase)': '#4caf50', # green
    'complex III (ubichinol cytochrome c reductase)': '#cddc39', # lime
    'complex IV (cytochrome c oxidase)': '#f48fb1', # pink
    'ATP synthase': '#8bc34a', # olive green
    'cytochrome c biogenesis': '#8e24aa', # purple
    'RNA polymerase': '#d32f2f', # red
    'ribosomal proteins (SSU)': '#d7ccc8', # light brown
    'ribosomal proteins (LSU)': '#8d6e63', # brown
    'maturases': '#ff9800', # orange
    'transfer RNAs': '#3f51b5', # dark blue
    'ribosomal RNAs': '#f44336', # bright red
    'ORFs': '#00bcd4', # cyan
    'other genes': '#ce93d8' # light purple
}

def get_gene_info(feature):
    gene = feature.qualifiers.get('gene', [''])[0].lower()
    if not gene:
        gene = feature.qualifiers.get('product', [''])[0].lower()
        
    if gene.startswith('nad') or gene.startswith('ndh'):
        return COLOR_DICT['complex I (NADH dehydrogenase)'], 'complex I (NADH dehydrogenase)'
    elif gene.startswith('sdh'):
        return COLOR_DICT['complex II (succinate dehydrogenase)'], 'complex II (succinate dehydrogenase)'
    elif gene.startswith('cob'):
        return COLOR_DICT['complex III (ubichinol cytochrome c reductase)'], 'complex III (ubichinol cytochrome c reductase)'
    elif gene.startswith('cox'):
        return COLOR_DICT['complex IV (cytochrome c oxidase)'], 'complex IV (cytochrome c oxidase)'
    elif gene.startswith('atp'):
        return COLOR_DICT['ATP synthase'], 'ATP synthase'
    elif gene.startswith('ccm'):
        return COLOR_DICT['cytochrome c biogenesis'], 'cytochrome c biogenesis'
    elif gene.startswith('rpo'):
        return COLOR_DICT['RNA polymerase'], 'RNA polymerase'
    elif gene.startswith('rps'):
        return COLOR_DICT['ribosomal proteins (SSU)'], 'ribosomal proteins (SSU)'
    elif gene.startswith('rpl'):
        return COLOR_DICT['ribosomal proteins (LSU)'], 'ribosomal proteins (LSU)'
    elif gene.startswith('mat'):
        return COLOR_DICT['maturases'], 'maturases'
    elif gene.startswith('trn') or feature.type == 'tRNA':
        return COLOR_DICT['transfer RNAs'], 'transfer RNAs'
    elif gene.startswith('rrn') or feature.type == 'rRNA':
        return COLOR_DICT['ribosomal RNAs'], 'ribosomal RNAs'
    elif gene.startswith('orf'):
        return COLOR_DICT['ORFs'], 'ORFs'
    else:
        return COLOR_DICT['other genes'], 'other genes'

def repel_positions(positions, max_pos, delta, max_iter=3000):
    n = len(positions)
    if n <= 1: return positions
    new_pos = np.array(positions, dtype=float)
    for i in range(n):
        new_pos[i] += random.uniform(-0.01, 0.01) * delta
        new_pos[i] %= max_pos
    for _ in range(max_iter):
        moved = False
        order = np.argsort(new_pos)
        for i in range(n):
            idx1 = order[i]
            idx2 = order[(i + 1) % n]
            p1 = new_pos[idx1]
            p2 = new_pos[idx2]
            diff = (p2 - p1) % max_pos
            if diff < delta:
                overlap = delta - diff
                new_pos[idx1] = (p1 - overlap / 2.0) % max_pos
                new_pos[idx2] = (p2 + overlap / 2.0) % max_pos
                moved = True
        if not moved:
            break
    return new_pos

def draw_labels_with_lines(ax, sector, features, max_pos, r_start, r_end, r_text, is_forward):
    valid_features, names, positions = [], [], []
    for f in features:
        name = f.qualifiers.get("gene", [""])[0]
        if not name: continue
        valid_features.append(f)
        names.append(name)
        positions.append((int(f.location.start) + int(f.location.end)) / 2.0)
        
    if not positions: return
    
    # Use a much larger delta. Inner track needs more degrees because circumference is smaller.
    if is_forward:
        delta = max_pos * (6.0 / 360.0) 
    else:
        delta = max_pos * (9.0 / 360.0) 
        
    new_positions = repel_positions(positions, max_pos, delta, max_iter=8000)
    
    for name, orig_pos, new_pos in zip(names, positions, new_positions):
        theta_orig = sector.x_to_rad(orig_pos)
        theta_new = sector.x_to_rad(new_pos)
        
        # Phase unwrapping to prevent lines crossing the circle center
        diff_theta = theta_new - theta_orig
        if diff_theta > np.pi:
            theta_new -= 2 * np.pi
        elif diff_theta < -np.pi:
            theta_new += 2 * np.pi
            
        ax.plot([theta_orig, theta_new], [r_start, r_end], color="grey", lw=0.6, zorder=1)
        
        deg = np.degrees(theta_new)
        if is_forward:
            ha = "left"
            if 90 < deg < 270:
                rot = deg + 180
                ha = "right"
            else:
                rot = deg
        else:
            ha = "right"
            if 90 < deg < 270:
                rot = deg + 180
                ha = "left"
            else:
                rot = deg
                
        ax.text(theta_new, r_text, name, rotation=rot, ha=ha, va="center", 
                rotation_mode="anchor", fontsize=7.5, color="black", zorder=5)

def plot_circular_genome(gbk_file, chrom_num):
    gbk = Genbank(gbk_file)
    circos = Circos(sectors={gbk.name: gbk.range_size})
    sector = circos.get_sector(gbk.name)
    
    title_text = (
        r"$\mathit{Parashorea\ chinensis}$" + "\n"
        f"mitochondrion chromosome {chrom_num},\n"
        "complete sequence\n"
        f"{gbk.range_size:,} bp"
    )
    circos.text(title_text, r=0, size=15)
    
    features = gbk.extract_features("CDS") + gbk.extract_features("tRNA") + gbk.extract_features("rRNA")
    fwd_features = [f for f in features if f.location.strand == 1]
    rev_features = [f for f in features if f.location.strand == -1]
    
    track_line = sector.add_track((85, 85))
    track_line.axis(ec="black", lw=1.5)
    
    track_fwd = sector.add_track((85, 91))
    for f in fwd_features:
        color = get_gene_info(f)[0]
        track_fwd.genomic_features([f], plotstyle="box", fc=color, ec="black", lw=0.6, zorder=3)
    
    track_rev = sector.add_track((79, 85))
    for f in rev_features:
        color = get_gene_info(f)[0]
        track_rev.genomic_features([f], plotstyle="box", fc=color, ec="black", lw=0.6, zorder=3)
        
    gc_track = sector.add_track((50, 70))
    gc_track.axis(fc="#e0e0e0", ec="none", alpha=1.0) # Thick light grey circle
    
    pos_list, gc_contents = gbk.calc_gc_content(window_size=1000, step_size=500)
    avg_gc = sum(gc_contents) / len(gc_contents)
    gc_track.line([0, gbk.range_size], [avg_gc, avg_gc], color="silver", lw=1)
    gc_track.fill_between(pos_list, gc_contents, avg_gc, vmin=0, vmax=100, color="gray", alpha=0.9)
    
    sector.text("➞", x=gbk.range_size*0.01, r=95, size=24, color="darkgrey")
    sector.text("⟵", x=gbk.range_size*0.99, r=75, size=24, color="darkgrey")

    return circos, sector, fwd_features, rev_features, gbk.range_size

def main():
    os.makedirs("results", exist_ok=True)
    fig = plt.figure(figsize=(26, 13))
    
    ax1 = fig.add_subplot(121, polar=True)
    circos1, sector1, fwd1, rev1, max1 = plot_circular_genome("data/mt_VN1_contig1.gb", "1")
    circos1.plotfig(ax=ax1)
    draw_labels_with_lines(ax1, sector1, fwd1, max1, 91, 96, 97, True)
    draw_labels_with_lines(ax1, sector1, rev1, max1, 79, 74, 73, False)

    ax2 = fig.add_subplot(122, polar=True)
    circos2, sector2, fwd2, rev2, max2 = plot_circular_genome("data/mt_VN1_contig2.gb", "2")
    circos2.plotfig(ax=ax2)
    draw_labels_with_lines(ax2, sector2, fwd2, max2, 91, 96, 97, True)
    draw_labels_with_lines(ax2, sector2, rev2, max2, 79, 74, 73, False)

    legend_elements = [
        patches.Patch(facecolor=color, edgecolor="black", label=label)
        for label, color in COLOR_DICT.items()
    ]
    fig.legend(handles=legend_elements, loc="lower left", fontsize=11, bbox_to_anchor=(0.02, 0.05),
               title_fontsize=13, frameon=False, labelspacing=0.6, handleheight=1.2, handlelength=1.2)

    out_png = "results/Fig3_Circular_Genomic_Map.png"
    out_pdf = "results/Fig3_Circular_Genomic_Map.pdf"
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.savefig(out_pdf, bbox_inches='tight')
    print(f"Success! Saved outputs to {out_png}")

if __name__ == "__main__":
    main()
