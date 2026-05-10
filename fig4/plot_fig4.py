import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

# 1. Load the data
csv_file = "results/P_chinensis_structural_blocks.csv"
df = pd.read_csv(csv_file)

# 2. Setup the Plot
genomes = ['VN1', 'CH1', 'CH2', 'Shorea', 'Thebroma']

# Mapping to long names for display (no underscores, Italics)
display_map = {
    'VN1': 'VN1 Parashorea chinensis',
    'CH1': 'CH1 Parashorea chinensis',
    'CH2': 'CH2 Parashorea chinensis',
    'Shorea': 'Shorea roxburghii',
    'Thebroma': 'Thebroma cacao (outgroup)'
}

fig, ax = plt.subplots(figsize=(16, 11))

import matplotlib.colors as mcolors

# Colors for blocks (Gradient across the genome)
colors = plt.cm.viridis(np.linspace(0, 1, len(df)))

# Colors for ribbons (Direct vs Inverted)
color_direct = '#800000' # Maroon
color_inverted = '#1f4e78' # Dark blue


# 3. Draw Genomes and Links
for i in range(len(genomes)):
    y_pos = len(genomes) - i
    genome = genomes[i]
    
    # --- DRAW GENOME LINE ---
    end_col = f'{genome}_End'
    max_val = pd.to_numeric(df[end_col], errors='coerce').max()
    ax.hlines(y_pos, 0, max_val, color='grey', linewidth=0.8, zorder=1)
    
    # Label with Bold Italic
    ax.text(-8000, y_pos, display_map[genome], va='center', ha='right', 
            fontweight='bold', style='italic', fontsize=12)

    # --- DRAW BLOCKS (Centered on the line, no offset) ---
    for idx, row in df.iterrows():
        start = row[f'{genome}_Start']
        end = row[f'{genome}_End']
        if start == '-': continue
        
        start, end = float(start), float(end)
        width = end - start
        
        # Check if this block has a connection to adjacent genomes
        has_up = (i > 0) and (row[f'{genomes[i-1]}_Start'] != '-')
        has_down = (i < len(genomes) - 1) and (row[f'{genomes[i+1]}_Start'] != '-')
        
        if has_up or has_down:
            block_color = colors[idx]
        else:
            block_color = 'white'
        
        # Draw the block rectangle
        rect = patches.Rectangle((start, y_pos - 0.05), width, 0.1, 
                                 edgecolor='none', facecolor=block_color, alpha=0.9, zorder=4)
        ax.add_patch(rect)

    # --- DRAW LINKS (Filled ribbons with solid borders) ---
    if i < len(genomes) - 1:
        next_genome = genomes[i+1]
        next_y = y_pos - 1
        
        for idx, row in df.iterrows():
            s1, e1, st1 = row[f'{genome}_Start'], row[f'{genome}_End'], row[f'{genome}_Strand']
            s2, e2, st2 = row[f'{next_genome}_Start'], row[f'{next_genome}_End'], row[f'{next_genome}_Strand']
            
            if s1 == '-' or s2 == '-': continue
            
            s1, e1, s2, e2 = float(s1), float(e1), float(s2), float(e2)
            
            if st1 == st2:
                # Direct links
                poly_coords = [[s1, y_pos - 0.05], [e1, y_pos - 0.05], [e2, next_y + 0.05], [s2, next_y + 0.05]]
                ribbon_color = color_direct
            else:
                # Inverted links (crossed)
                poly_coords = [[s1, y_pos - 0.05], [e1, y_pos - 0.05], [s2, next_y + 0.05], [e2, next_y + 0.05]]
                ribbon_color = color_inverted
                
            f_color = mcolors.to_rgba(ribbon_color, alpha=0.35)
            e_color = mcolors.to_rgba(ribbon_color, alpha=0.9)
                
            polygon = patches.Polygon(poly_coords, closed=True, 
                                      facecolor=f_color, edgecolor=e_color, 
                                      linewidth=0.5, zorder=2)
            ax.add_patch(polygon)

# 4. Final Touches
ax.set_ylim(0.5, len(genomes) + 0.5)
all_end_cols = [f'{g}_End' for g in genomes]
max_coord = df[all_end_cols].apply(pd.to_numeric, errors='coerce').max().max()
ax.set_xlim(-160000, max_coord * 1.05)

# Formatting axes
ax.set_yticks([])
ax.tick_params(axis='x', labelsize=14) # Increased X-axis size
ax.set_xlabel("Genomic Position (bp)", fontsize=15, labelpad=10)
ax.set_title("Parashorea chinensis Mitogenome Synteny (5-way comparison)", fontsize=18, fontweight='bold', pad=30)

# Legend
direct_patch = patches.Patch(facecolor=mcolors.to_rgba(color_direct, 0.4), edgecolor=color_direct, label='Direct Synteny (Same Strand)')
inverted_patch = patches.Patch(facecolor=mcolors.to_rgba(color_inverted, 0.4), edgecolor=color_inverted, label='Inverted Synteny (Crossed)')
ax.legend(handles=[direct_patch, inverted_patch], loc='lower right', fontsize=12)

plt.tight_layout()

# 5. Save results
plt.savefig("results/Fig4_Synteny_Map.png", dpi=300)
plt.savefig("results/Fig4_Synteny_Map.pdf")

print("Success! Updated Synteny Plot with separated links in results/Fig4_Synteny_Map.png")
