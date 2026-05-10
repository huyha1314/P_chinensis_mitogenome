import pysam
import matplotlib.pyplot as plt

# --- CONFIGURATION ---
bam_file = "data/breakpoint_169bp.bam"
contig = "contig_1"
start = 13000
end = 14000
breakpoint_start = 13522
breakpoint_end = 13689
# ---------------------

print(f"Reading BAM file: {bam_file}...")
bam = pysam.AlignmentFile(bam_file, "rb")

reads = []
for read in bam.fetch(contig, start, end):
    # Only keep long reads that map properly
    if not read.is_unmapped:
        reads.append((read.reference_start, read.reference_end, read.query_name))

bam.close()

# Sort reads by their start position for a clean waterfall look
reads.sort(key=lambda x: x[0])

print(f"Found {len(reads)} reads in this region. Drawing plot...")

# Set up the image
fig, ax = plt.subplots(figsize=(12, 6))

# Draw the breakpoint zone (The 169bp repeat) as a shaded red box
ax.axvspan(breakpoint_start, breakpoint_end, color='red', alpha=0.3, label="169bp Repeat Breakpoint")

# Draw each read as a horizontal line
y = 1
spanning_count = 0

# Define how many base pairs the read MUST extend past the repeat to count as proof
min_overhang = 50 

for r_start, r_end, name in reads:
    # STRICT SPANNING: Read extends deeply into unique DNA on both sides
    if (r_start <= breakpoint_start - min_overhang) and (r_end >= breakpoint_end + min_overhang):
        color = 'green'
        spanning_count += 1
    # WEAK SPANNING: Read crosses the boundary, but doesn't have enough overhang to be 100% sure
    elif (r_start < breakpoint_start) and (r_end > breakpoint_end):
        color = 'orange'
    # NO SPANNING: Read dies inside the repeat or doesn't touch it
    else:
        color = 'blue'
        
    ax.plot([r_start, r_end], [y, y], color=color, linewidth=2)
    y += 1

# Formatting - CẬP NHẬT KÍCH THƯỚC CHỮ TẠI ĐÂY
ax.set_title(f"ONT Long Reads Spanning Recombination Breakpoint ({contig})", fontsize=25) # Tăng size tiêu đề
ax.set_xlabel("Genomic Coordinate (bp)", fontsize=25) # Tăng size label trục X
ax.set_ylabel("Individual Sequenced Reads", fontsize=23) # Tăng size label trục Y
ax.tick_params(axis='x', labelsize=20) # Tăng size text (các con số) trên trục X
ax.set_yticks([]) # Hide Y axis numbers
ax.legend(loc="upper right", fontsize=16) # Tăng size text legend

plt.tight_layout()

# Save as high-resolution PNG
output_png = "results/FigureS1.png"
plt.savefig(output_png, dpi=600, bbox_inches='tight')

# Save as Vector PDF
output_pdf = "results/FigureS1.pdf"
plt.savefig(output_pdf, format="pdf", bbox_inches='tight')

print(f"\nSUCCESS! {spanning_count} reads completely spanned the recombination site.")
print(f"Images saved as {output_png} and {output_pdf}.")