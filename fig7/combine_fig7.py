import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.patches as mpatches
import numpy as np

def combine_figures():
    # Load the Circos output images
    try:
        img_a = mpimg.imread('results/Fig7A.png')
        img_b = mpimg.imread('results/Fig7B.png')
    except Exception as e:
        print(f"Error loading images: {e}")
        return

    # Create a figure with two subplots side-by-side
    fig, axes = plt.subplots(1, 2, figsize=(20, 10))
    
    # Plot Fig6A on the left
    ax_a = axes[0]
    ax_a.imshow(img_a)
    ax_a.axis('off')
    
    # Plot Fig6B on the right
    ax_b = axes[1]
    ax_b.imshow(img_b)
    ax_b.axis('off')

    # Add 'A' and 'B' labels
    ax_a.text(0.05, 0.98, 'A', transform=ax_a.transAxes, fontsize=35, va='top')
    ax_b.text(0.05, 0.98, 'B', transform=ax_b.transAxes, fontsize=35,  va='top')
    
    def place_label(ax, angle_deg, text, r=0.53):
        # Image is centered at (0.5, 0.5)
        theta = np.radians(angle_deg)
        x = 0.5 + r * np.cos(theta)
        y = 0.5 + r * np.sin(theta)
        
        # Tangential rotation
        rot = angle_deg - 90
        if rot < -90:
            rot += 180
        elif rot > 90:
            rot -= 180
            
        ax.text(x, y, text, transform=ax.transAxes, fontsize=25, fontweight='bold',
                ha='center', va='center', rotation=rot)

    # --- FIG 7A LABELS ---
    # Kept r=0.52 to keep them close but safe
    place_label(ax_a, 39, 'Cp DNA', r=0.52)
    place_label(ax_a, 278, 'Mt DNA', r=0.52)
    place_label(ax_a, 150, 'Mt DNA', r=0.52)

    # --- FIG 7B LABELS ---
    # Increased r=0.58 to push labels outward and prevent overlap with organelle names.
    # Adjusted angles slightly to spread them perfectly over the red/green bands.
    place_label(ax_b, 75, 'Cp DNA', r=0.54)     
    place_label(ax_b, 115, 'Mt DNA', r=0.54)    
    place_label(ax_b, 270, 'Nc DNA', r=0.54)

    # --- ADD CUSTOM LEGEND ---
    # Create colored boxes matching the Circos ideograms
    mt_patch = mpatches.Patch(color='tomato', label='Mt DNA')
    cp_patch = mpatches.Patch(color='mediumseagreen', label='Cp DNA')
    nc_patch = mpatches.Patch(color='black', label='Nc DNA')

    # Add the legend to the bottom center of the whole figure
    fig.legend(handles=[mt_patch, cp_patch, nc_patch],
               loc='lower center', 
               ncol=3,               # Arrange horizontally in 3 columns
               fontsize=25, 
               bbox_to_anchor=(0.5, 0.02), # Pin it near the bottom edge
               frameon=False)        # Remove the box border for a cleaner look

    # Apply tight layout, but leave room at the bottom for the new legend
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15) 
    
    # Base output path
    output_prefix = 'results/Fig7_Combined'
    
    # Export to PNG, TIFF, and PDF
    formats = ['png', 'tiff', 'pdf']
    for ext in formats:
        output_file = f"{output_prefix}.{ext}"
        if ext == 'tiff':
            plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white', format=ext, pil_kwargs={"compression": "tiff_lzw"})
        else:
            plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white', format=ext)
        print(f"Success! Combined figure saved to {output_file}")

if __name__ == "__main__":
    combine_figures()