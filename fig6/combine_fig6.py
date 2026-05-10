import matplotlib.pyplot as plt
import matplotlib.image as mpimg

def combine_figures():
    # Load the Circos output images
    try:
        img_a = mpimg.imread('results/Fig6A.png')
        img_b = mpimg.imread('results/Fig6B.png')
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
    ax_a.text(0.05, 0.95, 'A', transform=ax_a.transAxes, fontsize=30, fontweight='bold', va='top')
    ax_b.text(0.05, 0.95, 'B', transform=ax_b.transAxes, fontsize=30, fontweight='bold', va='top')
    
    import numpy as np

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
            
        ax.text(x, y, text, transform=ax.transAxes, fontsize=24, fontweight='bold',
                ha='center', va='center', rotation=rot)

    # Add exact mathematical labels for Fig 6A
    place_label(ax_a, 39, 'Cp DNA')
    place_label(ax_a, 278, 'Mt DNA')
    place_label(ax_a, 150, 'Mt DNA')

    # Add exact mathematical labels for Fig 6B
    place_label(ax_b, 82.8, 'Cp DNA')
    place_label(ax_b, 104.4, 'Mt DNA')
    place_label(ax_b, 277.2, 'Nc DNA')

    plt.tight_layout()
    
    output_file = 'results/Fig6_Combined.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Success! Combined figure saved to {output_file}")

if __name__ == "__main__":
    combine_figures()
