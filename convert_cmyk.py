#!/usr/bin/env python3
import os
import glob
from PIL import Image

# Output directory/file mapping
FIGURE_MAPPING = {
    "fig1/results/Fig1_Parashorea_Geographic_Map_Final.png": "fig1/results/Fig1_Parashorea_Geographic_Map_Final_CMYK.tiff",
    "fig2/results/Fig2.png": "fig2/results/Fig2_CMYK.tiff",
    "fig3/results/Fig3_Circular_Genomic_Map.png": "fig3/results/Fig3_Circular_Genomic_Map_CMYK.tiff",
    "fig4/results/Figure_4_Synteny.png": "fig4/results/Figure_4_Synteny_CMYK.tiff",
    "fig5/results/Figure_5_Combined.png": "fig5/results/Figure_5_Combined_CMYK.tiff",
    "fig6/results/Figure_6.png": "fig6/results/Figure_6_CMYK.tiff",
    "fig7/results/Fig7_Combined.png": "fig7/results/Fig7_Combined_CMYK.tiff",
    "fig8/results/Figure_8.png": "fig8/results/Figure_8_CMYK.tiff",
    "figS1/results/FigureS1.png": "figS1/results/FigureS1_CMYK.tiff"
}

def convert_to_cmyk_tiff(src, dest):
    if not os.path.exists(src):
        print(f"Warning: Source file {src} does not exist.")
        return False
    
    print(f"Converting {src} to {dest}...")
    try:
        im = Image.open(src)
        
        # Handle transparency/alpha channels (CMYK does not support alpha)
        if im.mode in ("RGBA", "LA") or (im.mode == "P" and "transparency" in im.info):
            alpha = im.convert("RGBA").split()[-1]
            bg = Image.new("RGBA", im.size, (255, 255, 255, 255))
            bg.paste(im, mask=alpha)
            im = bg.convert("RGB")
        elif im.mode != "RGB":
            im = im.convert("RGB")
            
        # Standard publication dimensions
        # Full text width: 17.0 cm = 6.6929 inches
        # Max height: 22.5 cm = 8.858 inches
        target_width_in = 6.6929
        max_height_in = 8.858
        
        aspect_ratio = im.height / im.width
        target_height_in = target_width_in * aspect_ratio
        
        if target_height_in > max_height_in:
            target_height_in = max_height_in
            target_width_in = target_height_in / aspect_ratio
            
        # Target resolution: 600 DPI for publication quality line art/figures
        dpi = 600
        target_width_px = int(target_width_in * dpi)
        target_height_px = int(target_height_in * dpi)
        
        # Resize using high-quality LANCZOS filter
        im_resized = im.resize((target_width_px, target_height_px), Image.Resampling.LANCZOS)
        
        # Convert to CMYK
        cmyk_im = im_resized.convert("CMYK")
        
        # Save as TIFF with LZW compression
        cmyk_im.save(dest, "TIFF", dpi=(dpi, dpi), compression="tiff_lzw")
        print(f"Successfully saved {dest} (Size: {target_width_px}x{target_height_px} px, {target_width_in:.2f}x{target_height_in:.2f} in, Mode: CMYK)")
        return True
    except Exception as e:
        print(f"Error converting {src}: {e}")
        return False

if __name__ == "__main__":
    import sys
    if len(sys.argv) == 3:
        src = sys.argv[1]
        dest = sys.argv[2]
        convert_to_cmyk_tiff(src, dest)
    else:
        success_count = 0
        for src, dest in FIGURE_MAPPING.items():
            if convert_to_cmyk_tiff(src, dest):
                success_count += 1
        print(f"\nCompleted: {success_count}/{len(FIGURE_MAPPING)} figures converted successfully.")
