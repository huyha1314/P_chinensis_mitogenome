import csv
from Bio import SeqIO

# 1. Define your files
gb_file = "data/VN1.gb"
csv_file = "results/P_chinensis_structural_blocks.csv"
output_file = "results/VN1_Gene_Block_Map.csv"

print(f"Loading LCB blocks from {csv_file}...")

# 2. Load the structural blocks for VN1
vn1_blocks = []
with open(csv_file, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        # Only process blocks that actually exist in VN1
        if row['VN1_Start'] != '-':
            vn1_blocks.append({
                'Block_ID': row['Block_ID'],
                'Start': int(row['VN1_Start']),
                'End': int(row['VN1_End']),
                'Strand': row['VN1_Strand']
            })

print(f"Loaded {len(vn1_blocks)} blocks. Mapping genes from {gb_file}...")

# 3. Parse the GenBank file and map genes
mapped_genes = []
try:
    record = SeqIO.read(gb_file, "genbank")
    
    for feature in record.features:
        # We only care about actual genes and RNAs, not source or generic regions
        if feature.type in ["CDS", "tRNA", "rRNA"]:
            # Extract the best possible name for the gene
            gene_name = feature.qualifiers.get('gene', [''])[0]
            if not gene_name:
                gene_name = feature.qualifiers.get('product', ['Unidentified'])[0]
            
            # Biopython uses 0-based indexing for start, so we add 1 to match Mauve's 1-based CSV
            gene_start = int(feature.location.start) + 1
            gene_end = int(feature.location.end)
            gene_strand = '+' if feature.location.strand == 1 else '-'
            
            # Find which block this gene falls into
            assigned_block = "Intergenic/Gap"
            block_strand = ""
            for block in vn1_blocks:
                # Check if the gene falls within the block's boundaries
                if gene_start >= block['Start'] and gene_end <= block['End']:
                    assigned_block = block['Block_ID']
                    block_strand = block['Strand']
                    break
            
            # Avoid duplicating the same gene if it has multiple CDS segments
            gene_entry = [gene_name, feature.type, gene_start, gene_end, gene_strand, assigned_block, block_strand]
            if gene_entry not in mapped_genes:
                mapped_genes.append(gene_entry)

    # 4. Write the results to a new CSV
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Gene_Name', 'Feature_Type', 'Gene_Start', 'Gene_End', 'Gene_Strand', 'VN1_Block_ID', 'Block_Orientation'])
        writer.writerows(mapped_genes)
        
    print(f"Success! Mapped {len(mapped_genes)} features.")
    print(f"Saved to {output_file}")

except Exception as e:
    print(f"An error occurred: {e}")