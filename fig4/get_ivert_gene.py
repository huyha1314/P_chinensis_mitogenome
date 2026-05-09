import csv

blocks_file = "results/P_chinensis_structural_blocks.csv"
genes_file = "results/VN1_Gene_Block_Map.csv"

# 1. Find all Block IDs that are Inverted (-) in the CH2 genome
ch2_inverted_blocks = set()
with open(blocks_file, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        # Check if the CH2 strand is '-'
        if row['CH2_Strand'] == '-':
            ch2_inverted_blocks.add(row['Block_ID'])

# 2. Match those Block IDs to the mapped genes
print(f"Found {len(ch2_inverted_blocks)} inverted structural blocks in CH2.\n")
print(f"{'Gene Name':<15} | {'Type':<10} | {'Located in Block':<15}")
print("-" * 45)

affected_genes = 0
with open(genes_file, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        # If the gene's block is in our list of inverted blocks, print it!
        if row['VN1_Block_ID'] in ch2_inverted_blocks:
            print(f"{row['Gene_Name']:<15} | {row['Feature_Type']:<10} | {row['VN1_Block_ID']:<15}")
            affected_genes += 1

print("-" * 45)
print(f"Total genes trapped inside the CH2 inversions: {affected_genes}")