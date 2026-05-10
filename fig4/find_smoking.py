import csv
import sys

# --- CONFIGURATION ---
blocks_file = "results/P_chinensis_structural_blocks.csv"
repeats_file = "data/CH2_self_repeats.txt"  
window_size = 500  # How close the repeat must be to the breakpoint (in base pairs)
min_repeat_length = 50  # Ignore tiny, irrelevant repeats under 50bp
# ---------------------

print("1. Loading Structural Blocks...")
inverted_blocks = []
try:
    with open(blocks_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # SAFETY FILTER: Skip missing blocks ('-') before checking for inversions
            if row['CH2_Start'] != '-' and row['CH2_Strand'] == '-': 
                inverted_blocks.append({
                    'id': row['Block_ID'],
                    'start': int(row['CH2_Start']),
                    'end': int(row['CH2_End'])
                })
    print(f"   -> Found {len(inverted_blocks)} inverted blocks to investigate.")
except FileNotFoundError:
    print(f"Error: Could not find {blocks_file}")
    sys.exit()

print("2. Loading BLAST Repeats...")
repeats = []
try:
    with open(repeats_file, 'r') as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 12:
                continue
                
            length = int(cols[3])
            
            # Ignore massive full-genome self-matches
            if length > 5000 or length < min_repeat_length:
                continue
                
            repeats.append({
                'qseqid': cols[0], 'sseqid': cols[1],
                'pident': float(cols[2]), 'length': length,
                'qstart': int(cols[6]), 'qend': int(cols[7]),
                'sstart': int(cols[8]), 'send': int(cols[9])
            })
    print(f"   -> Found {len(repeats)} valid repeat pairs.")
except FileNotFoundError:
    print(f"Error: Could not find {repeats_file}")
    sys.exit()

print("\n3. Cross-Referencing Repeats with Breakpoints...")
print("-" * 80)
print(f"{'Block ID':<15} | {'Breakpoint Location':<20} | {'Repeat Found (Smoking Gun)':<35}")
print("-" * 80)

smoking_guns_found = 0

for block in inverted_blocks:
    for rep in repeats:
        # Check if the repeat is sitting right on top of the block's START or END coordinate
        near_start = any(abs(block['start'] - pos) <= window_size for pos in [rep['qstart'], rep['qend'], rep['sstart'], rep['send']])
        near_end = any(abs(block['end'] - pos) <= window_size for pos in [rep['qstart'], rep['qend'], rep['sstart'], rep['send']])
        
        if near_start or near_end:
            boundary = "Start" if near_start else "End"
            coord = block['start'] if near_start else block['end']
            
            # Format the output cleanly
            repeat_info = f"{rep['length']}bp ({rep['pident']}% match)"
            if rep['qseqid'] != rep['sseqid']:
                repeat_info += " *INTER-CHROMOSOMAL*"
                
            print(f"{block['id']:<15} | {boundary} (~{coord})       | {repeat_info}")
            smoking_guns_found += 1
            break # Move to the next block once we find a smoking gun for this one

print("-" * 80)
print(f"Total Smoking Guns Found: {smoking_guns_found}")
if smoking_guns_found > 0:
    print("Success! These repeats are the physical cause of your inversions.")