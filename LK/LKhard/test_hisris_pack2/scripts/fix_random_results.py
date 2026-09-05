import sys
import random

def select_best_ancestry(ancestry_str):
    parts = [p.strip() for p in ancestry_str.split(',')]
    if any(x in ['NA', 'N', '', 'null'] for x in parts):
        return "NA"
    
    pairs = []
    for i in range(0, len(parts), 2):
        if i + 1 < len(parts):
            pairs.append(f"{parts[i]},{parts[i+1]}")
    
    if not pairs:
        return "NA"

    unique_pairs = set(pairs)
    
    priority_states = {"0,2", "1,1"}
    
    available_priority = unique_pairs.intersection(priority_states)
    
    if available_priority:
        return random.choice(list(available_priority))
    elif "2,0" in unique_pairs:
        return "2,0"
    else:
        return "NA"

def process_file(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            line = line.strip()
            if not line:
                continue
            
            columns = line.split('\t')
            if len(columns) < 4:
                continue

            chrom, start, end, ancestry = columns[0], columns[1], columns[2], columns[3]
            
            new_ancestry = select_best_ancestry(ancestry)

            outfile.write(f"{chrom}\t{start}\t{end}\t{new_ancestry}\n")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python3 02.do.rand.result.py <input> <output>")
        sys.exit(1)
    process_file(sys.argv[1], sys.argv[2])
