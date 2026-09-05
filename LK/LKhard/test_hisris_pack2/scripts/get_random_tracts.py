import sys
import random

def read_chrom_file(chrom_file):
    chrom_lengths = []
    with open(chrom_file, 'r') as f:
        for line in f:
            if line.lower().startswith('chrom') or not line.strip():
                continue
            parts = line.strip().split()
            chrom_lengths.append((parts[0], int(parts[2]))) 
    return chrom_lengths

def read_fragment_file(frag_file):
    fragments = []
    with open(frag_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 5: continue 
            ancestry = parts[3]
            length = int(parts[4])
            fragments.append((length, ancestry))
    return fragments

def generate_random_fragments(chrom_lengths, fragments):
    results = []
    for length, ancestry in fragments:
        eligible_chroms = [c for c in chrom_lengths if c[1] >= length]
        
        if not eligible_chroms:
            continue

        chrom, chrom_max = random.choice(eligible_chroms)
        random_start = random.randint(0, chrom_max - length)
        random_end = random_start + length
        results.append(f"{chrom}\t{random_start}\t{random_end}\t{ancestry}")
    return results

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 get_random_tracts.py chrom.txt fragments.bed")
    else:
        chroms = read_chrom_file(sys.argv[1])
        frags = read_fragment_file(sys.argv[2])
        random_res = generate_random_fragments(chroms, frags)
        for line in random_res:
            print(line)
