import sys

def process_bed_file(input_file, output_file):
    with open(input_file, 'r') as f:
        raw_lines = [line.strip() for line in f if line.strip()]
    
    if not raw_lines:
        return

    outfile = open(output_file, 'w')
    
    data_start_idx = 0
    if raw_lines[0].lower().startswith('chrom'):
        outfile.write(raw_lines[0] + '\n')
        data_start_idx = 1
    
    data_lines = raw_lines[data_start_idx:]
    num_lines = len(data_lines)
    
    parsed_rows = []
    for line in data_lines:
        parsed_rows.append(line.split())

    for i in range(num_lines):
        chrom = parsed_rows[i][0]
        start = int(parsed_rows[i][1])
        end = int(parsed_rows[i][2])
        ancestry = parsed_rows[i][3]

        if end > start and i < num_lines - 1:
            end -= 1

        if start == end and i < num_lines - 1:
            next_start = int(parsed_rows[i+1][1])
            parsed_rows[i+1][1] = str(next_start + 1)

        outfile.write(f"{chrom}\t{start}\t{end}\t{ancestry}\n")
    
    outfile.close()

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 fix_bed.py input.bed output.bed")
    else:
        process_bed_file(sys.argv[1], sys.argv[2])
