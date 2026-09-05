#!/bin/bash
input_bed=$1; CHROM_FILE=$2; ITERATIONS=$3; RESULT_DIR=$4; SCRIPT_PATH=$5
sName=$(basename "$input_bed" _ahmm2bed)
t_dir="temp_${sName}"; mkdir -p "$t_dir"

python3 "$SCRIPT_PATH/fix_bed.py" "$input_bed" "$t_dir/fixed.bed"
awk 'NR>1 && $4!="1,1" && $4!="0,2" {print $1"\t"$2"\t"$3"\t"$4}' "$t_dir/fixed.bed" > "$t_dir/bg.bed"
awk 'NR>1 && $4=="1,1" {print $1"\t"$2"\t"$3"\t1,1\t"($3-$2)}' "$t_dir/fixed.bed" > "$t_dir/hetro.bed"
awk 'NR>1 && $4=="0,2" {print $1"\t"$2"\t"$3"\t0,2\t"($3-$2)}' "$t_dir/fixed.bed" > "$t_dir/minor.bed"

for i in $(seq 1 $ITERATIONS); do
    python3 "$SCRIPT_PATH/get_random_tracts.py" "$CHROM_FILE" "$t_dir/hetro.bed" > "$t_dir/r_h.bed"
    python3 "$SCRIPT_PATH/get_random_tracts.py" "$CHROM_FILE" "$t_dir/minor.bed" > "$t_dir/r_m.bed"
    cat "$t_dir/bg.bed" "$t_dir/r_h.bed" "$t_dir/r_m.bed" | sort -k1,1 -k2,2n | \
    bedtools merge -i stdin -c 4 -o distinct | \
    python3 "$SCRIPT_PATH/fix_random_results.py" /dev/stdin "${RESULT_DIR}/${sName}.rand.$i.fixed.bed"
done
rm -rf "$t_dir"
