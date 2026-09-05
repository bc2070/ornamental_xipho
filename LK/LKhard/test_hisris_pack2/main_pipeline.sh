#!/usr/bin/bash
#SBATCH --job-name=Master_Check
#SBATCH -p long
#SBATCH -c 8
#SBATCH --mem=40G
#SBATCH -o logs/%j.out
#SBATCH -e logs/%j.err

set -e

module load bedtools
module load R

source config.ini

mkdir -p logs result_final final_rand_probs scripts temp_splits
mkdir -p "$PREPARED_POSTERIOR_DIR"

ROOT_PATH=$(pwd)
export SCRIPT_PATH="$ROOT_PATH/$SCRIPT_DIR"

check_status() {
    local file=$1
    local msg=$2

    if [ ! -f "$file" ] || [ ! -s "$file" ]; then
        echo "ERROR: $msg (File $file is missing or empty)"
        exit 1
    else
        local count=$(wc -l < "$file")
        echo "CHECK: $file created with $count lines."
    fi
}


# ==========================================================
# Step 0: Posterior preprocessing
# ==========================================================

echo "=========================================================="
echo "==== [$(date)] Step 0: Posterior Preprocessing"
echo "=========================================================="

rm -rf "$PREPARED_POSTERIOR_DIR"
mkdir -p "$PREPARED_POSTERIOR_DIR"

# ----------------------------------------------------------
# 清理当前目录中旧的 posterior 和 ahmm2bed
# ----------------------------------------------------------

rm -f ./*.posterior
rm -f ./*_ahmm2bed

# ----------------------------------------------------------
# Process every posterior listed in file_list.txt
# ----------------------------------------------------------

while read -r input_file; do

    # 跳过空行
    [ -z "$input_file" ] && continue

    # 检查原始 posterior
    if [ ! -f "$input_file" ]; then
        echo "ERROR: Input posterior not found:"
        echo "       $input_file"
        exit 1
    fi

    fname=$(basename "$input_file")

    output_file="$PREPARED_POSTERIOR_DIR/$fname"

    ranking_file="$PREPARED_POSTERIOR_DIR/${fname}.ranking.txt"

    echo ""
    echo "----------------------------------------------------------"
    echo "Input : $input_file"
    echo "Output: $output_file"
    echo "----------------------------------------------------------"

    python3 "$SCRIPT_PATH/00.prepare_posterior.py" \
        "$input_file" \
        "$output_file" \
        "$ranking_file"

    # 检查输出
    if [ ! -s "$output_file" ]; then
        echo "ERROR: Prepared posterior is missing:"
        echo "       $output_file"
        exit 1
    fi

    # 建立链接
    ln -sf \
        "$output_file" \
        "./$fname"

    echo "Prepared posterior linked: $fname"

done < "$POSTERIOR_LIST"


echo ""
echo "=========================================================="
echo "Prepared posterior files"
echo "=========================================================="

ls -lh ./*.posterior


echo ""
echo "=========================================================="
echo "Ancestry ranking files"
echo "=========================================================="

ls -lh "$PREPARED_POSTERIOR_DIR"/*.ranking.txt

# ==========================================================
# Step 1: Original two-ancestry pipeline
# ==========================================================

echo "=========================================================="
echo "==== [$(date)] Step 1: Baseline Pipeline"
echo "=========================================================="

python3 "$SCRIPT_PATH/ahmm_bed.2.py" \
    "./" \
    "./" \
    "$CHROM_FILE" \
    "$THRESHOLD_P1"

check_status \
    "hwe_proportions_summary.txt" \
    "Step 1 HWE summary"

HWE_SUMMARY=$(readlink -f "hwe_proportions_summary.txt")

CHROMS=$(awk 'NR>1 {print $1}' "$CHROM_FILE" | sort -uV)


echo "Processing baseline probabilities chromosome by chromosome..."

[ -f "$FINAL_PROB_BED" ] && rm "$FINAL_PROB_BED"


for chr in $CHROMS; do

    rm -f temp_splits/*.${chr}

    for f in *_ahmm2bed; do

        awk -v c="$chr" \
            'NR>1 && $1==c {print $0}' \
            "$f" \
            > "temp_splits/${f}.${chr}" || true

    done

    if [ "$(ls temp_splits/*_ahmm2bed.${chr} 2>/dev/null | wc -l)" -gt 0 ]; then

        python3 "$SCRIPT_PATH/05.probility.py" \
            "$HWE_SUMMARY" \
            "temp_splits/${chr}.prob.tmp" \
            temp_splits/*_ahmm2bed."${chr}"

        cat "temp_splits/${chr}.prob.tmp" >> "$FINAL_PROB_BED"

        echo "Chromosome $chr: Done."

    else

        echo "Chromosome $chr: No data found, skipping."

    fi

done


check_status \
    "$FINAL_PROB_BED" \
    "Baseline probability calculation"


echo "==== [$(date)] Step 2: Randomization ===="
cat << 'EOF' > scripts/worker_rand.sh
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
EOF
chmod +x scripts/worker_rand.sh

ls *_ahmm2bed | xargs -I {} -P $THREADS ./scripts/worker_rand.sh {} "$CHROM_FILE" "$ITERATIONS" "result_final" "$SCRIPT_PATH"
echo "Randomization files generated: $(ls result_final/*.fixed.bed | wc -l) files."

echo "==== [$(date)] Step 2.2: Random Probabilities & Threshold ===="
rm -rf final_rand_probs/*.prob.bed

for chr in $CHROMS; do
    mkdir -p "temp_splits/$chr"
    for f in result_final/*.fixed.bed; do
        fname=$(basename "$f")
        awk -v c="$chr" '$1==c {print $0}' "$f" > "temp_splits/$chr/${fname}.${chr}"
    done
    
    for i in $(seq 1 $ITERATIONS); do
        iter_files=( temp_splits/"$chr"/*.rand."$i".*.${chr} )
        if [ -f "${iter_files[0]}" ]; then
            python3 "$SCRIPT_PATH/05.1.probility.py" "$HWE_SUMMARY" "temp_splits/rand.${i}.${chr}.tmp" "${iter_files[@]}"
            cat "temp_splits/rand.${i}.${chr}.tmp" >> "final_rand_probs/iter_${i}.prob.bed"
        fi
    done
    rm -rf "temp_splits/$chr"
    echo "Random prob calculation for $chr: Done."
done

echo "Computing 97.5% threshold..."
# 检查输入给 Python 的数据是否存在
PROB_DATA_COUNT=$(find final_rand_probs -name "iter_*.prob.bed" | xargs cat | awk '{print $5}' | wc -l)
echo "Total probability data points for threshold: $PROB_DATA_COUNT"

if [ "$PROB_DATA_COUNT" -eq 0 ]; then
    echo "ERROR: No probability data collected for threshold calculation!"
    exit 1
fi

find final_rand_probs -name "iter_*.prob.bed" | xargs cat | awk '{if($5!="") print $5}' | \
python3 -c "
import sys, math
try:
    data = [float(l.strip()) for l in sys.stdin if l.strip()]
    data.sort(); n=len(data)
    idx=math.ceil(n*0.975)-1
    t=data[max(0,min(idx,n-1))]
    print(f'Threshold calculated: {t}')
    with open('$THRESHOLD_RESULT','w') as f: f.write(f'Threshold_97.5: {t}\n')
except Exception as e:
    print(f'Python Threshold Error: {e}')
    sys.exit(1)
"
check_status "$THRESHOLD_RESULT" "Threshold calculation"

echo "==== [$(date)] Step 3: Plotting ===="
PLOT_INPUT="plot_ready.bed"
# 增加一个检查：确保 $4 是数字
awk 'BEGIN{OFS="\t"} {
    if($4 ~ /^[0-9.]+$/) {
        lp=($4<=0)?20:-log($4)/log(10); 
        print $1,$2,$3,$4,lp
    }
}' "$FINAL_PROB_BED" > "$PLOT_INPUT"

check_status "$PLOT_INPUT" "Plotting input data"
Rscript "$SCRIPT_PATH/plot_manhattan.R" "$PLOT_INPUT" "$THRESHOLD_RESULT" "$PLOT_PREFIX"

echo "==== [$(date)] All Steps Completed! ===="

