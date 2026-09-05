import sys
import os
import math

def load_hwe_proportions(prop_file):
    """读取比例文件，建立 ID -> [p20, p11, p02] 映射"""
    hwe_dict = {}
    if not os.path.exists(prop_file):
        return hwe_dict
    with open(prop_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if not parts: continue
            sample_id = parts[0].split('.')[0].split('_')[0]
            try:
                probs = [float(x) for x in parts[1:]]
                hwe_dict[sample_id] = probs
            except ValueError:
                continue
    return hwe_dict

def get_sample_id_from_filename(filename):
    """提取ID: APS1_ahmm2bed_1 -> APS1"""
    return os.path.basename(filename).split('_')[0]

def process_chromosome(input_files, prop_file, output_file):
    hwe_lookup = load_hwe_proportions(prop_file)
    state_map = {'2,0': 0, '1,1': 1, '0,2': 2}
    
    all_breakpoints = set()
    samples_data = []
    chrom_name = ""

    for f_path in input_files:
        if not os.path.exists(f_path): continue
        sid = get_sample_id_from_filename(f_path)
        sprobs = hwe_lookup.get(sid, [1.0, 1.0, 1.0])
        
        segments = []
        with open(f_path, 'r') as f:
            for line in f:
                p = line.strip().split('\t')
                if len(p) < 4: continue
                if not chrom_name: chrom_name = p[0]
                start, end, state = int(p[1]), int(p[2]), p[3]
                segments.append((start, end, state))
                all_breakpoints.add(start)
                all_breakpoints.add(end)
        
        samples_data.append({
            'id': sid,
            'probs': sprobs,
            'segments': segments,
            'ptr': 0
        })

    if not samples_data: return

    sorted_breakpoints = sorted(list(all_breakpoints))

    with open(output_file, 'w') as out:
        for i in range(len(sorted_breakpoints) - 1):
            curr_start = sorted_breakpoints[i]
            curr_end = sorted_breakpoints[i+1]
            
            total_prob = 1.0
            
            for s_info in samples_data:
                state = "NA"
                while s_info['ptr'] < len(s_info['segments']):
                    seg_s, seg_e, seg_st = s_info['segments'][s_info['ptr']]
                    if seg_e <= curr_start:
                        s_info['ptr'] += 1
                        continue
                    if seg_s <= curr_start and seg_e >= curr_end:
                        state = seg_st
                    break
                
                if state in state_map:
                    idx = state_map[state]
                    total_prob *= s_info['probs'][idx]

            if total_prob > 0:
                neg_log10 = -math.log10(total_prob)
            else:
                neg_log10 = 999.0 
            
            # 输出格式：chrom, start, end, total_prob, neg_log10
            out.write(f"{chrom_name}\t{curr_start}\t{curr_end}\t{total_prob:.12f}\t{neg_log10:.6f}\n")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.exit(1)
    
    prop_f = sys.argv[1]
    out_f = sys.argv[2]
    in_fs = sys.argv[3:]
    
    process_chromosome(in_fs, prop_f, out_f)
