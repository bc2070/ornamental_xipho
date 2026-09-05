import pandas as pd
import glob
import os

def calculate_hwe_proportions(file_pattern):
    files = glob.glob(file_pattern)

    if not files:
        print("未找到匹配的文件")
        return

    output_filename = "hwe_proportions_summary.txt"

    with open(output_filename, 'w') as f_out:
        for file in files:
            try:
                df = pd.read_csv(
                    file,
                    sep='\s+',
                    usecols=['2,0', '1,1', '0,2'],
                    dtype={'2,0': float, '1,1': float, '0,2': float}
                )

                s20 = df['2,0'].sum()
                s11 = df['1,1'].sum()
                s02 = df['0,2'].sum()
                total = s20 + s11 + s02

                if total == 0:
                    continue

                A = (s20 + 0.5 * s11) / total

                a = 1 - A

                prop_20 = A ** 2
                prop_11 = 2 * A * a
                prop_02 = a ** 2

                filename = os.path.basename(file)

                result_line = f"{filename}\t{prop_20:.6f}\t{prop_11:.6f}\t{prop_02:.6f}\n"

                f_out.write(result_line)

            except Exception as e:
                print(f"处理文件 {file} 时出错: {e}")

    print(f"处理完成，结果已保存至: {output_filename}")

if __name__ == "__main__":
    calculate_hwe_proportions("*.posterior")
