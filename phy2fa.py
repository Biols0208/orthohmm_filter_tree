#!/dellfsqd2/ST_OCEAN/USER/lishuo11/09_test/zz_tmp/home_micromamba/envs/xgboost/bin/python3.9
import sys

def phy_to_fasta(input_phy_file, output_fasta_file):
    """将 phy 文件转换为 fasta 文件格式"""
    with open(input_phy_file, 'r') as infile, open(output_fasta_file, 'w') as outfile:
        first_line = True  # 用于标记是否是第一行
        for line in infile:
            line = line.strip()
            if not line:  # 跳过空行
                continue

            if first_line:
                first_line = False  # 跳过第一行
                continue

            # 按空格分割物种名称和序列
            parts = line.split()
            species_name = parts[0]
            sequence = ''.join(parts[1:])

            # 写入 fasta 格式，物种名称作为标识符
            outfile.write(f">{species_name}\n{sequence}\n")

if __name__ == "__main__":
    # 获取命令行参数
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_phy_file> <output_fasta_file>")
        sys.exit(1)
    
    input_phy_file = sys.argv[1]
    output_fasta_file = sys.argv[2]

    # 执行转换
    phy_to_fasta(input_phy_file, output_fasta_file)
    
