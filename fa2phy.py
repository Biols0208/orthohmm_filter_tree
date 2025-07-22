#!/usr/bin/env python3

"""
FASTA到PHYLIP格式转换工具 - 性能优化版
支持处理大型序列比对文件，Python 3版本
"""

import sys
import os
import argparse
from datetime import datetime


# 简单用法信息
USAGE = """
FASTA转PHYLIP格式转换工具

用法: fa2phy.py [input FASTA] [output PHYLIP]

示例: fa2phy.py input.fasta output.phy
"""

def msg(*args, **kwargs):
    """向stderr输出信息"""
    print(*args, file=sys.stderr, **kwargs)

def err(*args, **kwargs):
    """向stderr输出错误信息并退出"""
    msg("ERROR:", *args, **kwargs)
    sys.exit(1)

def parse_fasta_efficient(filename):
    """
    高效解析FASTA文件
    使用生成器避免将整个文件加载到内存
    """
    sequence_dict = {}
    sequence_list = []  # 保持序列顺序
    current_id = None
    
    try:
        with open(filename, 'r') as fh:
            seq_parts = []  # 收集序列片段
            
            for line in fh:
                line = line.rstrip()
                if not line:  # 跳过空行
                    continue
                    
                if line[0] == '>':
                    # 处理前一个序列
                    if current_id:
                        sequence_dict[current_id] = ''.join(seq_parts)
                    
                    # 提取新序列ID
                    header = line[1:].strip()
                    current_id = header.split()[0]
                    sequence_list.append(current_id)
                    seq_parts = []  # 重置序列片段列表
                else:
                    seq_parts.append(line)
            
            # 处理最后一个序列
            if current_id and seq_parts:
                sequence_dict[current_id] = ''.join(seq_parts)
        
        return sequence_dict, sequence_list
                
    except Exception as e:
        err(f"解析FASTA文件时出错: {e}")

def write_phylip(sequence_dict, sequence_list, outfile):
    """
    高效写入PHYLIP格式文件
    """
    # 检查比对长度
    alignment_length = 0
    for gene in sequence_dict:
        if alignment_length == 0:
            alignment_length = len(sequence_dict[gene])
        elif len(sequence_dict[gene]) != alignment_length:
            err(f"比对长度错误: {gene}序列长度({len(sequence_dict[gene])})与其他序列长度({alignment_length})不一致")
    
    # 找出最长的序列ID
    if sequence_list:
        longest_id_len = max(len(id) for id in sequence_list)
    else:
        err("没有找到有效的序列")
    
    # 写入PHYLIP文件
    try:
        with open(outfile, "w") as phyfile:
            # 写入序列数量和比对长度
            phyfile.write(f"{len(sequence_dict)} {alignment_length}\n")
            
            # 写入序列
            for gene in sequence_list:
                phyfile.write(f"{gene.ljust(longest_id_len)}   {sequence_dict[gene]}\n")
                
        msg(f"成功写入PHYLIP文件: {outfile}")
        msg(f"  序列数量: {len(sequence_dict)}")
        msg(f"  比对长度: {alignment_length}")
            
    except Exception as e:
        err(f"写入PHYLIP文件时出错: {e}")

def main():
    """主函数"""
    # 处理命令行参数
    if len(sys.argv) != 3:
        print(USAGE)
        sys.exit(0)
    
    fasta_file = sys.argv[1]
    phylip_file = sys.argv[2]
    
    # 验证输入文件
    if not os.path.isfile(fasta_file):
        err(f"找不到输入文件: {fasta_file}")
    
    if os.path.exists(phylip_file):
        msg(f"警告: 输出文件 {phylip_file} 已存在，将被覆盖")
    
    # 解析FASTA文件
    sequence_dict, sequence_list = parse_fasta_efficient(fasta_file)
    
    # 写入PHYLIP文件
    write_phylip(sequence_dict, sequence_list, phylip_file)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
