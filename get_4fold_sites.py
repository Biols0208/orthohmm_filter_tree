#!/usr/bin/env python3
"""
从PHYLIP或FASTA格式比对文件中提取4重简并位点，处理缺口（gaps）
用法: python extract_4d_sites_universal.py input_file output_file [--input-format FORMAT] [--output-format FORMAT] [--include-gaps] [--mask-gaps]

FORMAT选项可以是: phylip, fasta (自动检测)
"""

import sys
import re
import os
import argparse
from collections import OrderedDict

class SequenceParser:
    """序列文件解析的基类"""
    
    def __init__(self, filename):
        self.filename = filename
        self.sequences = OrderedDict()  # 使用有序字典保持序列顺序
        self.seq_length = 0
    
    def get_sequences(self):
        """返回解析的序列"""
        seq_order = list(self.sequences.keys())
        return self.sequences, seq_order, self.seq_length

class PhylipParser(SequenceParser):
    """解析PHYLIP格式文件的类"""
    
    def __init__(self, filename):
        super().__init__(filename)
        self.parse()
    
    def parse(self):
        """解析PHYLIP文件"""
        try:
            with open(self.filename, 'r') as f:
                lines = [l.strip() for l in f.readlines() if l.strip()]
                
            # 第一行包含序列数量和长度
            if not lines:
                raise ValueError("文件为空")
                
            header = re.split(r'\s+', lines[0].strip())
            if len(header) < 2:
                raise ValueError("PHYLIP文件头格式错误")
            
            num_seqs = int(header[0])
            self.seq_length = int(header[1])
            
            # 确定文件格式（交错式或序列式）
            if len(lines) >= num_seqs + 1:
                # 检查是否为序列式格式
                first_seq_len = len(''.join(re.split(r'\s+', lines[1], 1)[1:]))
                if first_seq_len == self.seq_length:
                    self.parse_sequential(lines, num_seqs)
                else:
                    # 假设为交错式格式
                    self.parse_interleaved(lines, num_seqs)
            else:
                raise ValueError("PHYLIP文件格式无效或序列数量不匹配")
                
        except Exception as e:
            raise Exception(f"解析PHYLIP文件时出错: {str(e)}")
    
    def parse_sequential(self, lines, num_seqs):
        """解析序列式PHYLIP格式"""
        for i in range(1, num_seqs + 1):
            if i < len(lines):
                parts = re.split(r'\s+', lines[i], 1)
                if len(parts) < 2:
                    raise ValueError(f"第{i}行格式错误")
                
                seq_id = parts[0]
                sequence = parts[1].replace(" ", "")
                
                self.sequences[seq_id] = sequence
    
    def parse_interleaved(self, lines, num_seqs):
        """解析交错式PHYLIP格式"""
        # 处理第一个块（包含序列ID）
        seq_order = []
        for i in range(1, num_seqs + 1):
            if i < len(lines):
                parts = re.split(r'\s+', lines[i], 1)
                if len(parts) < 2:
                    raise ValueError(f"第{i}行格式错误")
                
                seq_id = parts[0]
                sequence = parts[1].replace(" ", "")
                
                self.sequences[seq_id] = sequence
                seq_order.append(seq_id)
        
        # 处理剩余的块
        current_block = num_seqs + 1
        while current_block < len(lines):
            for i, seq_id in enumerate(seq_order):
                if current_block + i < len(lines):
                    # 移除行中所有空格
                    seq_data = lines[current_block + i].replace(" ", "")
                    self.sequences[seq_id] += seq_data
            
            current_block += num_seqs

class FastaParser(SequenceParser):
    """解析FASTA格式文件的类"""
    
    def __init__(self, filename):
        super().__init__(filename)
        self.parse()
    
    def parse(self):
        """解析FASTA文件"""
        try:
            with open(self.filename, 'r') as f:
                current_id = None
                current_seq = ""
                
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    
                    if line.startswith(">"):
                        # 如果已有序列，保存它
                        if current_id is not None:
                            self.sequences[current_id] = current_seq
                        
                        # 开始新序列
                        current_id = line[1:].split()[0]  # 取ID部分
                        current_seq = ""
                    else:
                        # 添加序列数据
                        current_seq += line.replace(" ", "")
                
                # 保存最后一个序列
                if current_id is not None:
                    self.sequences[current_id] = current_seq
            
            # 确认所有序列长度相同
            if self.sequences:
                seq_lengths = [len(seq) for seq in self.sequences.values()]
                if len(set(seq_lengths)) > 1:
                    print(f"警告: 不同序列的长度不同 {seq_lengths}")
                self.seq_length = seq_lengths[0]
            else:
                raise ValueError("没有找到有效的FASTA序列")
                
        except Exception as e:
            raise Exception(f"解析FASTA文件时出错: {str(e)}")

def detect_file_format(filename):
    """检测文件格式（PHYLIP 或 FASTA）"""
    try:
        with open(filename, 'r') as f:
            first_line = f.readline().strip()
            
        if first_line.startswith(">"):
            return "fasta"
        else:
            try:
                # 尝试解析PHYLIP头部
                parts = re.split(r'\s+', first_line)
                if len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit():
                    return "phylip"
            except:
                pass
        
        # 如果无法确定，默认为FASTA
        return "fasta"
    except:
        raise Exception(f"无法打开文件或检测格式: {filename}")

def identify_4fold_sites(sequences, seq_length, handle_gaps="skip"):
    """
    识别4重简并位点
    
    参数:
    - sequences: 序列字典
    - seq_length: 序列长度
    - handle_gaps: 处理缺口的方式
      - "skip": 跳过含有缺口的密码子（默认）
      - "include": 在某些序列有缺口时仍然考虑该位点
      - "mask": 将含有缺口的位点标记为特殊字符（如'X'）
    """
    
    # 定义4重简并密码子的首两位
    four_fold_codons = {
        'CT': True,  # 亮氨酸 (Leu)
        'GT': True,  # 缬氨酸 (Val) 
        'TC': True,  # 丝氨酸 (Ser)
        'CC': True,  # 脯氨酸 (Pro)
        'AC': True,  # 苏氨酸 (Thr)
        'GC': True,  # 丙氨酸 (Ala)
        'CG': True,  # 精氨酸 (Arg)
        'GG': True   # 甘氨酸 (Gly)
    }
    
    # 确认所有序列长度相同且是3的倍数
    for seq_id, seq in sequences.items():
        if len(seq) != seq_length:
            raise ValueError(f"序列 {seq_id} 长度 ({len(seq)}) 与指定长度 ({seq_length}) 不符")
    
    if seq_length % 3 != 0:
        print(f"警告: 序列长度 ({seq_length}) 不是3的倍数，最后 {seq_length % 3} 个位置将被忽略")
    
    # 找出4重简并位点位置
    four_fold_sites = []
    four_fold_sites_with_gaps = []  # 含有缺口但仍然是4重简并位点的位置
    
    # 遍历每个密码子位置
    for i in range(0, seq_length - 2, 3):
        # 检查所有序列在该位置是否都是4重简并密码子
        all_4fold = True
        has_gaps = False
        
        # 获取所有序列在当前密码子位置的前两个碱基
        first_seq_id = list(sequences.keys())[0]
        codon_prefix = sequences[first_seq_id][i:i+2].upper()
        
        # 如果第一个序列在这个位置有缺口，根据处理方式决定是否跳过
        if '-' in codon_prefix or 'N' in codon_prefix:
            if handle_gaps == "skip":
                continue
            else:
                has_gaps = True
                
        # 如果前两个碱基不在4重简并名单中且不是缺口，跳过
        if not has_gaps and codon_prefix not in four_fold_codons:
            continue
            
        # 检查所有序列在此位置的密码子前两位
        for seq_id, seq in sequences.items():
            if i+2 >= len(seq):  # 确保不会超出序列范围
                all_4fold = False
                break
                
            current_prefix = seq[i:i+2].upper()
            
            # 检查间隔和未知碱基
            if '-' in current_prefix or 'N' in current_prefix:
                has_gaps = True
                if handle_gaps == "skip":
                    all_4fold = False
                    break
            elif not has_gaps:  # 只在没有缺口的情况下比较前缀
                # 如果当前前缀不是4重简并的，或者与第一个序列不同
                if current_prefix not in four_fold_codons or current_prefix != codon_prefix:
                    all_4fold = False
                    break
        
        # 如果是有效的4重简并位点
        if all_4fold:
            if has_gaps and handle_gaps == "mask":
                four_fold_sites_with_gaps.append(i + 2)
            else:
                four_fold_sites.append(i + 2)
    
    return four_fold_sites, four_fold_sites_with_gaps

def extract_4fold_sites(sequences, sites, sites_with_gaps, seq_order, handle_gaps="skip"):
    """从序列中提取指定位置的碱基"""
    extracted_sequences = OrderedDict()
    
    for seq_id in seq_order:
        if handle_gaps == "mask":
            # 提取4重简并位点的碱基，将含有缺口的位点标记为'X'
            extracted_seq = ''.join([sequences[seq_id][i] if i in sites else 'X' for i in sorted(sites + sites_with_gaps)])
        else:
            # 只提取不含缺口的4重简并位点
            extracted_seq = ''.join([sequences[seq_id][i] for i in sorted(sites)])
        
        extracted_sequences[seq_id] = extracted_seq
    
    return extracted_sequences, len(extracted_sequences[seq_order[0]])

def write_phylip(sequences, seq_order, seq_length, output_file):
    """将序列写入PHYLIP格式文件"""
    try:
        with open(output_file, 'w') as f:
            # 写入头部
            f.write(f"{len(sequences)} {seq_length}\n")
            
            # 写入序列
            for seq_id in seq_order:
                # PHYLIP格式通常要求序列ID左对齐且固定宽度
                f.write(f"{seq_id.ljust(10)} {sequences[seq_id]}\n")
                
    except Exception as e:
        raise Exception(f"写入PHYLIP文件时出错: {str(e)}")

def write_fasta(sequences, seq_order, output_file):
    """将序列写入FASTA格式文件"""
    try:
        with open(output_file, 'w') as f:
            for seq_id in seq_order:
                # 写入FASTA格式
                f.write(f">{seq_id}\n")
                
                # 将序列分成固定长度的行（通常是80个字符）
                seq = sequences[seq_id]
                for i in range(0, len(seq), 80):
                    f.write(f"{seq[i:i+80]}\n")
                
    except Exception as e:
        raise Exception(f"写入FASTA文件时出错: {str(e)}")

def main():
    # 解析命令行参数
    parser = argparse.ArgumentParser(description='从序列比对文件中提取4重简并位点')
    parser.add_argument('input_file', help='输入序列文件')
    parser.add_argument('output_file', help='输出序列文件')
    parser.add_argument('--input-format', choices=['phylip', 'fasta', 'auto'], 
                        default='auto', help='输入文件格式 (默认: auto)')
    parser.add_argument('--output-format', choices=['phylip', 'fasta'], 
                        default=None, help='输出文件格式 (默认: 与输入格式相同)')
    parser.add_argument('--include-gaps', action='store_true', 
                        help='在某些序列有缺口时仍然考虑该位点')
    parser.add_argument('--mask-gaps', action='store_true',
                        help='将含有缺口的位点标记为X')
    
    args = parser.parse_args()
    
    # 确定如何处理缺口
    if args.include_gaps:
        handle_gaps = "include"
    elif args.mask_gaps:
        handle_gaps = "mask"
    else:
        handle_gaps = "skip"
    
    try:
        # 检测或使用指定的输入格式
        input_format = args.input_format
        if input_format == 'auto':
            input_format = detect_file_format(args.input_file)
        
        # 使用适当的解析器
        if input_format == 'phylip':
            parser = PhylipParser(args.input_file)
        else:  # fasta
            parser = FastaParser(args.input_file)
        
        # 获取序列
        sequences, seq_order, seq_length = parser.get_sequences()
        
        # 确定输出格式
        output_format = args.output_format if args.output_format else input_format
        
        # 识别4重简并位点
        four_fold_sites, four_fold_sites_with_gaps = identify_4fold_sites(
            sequences, seq_length, handle_gaps)
        
        # 提取4重简并位点
        extracted_sequences, new_length = extract_4fold_sites(
            sequences, four_fold_sites, four_fold_sites_with_gaps, seq_order, handle_gaps)
        
        # 写入输出文件
        if output_format == 'phylip':
            write_phylip(extracted_sequences, seq_order, new_length, args.output_file)
        else:  # fasta
            write_fasta(extracted_sequences, seq_order, args.output_file)
        
        print(f"从 {seq_length} 个位点中识别出 {len(four_fold_sites)} 个不含缺口的4重简并位点")
        if handle_gaps != "skip":
            print(f"另有 {len(four_fold_sites_with_gaps)} 个含缺口的4重简并位点被{('包含' if handle_gaps == 'include' else '标记')}")
        print(f"4重简并位点已保存到: {args.output_file}")
        
        # 输出统计信息
        stats_file = f"{os.path.splitext(args.output_file)[0]}_stats.txt"
        with open(stats_file, 'w') as f:
            f.write(f"输入文件格式: {input_format}\n")
            f.write(f"输出文件格式: {output_format}\n")
            f.write(f"总序列数: {len(sequences)}\n")
            f.write(f"原始序列长度: {seq_length}\n")
            f.write(f"不含缺口的4重简并位点数: {len(four_fold_sites)}\n")
            f.write(f"含缺口的4重简并位点数: {len(four_fold_sites_with_gaps)}\n")
            f.write(f"处理缺口的方式: {handle_gaps}\n")
            f.write(f"不含缺口的4重简并位点位置:\n{', '.join(map(str, sorted(four_fold_sites)))}\n")
            if handle_gaps != "skip":
                f.write(f"含缺口的4重简并位点位置:\n{', '.join(map(str, sorted(four_fold_sites_with_gaps)))}\n")
        
        print(f"统计信息已保存到: {stats_file}")
        
    except Exception as e:
        print(f"错误: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
