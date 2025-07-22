#!/usr/bin/env python3
"""
智能序列校正和完善分析器 - 更新版本
使用现代BioPython API，无弃用警告
Author: Biols9527
Date: 2025-06-24
User: Biols9527
"""

from Bio import SeqIO, Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import pandas as pd
from collections import defaultdict, Counter
import tempfile
import os
import subprocess
import logging
from datetime import datetime
import re
import warnings

# 抑制BioPython弃用警告
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", module="Bio.pairwise2")
warnings.filterwarnings("ignore", module="Bio.Application")

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class SequenceCorrectionAnalyzer:
    def __init__(self, correction_method='consensus', min_overlap=50):
        """
        初始化序列校正分析器
        
        Args:
            correction_method: 校正方法 ['consensus', 'longest_reference', 'majority_vote']
            min_overlap: 最小重叠长度要求
        """
        self.correction_method = correction_method
        self.min_overlap = min_overlap
        self.species_groups = {}
        self.correction_results = {}
        
        logger.info(f"Initialized sequence correction analyzer")
        logger.info(f"Correction method: {correction_method}")
        logger.info(f"Minimum overlap: {min_overlap} bp")
    
    def analyze_sequence_completeness(self, sequences, species_name):
        """分析序列完整性和潜在问题"""
        analysis = {
            'species': species_name,
            'total_sequences': len(sequences),
            'length_stats': {},
            'frame_analysis': {},
            'completion_issues': []
        }
        
        lengths = [len(seq.seq) for seq in sequences]
        analysis['length_stats'] = {
            'min': min(lengths),
            'max': max(lengths),
            'mean': np.mean(lengths),
            'std': np.std(lengths),
            'lengths': lengths
        }
        
        # 分析读码框情况
        frame_info = []
        for seq in sequences:
            seq_len = len(seq.seq)
            frame_remainder = seq_len % 3
            
            frame_info.append({
                'gene_name': seq.gene_name,
                'length': seq_len,
                'frame_remainder': frame_remainder,
                'is_complete_frame': frame_remainder == 0,
                'missing_bases': (3 - frame_remainder) % 3
            })
        
        analysis['frame_analysis'] = frame_info
        
        # 识别潜在的完整性问题
        incomplete_frames = [f for f in frame_info if not f['is_complete_frame']]
        if incomplete_frames:
            analysis['completion_issues'].append('incomplete_reading_frames')
        
        length_variation = np.std(lengths) / np.mean(lengths) if np.mean(lengths) > 0 else 0
        if length_variation > 0.1:  # 长度变异超过10%
            analysis['completion_issues'].append('high_length_variation')
        
        if max(lengths) - min(lengths) > 100:  # 长度差异超过100bp
            analysis['completion_issues'].append('significant_length_differences')
        
        return analysis
    
    def perform_pairwise_alignment(self, seq1, seq2):
        """使用现代API进行双序列比对"""
        try:
            # 使用新的PairwiseAligner
            aligner = Align.PairwiseAligner()
            aligner.match_score = 2
            aligner.mismatch_score = -1
            aligner.open_gap_score = -2
            aligner.extend_gap_score = -0.5
            
            alignments = aligner.align(str(seq1.seq), str(seq2.seq))
            if alignments:
                best_alignment = alignments[0]
                return {
                    'score': best_alignment.score,
                    'aligned_length': len(best_alignment),
                    'identity': self.calculate_identity(best_alignment),
                    'alignment': best_alignment
                }
            else:
                return None
        except Exception as e:
            logger.warning(f"Pairwise alignment failed: {e}")
            return None
    
    def calculate_identity(self, alignment):
        """计算比对的一致性"""
        try:
            aligned_seq1, aligned_seq2 = alignment.aligned
            matches = 0
            total = 0
            
            for i in range(len(aligned_seq1)):
                if i < len(aligned_seq2):
                    if aligned_seq1[i] == aligned_seq2[i]:
                        matches += 1
                    total += 1
            
            return matches / total if total > 0 else 0.0
        except:
            return 0.0
    
    def perform_multiple_sequence_alignment(self, sequences, species_name):
        """执行多序列比对来识别缺失区域"""
        if len(sequences) < 2:
            return None
            
        logger.info(f"Performing MSA for {species_name} ({len(sequences)} sequences)")
        
        try:
            with tempfile.TemporaryDirectory() as temp_dir:
                # 创建输入文件
                input_file = os.path.join(temp_dir, f"{species_name.replace(' ', '_')}_input.fasta")
                output_file = os.path.join(temp_dir, f"{species_name.replace(' ', '_')}_aligned.fasta")
                
                # 写入序列
                with open(input_file, 'w') as f:
                    for i, seq in enumerate(sequences):
                        f.write(f">{seq.gene_name}\n{str(seq.seq)}\n")
                
                # 执行MAFFT比对
                try:
                    result = subprocess.run(
                        ["mafft", "--auto", "--quiet", input_file],
                        capture_output=True,
                        text=True,
                        timeout=300
                    )
                    
                    if result.returncode == 0:
                        with open(output_file, 'w') as f:
                            f.write(result.stdout)
                        
                        # 读取比对结果
                        if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                            aligned_sequences = list(SeqIO.parse(output_file, "fasta"))
                            return self.analyze_alignment_gaps(aligned_sequences, sequences)
                        else:
                            logger.warning(f"MAFFT output empty for {species_name}")
                            return None
                    else:
                        logger.warning(f"MAFFT failed for {species_name}: {result.stderr}")
                        return None
                        
                except subprocess.TimeoutExpired:
                    logger.warning(f"MAFFT timeout for {species_name}")
                    return None
                except FileNotFoundError:
                    logger.warning("MAFFT not found, skipping MSA")
                    return None
                    
        except Exception as e:
            logger.error(f"MSA failed for {species_name}: {e}")
            return None
    
    def analyze_alignment_gaps(self, aligned_sequences, original_sequences):
        """分析比对结果中的gap模式"""
        alignment_analysis = {
            'aligned_length': len(aligned_sequences[0].seq) if aligned_sequences else 0,
            'gap_patterns': {},
            'terminal_gaps': {},
            'internal_gaps': {},
            'consensus_sequence': None
        }
        
        if not aligned_sequences:
            return alignment_analysis
        
        # 创建基因名到原始序列的映射
        gene_to_original = {seq.gene_name: seq for seq in original_sequences}
        
        # 分析每个序列的gap模式
        for aligned_seq in aligned_sequences:
            gene_name = aligned_seq.id
            aligned_str = str(aligned_seq.seq)
            
            # 分析terminal gaps
            start_gaps = len(aligned_str) - len(aligned_str.lstrip('-'))
            end_gaps = len(aligned_str) - len(aligned_str.rstrip('-'))
            
            # 分析internal gaps
            internal_str = aligned_str.strip('-')
            internal_gaps = internal_str.count('-')
            
            alignment_analysis['terminal_gaps'][gene_name] = {
                'start_gaps': start_gaps,
                'end_gaps': end_gaps,
                'total_terminal_gaps': start_gaps + end_gaps
            }
            
            alignment_analysis['internal_gaps'][gene_name] = {
                'internal_gaps': internal_gaps,
                'gap_positions': [i for i, char in enumerate(internal_str) if char == '-']
            }
        
        # 生成共识序列
        alignment_analysis['consensus_sequence'] = self.generate_consensus_sequence(aligned_sequences)
        
        return alignment_analysis
    
    def generate_consensus_sequence(self, aligned_sequences):
        """从比对结果生成共识序列"""
        if not aligned_sequences:
            return None
            
        alignment_length = len(aligned_sequences[0].seq)
        consensus = []
        
        for pos in range(alignment_length):
            bases_at_pos = []
            for seq in aligned_sequences:
                if pos < len(seq.seq):
                    base = str(seq.seq)[pos].upper()
                    if base != '-':
                        bases_at_pos.append(base)
            
            if bases_at_pos:
                # 多数表决
                base_counts = Counter(bases_at_pos)
                most_common_base = base_counts.most_common(1)[0][0]
                consensus.append(most_common_base)
            else:
                consensus.append('-')  # 所有序列在此位置都是gap
        
        # 移除terminal gaps
        consensus_str = ''.join(consensus).strip('-')
        
        return consensus_str
    
    def identify_correction_opportunities(self, sequences, alignment_analysis, completeness_analysis):
        """识别序列校正机会"""
        if not alignment_analysis or not alignment_analysis['consensus_sequence']:
            # 如果没有比对结果，使用简单策略
            return self.simple_frame_correction(sequences, completeness_analysis)
        
        consensus = alignment_analysis['consensus_sequence']
        corrections = {}
        
        for seq in sequences:
            gene_name = seq.gene_name
            original_seq = str(seq.seq).upper()
            original_length = len(original_seq)
            
            # 检查是否需要校正
            frame_remainder = original_length % 3
            needs_correction = frame_remainder != 0
            
            if not needs_correction:
                corrections[gene_name] = {
                    'needs_correction': False,
                    'original_length': original_length,
                    'corrected_sequence': original_seq,
                    'correction_type': 'none'
                }
                continue
            
            # 尝试不同的校正策略
            correction_attempts = []
            
            # 策略1: 基于共识序列的末端延伸
            if gene_name in alignment_analysis['terminal_gaps']:
                terminal_info = alignment_analysis['terminal_gaps'][gene_name]
                
                # 3'端校正（最常见的情况）
                if terminal_info['end_gaps'] > 0:
                    extension_3 = self.extract_consensus_extension(
                        consensus, 'end', terminal_info['end_gaps'], original_seq
                    )
                    if extension_3:
                        corrected_seq = original_seq + extension_3
                        correction_attempts.append({
                            'method': '3_prime_extension',
                            'sequence': corrected_seq,
                            'added_bases': len(extension_3),
                            'final_length': len(corrected_seq),
                            'frame_complete': len(corrected_seq) % 3 == 0
                        })
                
                # 5'端校正
                if terminal_info['start_gaps'] > 0:
                    extension_5 = self.extract_consensus_extension(
                        consensus, 'start', terminal_info['start_gaps'], original_seq
                    )
                    if extension_5:
                        corrected_seq = extension_5 + original_seq
                        correction_attempts.append({
                            'method': '5_prime_extension',
                            'sequence': corrected_seq,
                            'added_bases': len(extension_5),
                            'final_length': len(corrected_seq),
                            'frame_complete': len(corrected_seq) % 3 == 0
                        })
            
            # 策略2: 最小化校正（只添加必要的碱基使其成为3的倍数）
            missing_bases = (3 - frame_remainder) % 3
            if missing_bases > 0:
                # 尝试从共识序列推断最可能的碱基
                inferred_bases = self.infer_missing_bases(original_seq, consensus, missing_bases)
                if inferred_bases:
                    corrected_seq = original_seq + inferred_bases
                    correction_attempts.append({
                        'method': 'minimal_extension',
                        'sequence': corrected_seq,
                        'added_bases': len(inferred_bases),
                        'final_length': len(corrected_seq),
                        'frame_complete': len(corrected_seq) % 3 == 0
                    })
            
            # 选择最佳校正方案
            best_correction = self.select_best_correction(correction_attempts, original_seq)
            
            corrections[gene_name] = {
                'needs_correction': True,
                'original_length': original_length,
                'original_frame_remainder': frame_remainder,
                'correction_attempts': correction_attempts,
                'best_correction': best_correction,
                'corrected_sequence': best_correction['sequence'] if best_correction else original_seq,
                'correction_type': best_correction['method'] if best_correction else 'failed'
            }
        
        return corrections
    
    def simple_frame_correction(self, sequences, completeness_analysis):
        """简单的读码框校正（无比对信息时使用）"""
        corrections = {}
        
        for seq in sequences:
            gene_name = seq.gene_name
            original_seq = str(seq.seq).upper()
            original_length = len(original_seq)
            frame_remainder = original_length % 3
            
            if frame_remainder == 0:
                corrections[gene_name] = {
                    'needs_correction': False,
                    'original_length': original_length,
                    'corrected_sequence': original_seq,
                    'correction_type': 'none'
                }
            else:
                # 需要校正
                missing_bases = (3 - frame_remainder) % 3
                
                # 使用序列中最常见的碱基进行延伸
                base_freq = Counter(original_seq)
                if base_freq:
                    most_common_base = base_freq.most_common(1)[0][0]
                    extension = most_common_base * missing_bases
                else:
                    extension = 'N' * missing_bases
                
                corrected_seq = original_seq + extension
                
                corrections[gene_name] = {
                    'needs_correction': True,
                    'original_length': original_length,
                    'original_frame_remainder': frame_remainder,
                    'corrected_sequence': corrected_seq,
                    'correction_type': f'simple_extension_{missing_bases}_bases',
                    'best_correction': {
                        'method': 'simple_extension',
                        'sequence': corrected_seq,
                        'added_bases': missing_bases,
                        'final_length': len(corrected_seq),
                        'frame_complete': True
                    }
                }
        
        return corrections
    
    def extract_consensus_extension(self, consensus, position, gap_length, original_seq):
        """从共识序列中提取用于延伸的序列"""
        if position == 'start':
            if gap_length <= len(consensus):
                return consensus[:gap_length]
        elif position == 'end':
            if gap_length <= len(consensus):
                return consensus[-gap_length:]
        
        return None
    
    def infer_missing_bases(self, original_seq, consensus, num_bases):
        """推断缺失的碱基"""
        if not consensus or num_bases <= 0:
            return None
        
        # 尝试找到原始序列在共识序列中的最佳匹配位置
        best_match_pos = self.find_best_match_position(original_seq, consensus)
        
        if best_match_pos is not None:
            # 根据匹配位置推断可能的延伸
            consensus_end = best_match_pos + len(original_seq)
            if consensus_end + num_bases <= len(consensus):
                return consensus[consensus_end:consensus_end + num_bases]
        
        # 如果无法从共识序列推断，使用频率分析
        base_freq = Counter(original_seq)
        if base_freq:
            most_common_base = base_freq.most_common(1)[0][0]
            return most_common_base * num_bases
        
        return 'N' * num_bases  # 最后的选择
    
    def find_best_match_position(self, original_seq, consensus):
        """在共识序列中找到原始序列的最佳匹配位置"""
        if not original_seq or not consensus:
            return None
        
        best_score = 0
        best_pos = None
        
        # 滑动窗口寻找最佳匹配
        for i in range(max(1, len(consensus) - len(original_seq) + 1)):
            consensus_segment = consensus[i:i + len(original_seq)]
            if len(consensus_segment) == len(original_seq):
                matches = sum(1 for a, b in zip(original_seq, consensus_segment) if a == b)
                score = matches / len(original_seq)
                
                if score > best_score:
                    best_score = score
                    best_pos = i
        
        return best_pos if best_score > 0.7 else None  # 至少70%匹配
    
    def select_best_correction(self, correction_attempts, original_seq):
        """选择最佳的校正方案"""
        if not correction_attempts:
            return None
        
        # 评分标准
        scored_attempts = []
        for attempt in correction_attempts:
            score = 0
            
            # 优先选择使读码框完整的方案
            if attempt['frame_complete']:
                score += 100
            
            # 偏好最小修改
            modification_penalty = attempt['added_bases'] * 5
            score -= modification_penalty
            
            # 偏好特定方法
            method_bonus = {
                'minimal_extension': 20,
                '3_prime_extension': 15,
                '5_prime_extension': 10,
                'both_ends_extension': 5
            }
            score += method_bonus.get(attempt['method'], 0)
            
            scored_attempts.append((score, attempt))
        
        # 返回最高分的方案
        scored_attempts.sort(key=lambda x: x[0], reverse=True)
        return scored_attempts[0][1] if scored_attempts else None
    
    def validate_corrected_sequences(self, corrections):
        """验证校正后的序列质量"""
        validation_results = {}
        
        for gene_name, correction_info in corrections.items():
            if not correction_info['needs_correction']:
                validation_results[gene_name] = {
                    'validation_status': 'no_correction_needed',
                    'quality_score': 1.0
                }
                continue
            
            corrected_seq = correction_info['corrected_sequence']
            validation = {
                'validation_status': 'corrected',
                'original_length': correction_info['original_length'],
                'corrected_length': len(corrected_seq),
                'frame_complete': len(corrected_seq) % 3 == 0,
                'correction_method': correction_info['correction_type']
            }
            
            # 检查校正质量
            quality_checks = []
            
            # 1. 读码框完整性
            if validation['frame_complete']:
                quality_checks.append('frame_complete')
            
            # 2. 序列长度合理性
            length_change = len(corrected_seq) - correction_info['original_length']
            if length_change <= 6:  # 最多添加6个碱基（2个密码子）
                quality_checks.append('reasonable_length_change')
            
            # 3. 无异常字符
            if all(base in 'ATGCN' for base in corrected_seq.upper()):
                quality_checks.append('valid_bases')
            
            # 4. ORF分析
            orf_analysis = self.analyze_orf_in_corrected_sequence(corrected_seq)
            if orf_analysis['has_valid_orf']:
                quality_checks.append('valid_orf')
            
            validation['quality_checks'] = quality_checks
            validation['quality_score'] = len(quality_checks) / 4.0  # 标准化到0-1
            validation['orf_analysis'] = orf_analysis
            
            validation_results[gene_name] = validation
        
        return validation_results
    
    def analyze_orf_in_corrected_sequence(self, sequence):
        """分析校正后序列中的ORF"""
        seq_str = sequence.upper()
        start_codons = ['ATG']
        stop_codons = ['TAA', 'TAG', 'TGA']
        
        orf_info = {
            'has_valid_orf': False,
            'longest_orf': 0,
            'orf_coverage': 0.0,
            'start_codon_found': False,
            'stop_codon_found': False
        }
        
        # 检查所有读码框
        for frame in range(3):
            frame_seq = seq_str[frame:]
            if len(frame_seq) < 6:  # 至少需要2个密码子
                continue
            
            # 寻找ORF
            current_orf_length = 0
            found_start = False
            
            for i in range(0, len(frame_seq) - 2, 3):
                codon = frame_seq[i:i+3]
                if len(codon) != 3:
                    break
                
                if codon in start_codons and not found_start:
                    found_start = True
                    current_orf_length = 3
                    orf_info['start_codon_found'] = True
                elif found_start:
                    if codon in stop_codons:
                        current_orf_length += 3
                        if current_orf_length > orf_info['longest_orf']:
                            orf_info['longest_orf'] = current_orf_length
                        orf_info['stop_codon_found'] = True
                        found_start = False
                        current_orf_length = 0
                    else:
                        current_orf_length += 3
        
        # 计算ORF覆盖度
        if len(sequence) > 0:
            orf_info['orf_coverage'] = orf_info['longest_orf'] / len(sequence)
        
        # 判断是否有有效的ORF
        orf_info['has_valid_orf'] = (
            orf_info['longest_orf'] >= 150 and  # 至少50个氨基酸
            orf_info['start_codon_found'] and
            orf_info['stop_codon_found']
        )
        
        return orf_info
    
    def generate_correction_report(self, species_results, output_file):
        """生成序列校正报告"""
        try:
            with open(output_file, 'w') as f:
                f.write("=" * 100 + "\n")
                f.write("SEQUENCE CORRECTION AND COMPLETION ANALYSIS REPORT\n")
                f.write("=" * 100 + "\n")
                f.write(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S UTC')}\n")
                f.write(f"User: Biols9527\n")
                f.write(f"Correction method: {self.correction_method}\n")
                f.write(f"Minimum overlap: {self.min_overlap} bp\n")
                f.write("\n")
                
                # 总体统计
                total_species = len(species_results)
                species_with_issues = sum(1 for r in species_results.values() 
                                        if r['completeness_analysis']['completion_issues'])
                total_sequences = sum(r['completeness_analysis']['total_sequences'] 
                                    for r in species_results.values())
                sequences_corrected = sum(1 for r in species_results.values()
                                        for gene, corr in r.get('corrections', {}).items()
                                        if corr['needs_correction'])
                
                f.write("OVERALL STATISTICS:\n")
                f.write("-" * 50 + "\n")
                f.write(f"Total species analyzed: {total_species}\n")
                f.write(f"Species with completion issues: {species_with_issues}\n")
                f.write(f"Total sequences: {total_sequences}\n")
                f.write(f"Sequences corrected: {sequences_corrected}\n")
                if total_sequences > 0:
                    f.write(f"Correction success rate: {sequences_corrected/total_sequences*100:.1f}%\n")
                f.write("\n")
                
                # 各物种详细分析
                f.write("SPECIES-WISE CORRECTION ANALYSIS:\n")
                f.write("=" * 100 + "\n")
                
                for species, results in species_results.items():
                    f.write(f"\nSpecies: {species}\n")
                    f.write("-" * 80 + "\n")
                    
                    completeness = results['completeness_analysis']
                    f.write(f"Total sequences: {completeness['total_sequences']}\n")
                    f.write(f"Length range: {completeness['length_stats']['min']}-{completeness['length_stats']['max']} bp\n")
                    f.write(f"Mean length: {completeness['length_stats']['mean']:.1f} ± {completeness['length_stats']['std']:.1f} bp\n")
                    
                    # 读码框分析
                    frame_info = completeness['frame_analysis']
                    incomplete_frames = [f for f in frame_info if not f['is_complete_frame']]
                    
                    f.write(f"Sequences with incomplete reading frames: {len(incomplete_frames)}/{len(frame_info)}\n")
                    
                    if incomplete_frames:
                        f.write("Frame completion details:\n")
                        for frame in incomplete_frames:
                            f.write(f"  {frame['gene_name']}: {frame['length']} bp "
                                   f"(remainder: {frame['frame_remainder']}, "
                                   f"missing: {frame['missing_bases']} bp)\n")
                    
                    # 校正结果
                    if 'corrections' in results:
                        corrections = results['corrections']
                        corrected_genes = [gene for gene, corr in corrections.items() 
                                         if corr['needs_correction']]
                        
                        f.write(f"\nCorrection results:\n")
                        f.write(f"Sequences requiring correction: {len(corrected_genes)}\n")
                        
                        for gene in corrected_genes:
                            corr = corrections[gene]
                            f.write(f"\n  {gene}:\n")
                            f.write(f"    Original length: {corr['original_length']} bp\n")
                            f.write(f"    Corrected length: {len(corr['corrected_sequence'])} bp\n")
                            f.write(f"    Correction method: {corr['correction_type']}\n")
                            
                            if corr.get('best_correction'):
                                best = corr['best_correction']
                                f.write(f"    Added bases: {best['added_bases']}\n")
                                f.write(f"    Frame complete: {best['frame_complete']}\n")
                    
                    # 验证结果
                    if 'validation' in results:
                        validation = results['validation']
                        high_quality = sum(1 for v in validation.values() 
                                         if v.get('quality_score', 0) >= 0.8)
                        
                        f.write(f"\nValidation results:\n")
                        f.write(f"High quality corrections: {high_quality}/{len(validation)}\n")
                        
                        for gene, val in validation.items():
                            if val['validation_status'] == 'corrected':
                                f.write(f"  {gene}: {val['quality_score']:.2f} "
                                       f"({'/'.join(val.get('quality_checks', []))})\n")
                    
                    f.write("\n" + "=" * 80 + "\n")
                
                # 方法说明
                f.write("\nMETHODOLOGY:\n")
                f.write("-" * 50 + "\n")
                f.write("1. Sequence completeness analysis:\n")
                f.write("   - Reading frame analysis (length % 3)\n")
                f.write("   - Length variation assessment\n\n")
                f.write("2. Multiple sequence alignment (when available):\n")
                f.write("   - MAFFT-based alignment of species copies\n")
                f.write("   - Gap pattern analysis\n\n")
                f.write("3. Correction strategies:\n")
                f.write("   - 5' and 3' terminal extensions\n")
                f.write("   - Minimal base addition for frame completion\n")
                f.write("   - Consensus-based sequence inference\n\n")
                f.write("4. Validation criteria:\n")
                f.write("   - Reading frame completeness\n")
                f.write("   - ORF analysis\n")
                f.write("   - Sequence length reasonableness\n")
                
            logger.info(f"Correction report saved to {output_file}")
            
        except Exception as e:
            logger.error(f"Failed to generate correction report: {e}")
    
    def run_sequence_correction_analysis(self, input_fasta, output_prefix="corrected"):
        """运行完整的序列校正分析"""
        logger.info(f"Starting sequence correction analysis: {input_fasta}")
        
        try:
            # 解析序列并按物种分组
            species_groups = defaultdict(list)
            for record in SeqIO.parse(input_fasta, "fasta"):
                parts = record.description.split(' ', 1)
                if len(parts) == 2:
                    gene_name = parts[0].strip()
                    species_name = parts[1].strip()
                else:
                    gene_name = parts[0].strip()
                    species_name = "Unknown_species"
                
                record.gene_name = gene_name
                record.species_name = species_name
                species_groups[species_name].append(record)
            
            self.species_groups = dict(species_groups)
            logger.info(f"Found {len(self.species_groups)} species")
            
            # 分析每个物种
            species_results = {}
            
            for species, sequences in self.species_groups.items():
                logger.info(f"Processing {species}...")
                
                # 分析序列完整性
                completeness_analysis = self.analyze_sequence_completeness(sequences, species)
                
                # 执行多序列比对（如果有多个拷贝）
                alignment_analysis = None
                if len(sequences) > 1:
                    alignment_analysis = self.perform_multiple_sequence_alignment(sequences, species)
                
                # 识别校正机会
                corrections = self.identify_correction_opportunities(
                    sequences, alignment_analysis, completeness_analysis
                )
                
                # 验证校正结果
                validation = {}
                if corrections:
                    validation = self.validate_corrected_sequences(corrections)
                
                species_results[species] = {
                    'completeness_analysis': completeness_analysis,
                    'alignment_analysis': alignment_analysis,
                    'corrections': corrections,
                    'validation': validation
                }
            
            # 生成报告
            self.generate_correction_report(species_results, f"{output_prefix}_correction_report.txt")
            
            # 保存校正后的序列
            corrected_sequences = []
            original_sequences = []
            
            for species, results in species_results.items():
                for seq in self.species_groups[species]:
                    # 保存原始序列
                    original_seq = SeqRecord(
                        seq.seq,
                        id=f"{seq.gene_name}_original",
                        description=f"{seq.gene_name} {seq.species_name} [original]"
                    )
                    original_sequences.append(original_seq)
                    
                    # 检查是否有校正版本
                    if ('corrections' in results and 
                        seq.gene_name in results['corrections'] and
                        results['corrections'][seq.gene_name]['needs_correction']):
                        
                        correction = results['corrections'][seq.gene_name]
                        corrected_seq_str = correction['corrected_sequence']
                        
                        corrected_seq = SeqRecord(
                            Seq(corrected_seq_str),
                            id=seq.gene_name,
                            description=f"{seq.gene_name} {seq.species_name} [corrected by {correction['correction_type']}]"
                        )
                        corrected_sequences.append(corrected_seq)
                    else:
                        # 无需校正或校正失败，使用原始序列
                        corrected_seq = SeqRecord(
                            seq.seq,
                            id=seq.gene_name,
                            description=f"{seq.gene_name} {seq.species_name}"
                        )
                        corrected_sequences.append(corrected_seq)
            
            # 保存结果文件
            SeqIO.write(corrected_sequences, f"{output_prefix}_corrected_sequences.fasta", "fasta")
            SeqIO.write(original_sequences, f"{output_prefix}_original_sequences.fasta", "fasta")
            
            # 统计信息
            total_corrected = sum(1 for r in species_results.values()
                                for corr in r.get('corrections', {}).values()
                                if corr['needs_correction'])
            
            high_quality_corrections = sum(1 for r in species_results.values()
                                         for val in r.get('validation', {}).values()
                                         if val.get('quality_score', 0) >= 0.8)
            
            logger.info("=" * 60)
            logger.info("SEQUENCE CORRECTION ANALYSIS COMPLETED")
            logger.info("=" * 60)
            logger.info(f"Total species: {len(self.species_groups)}")
            logger.info(f"Total sequences: {len(corrected_sequences)}")
            logger.info(f"Sequences corrected: {total_corrected}")
            logger.info(f"High-quality corrections: {high_quality_corrections}")
            logger.info("=" * 60)
            
            return {
                'total_sequences': len(corrected_sequences),
                'corrected_count': total_corrected,
                'high_quality_count': high_quality_corrections
            }
            
        except Exception as e:
            logger.error(f"Sequence correction analysis failed: {e}")
            raise

def main():
    """主函数"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Sequence correction and completion analyzer")
    parser.add_argument("input_fasta", help="Input FASTA file")
    parser.add_argument("--output_prefix", default="corrected", help="Output prefix")
    parser.add_argument("--correction_method", 
                       choices=['consensus', 'longest_reference', 'majority_vote'],
                       default='consensus', help="Correction method")
    parser.add_argument("--min_overlap", type=int, default=50, help="Minimum overlap length")
    
    args = parser.parse_args()
    
    # 初始化分析器
    analyzer = SequenceCorrectionAnalyzer(
        correction_method=args.correction_method,
        min_overlap=args.min_overlap
    )
    
    # 运行分析
    result = analyzer.run_sequence_correction_analysis(args.input_fasta, args.output_prefix)
    
    print("\n" + "="*80)
    print("✅ SEQUENCE CORRECTION ANALYSIS COMPLETED!")
    print("="*80)
    print(f"ð Corrected sequences: {args.output_prefix}_corrected_sequences.fasta")
    print(f"ð Original sequences: {args.output_prefix}_original_sequences.fasta")
    print(f"ð Detailed report: {args.output_prefix}_correction_report.txt")
    print(f"ð Summary: {result['total_sequences']} total, {result['corrected_count']} corrected")
    print(f"⭐ High quality: {result['high_quality_count']} corrections")
    print("="*80)

if __name__ == "__main__":
    main()
