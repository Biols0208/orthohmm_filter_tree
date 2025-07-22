#!/usr/bin/env python3
"""
全局BLAST质量分析器
针对同一基因家族的序列进行全局比对分析
Author: Biols9527
Date: 2025-06-24
User: Biols9527
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction
import numpy as np
import pandas as pd
from collections import defaultdict, Counter
import tempfile
import os
import subprocess
import logging
from datetime import datetime
import time
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
import warnings

# 抑制BioPython弃用警告
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", module="Bio.Application")

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class GlobalBlastQualityAnalyzer:
    def __init__(self, sequence_type='nucleotide', weights=None, blast_params=None):
        """
        初始化全局BLAST质量分析器
        
        Args:
            sequence_type: 'nucleotide' 或 'protein'
            weights: 各质量指标的权重字典
            blast_params: BLAST参数配置
        """
        self.sequence_type = sequence_type.lower()
        
        # 默认权重配置（强化BLAST一致性权重）
        if weights is None:
            if self.sequence_type == 'nucleotide':
                self.weights = {
                    'length_quality': 0.15,         # 长度质量
                    'composition_quality': 0.10,     # 组成质量
                    'sequence_integrity': 0.20,      # 序列完整性
                    'coding_potential': 0.10,        # 编码潜力
                    'complexity': 0.10,              # 序列复杂度
                    'global_blast_consistency': 0.25, # 全局BLAST一致性
                    'phylogenetic_position': 0.10    # 系统发育位置
                }
            else:  # protein
                self.weights = {
                    'length_quality': 0.15,
                    'composition_quality': 0.15,
                    'sequence_integrity': 0.20,
                    'structural_features': 0.15,
                    'global_blast_consistency': 0.25,
                    'phylogenetic_position': 0.10
                }
        else:
            self.weights = weights
            
        # 标准化权重
        total_weight = sum(self.weights.values())
        self.weights = {k: v/total_weight for k, v in self.weights.items()}
        
        # BLAST参数配置
        if blast_params is None:
            self.blast_params = {
                'evalue': 1e-10,  # 更严格的阈值
                'word_size': 11 if sequence_type == 'nucleotide' else 3,
                'max_target_seqs': 1000,  # 允许更多目标序列
                'outfmt': 6,  # 表格格式
                'num_threads': min(mp.cpu_count(), 8)
            }
        else:
            self.blast_params = blast_params
            
        self.all_sequences = []
        self.species_groups = {}
        self.global_blast_results = {}
        
        logger.info(f"Initialized global BLAST analyzer for {sequence_type} sequences")
        logger.info(f"Quality weights: {self.weights}")
        logger.info(f"BLAST parameters: {self.blast_params}")
    
    def create_global_blast_database(self, all_sequences, temp_dir):
        """为所有序列创建全局BLAST数据库"""
        logger.info(f"Creating global BLAST database with {len(all_sequences)} sequences")
        
        # 创建全局数据库文件
        db_file = os.path.join(temp_dir, "global_gene_family_db")
        fasta_file = f"{db_file}.fasta"
        
        # 写入所有序列
        with open(fasta_file, 'w') as f:
            for i, seq in enumerate(all_sequences):
                # 使用包含物种信息的ID
                seq_id = f"seq_{i:04d}_{seq.gene_name}_{seq.species_name.replace(' ', '_')}"
                f.write(f">{seq_id}\n{str(seq.seq)}\n")
        
        # 创建BLAST数据库
        try:
            if self.sequence_type == 'nucleotide':
                cmd = [
                    "makeblastdb",
                    "-in", fasta_file,
                    "-dbtype", "nucl",
                    "-out", db_file,
                    "-title", "Global_Gene_Family_Database"
                ]
            else:
                cmd = [
                    "makeblastdb",
                    "-in", fasta_file,
                    "-dbtype", "prot",
                    "-out", db_file,
                    "-title", "Global_Gene_Family_Database"
                ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            
            if result.returncode == 0:
                expected_file = f"{db_file}.nin" if self.sequence_type == 'nucleotide' else f"{db_file}.pin"
                if os.path.exists(expected_file):
                    logger.info(f"Created global BLAST database successfully")
                    return db_file, fasta_file
                else:
                    logger.error(f"Global BLAST database files not found")
                    return None
            else:
                logger.error(f"makeblastdb failed: {result.stderr}")
                return None
                
        except FileNotFoundError:
            logger.error("makeblastdb not found. Please install BLAST+")
            return None
        except subprocess.TimeoutExpired:
            logger.error(f"makeblastdb timeout")
            return None
        except Exception as e:
            logger.error(f"Error creating global BLAST database: {e}")
            return None
    
    def perform_global_blast_analysis(self, all_sequences, temp_dir):
        """执行全局BLAST分析"""
        logger.info(f"Performing global BLAST analysis on {len(all_sequences)} sequences")
        
        db_info = self.create_global_blast_database(all_sequences, temp_dir)
        if not db_info:
            return {}
            
        db_file, fasta_file = db_info
        blast_results = {}
        
        try:
            # 执行全局BLAST比对
            blast_output = os.path.join(temp_dir, "global_blast_results.txt")
            
            if self.sequence_type == 'nucleotide':
                blast_cmd = [
                    "blastn",
                    "-query", fasta_file,
                    "-db", db_file,
                    "-evalue", str(self.blast_params['evalue']),
                    "-outfmt", str(self.blast_params['outfmt']),
                    "-out", blast_output,
                    "-word_size", str(self.blast_params['word_size']),
                    "-max_target_seqs", str(self.blast_params['max_target_seqs']),
                    "-num_threads", str(self.blast_params['num_threads']),
                    "-perc_identity", "50"  # 最低50%相似性
                ]
            else:
                blast_cmd = [
                    "blastp",
                    "-query", fasta_file,
                    "-db", db_file,
                    "-evalue", str(self.blast_params['evalue']),
                    "-outfmt", str(self.blast_params['outfmt']),
                    "-out", blast_output,
                    "-word_size", str(self.blast_params['word_size']),
                    "-max_target_seqs", str(self.blast_params['max_target_seqs']),
                    "-num_threads", str(self.blast_params['num_threads'])
                ]
            
            result = subprocess.run(blast_cmd, capture_output=True, text=True, timeout=600)
            
            if result.returncode == 0 and os.path.exists(blast_output):
                blast_results = self.parse_global_blast_results(blast_output, all_sequences)
                logger.info(f"Completed global BLAST analysis")
            else:
                logger.error(f"Global BLAST failed: {result.stderr}")
                
        except FileNotFoundError:
            logger.error("BLAST not found. Please install BLAST+")
        except subprocess.TimeoutExpired:
            logger.error(f"Global BLAST timeout")
        except Exception as e:
            logger.error(f"Global BLAST analysis failed: {e}")
            
        return blast_results
    
    def parse_global_blast_results(self, blast_output, all_sequences):
        """解析全局BLAST结果"""
        blast_results = defaultdict(list)
        
        # 创建序列ID映射
        seq_id_map = {}
        for i, seq in enumerate(all_sequences):
            seq_id = f"seq_{i:04d}_{seq.gene_name}_{seq.species_name.replace(' ', '_')}"
            seq_id_map[seq_id] = {
                'gene_name': seq.gene_name,
                'species_name': seq.species_name,
                'sequence': seq
            }
        
        try:
            with open(blast_output, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    
                    fields = line.strip().split('\t')
                    if len(fields) >= 12:
                        # BLAST表格格式字段
                        query_id = fields[0]
                        subject_id = fields[1]
                        identity = float(fields[2]) / 100.0
                        align_length = int(fields[3])
                        mismatches = int(fields[4])
                        gap_opens = int(fields[5])
                        query_start = int(fields[6])
                        query_end = int(fields[7])
                        subject_start = int(fields[8])
                        subject_end = int(fields[9])
                        evalue = float(fields[10])
                        bit_score = float(fields[11])
                        
                        if query_id not in seq_id_map or subject_id not in seq_id_map:
                            continue
                        
                        query_info = seq_id_map[query_id]
                        subject_info = seq_id_map[subject_id]
                        
                        # 跳过自比对
                        if query_id == subject_id:
                            continue
                        
                        # 计算覆盖度
                        query_seq = query_info['sequence']
                        subject_seq = subject_info['sequence']
                        
                        query_coverage = align_length / len(query_seq.seq)
                        subject_coverage = align_length / len(subject_seq.seq)
                        
                        # 判断是否为同物种比对
                        is_intraspecies = query_info['species_name'] == subject_info['species_name']
                        
                        blast_results[query_info['gene_name']].append({
                            'hit_gene': subject_info['gene_name'],
                            'hit_species': subject_info['species_name'],
                            'query_species': query_info['species_name'],
                            'is_intraspecies': is_intraspecies,
                            'identity': identity,
                            'coverage_query': query_coverage,
                            'coverage_subject': subject_coverage,
                            'align_length': align_length,
                            'evalue': evalue,
                            'bit_score': bit_score,
                            'gaps': gap_opens,
                            'query_start': query_start,
                            'query_end': query_end,
                            'subject_start': subject_start,
                            'subject_end': subject_end
                        })
                        
        except Exception as e:
            logger.error(f"Error parsing global BLAST results: {e}")
            
        return dict(blast_results)
    
    def analyze_global_blast_consistency(self, gene_name, blast_results, all_sequences):
        """分析序列在整个基因家族中的BLAST一致性"""
        if gene_name not in blast_results or not blast_results[gene_name]:
            return {
                'global_consistency_score': 0.3,
                'mean_identity': 0.0,
                'mean_coverage': 0.0,
                'num_significant_hits': 0,
                'intraspecies_hits': 0,
                'interspecies_hits': 0,
                'consistency_rank': 'isolated'
            }
        
        hits = blast_results[gene_name]
        
        # 过滤显著比对
        significant_hits = [hit for hit in hits if hit['evalue'] <= self.blast_params['evalue']]
        
        if not significant_hits:
            return {
                'global_consistency_score': 0.2,
                'mean_identity': 0.0,
                'mean_coverage': 0.0,
                'num_significant_hits': 0,
                'intraspecies_hits': 0,
                'interspecies_hits': 0,
                'consistency_rank': 'poor'
            }
        
        # 分类比对：种内 vs 种间
        intraspecies_hits = [hit for hit in significant_hits if hit['is_intraspecies']]
        interspecies_hits = [hit for hit in significant_hits if not hit['is_intraspecies']]
        
        # 计算各类比对的统计量
        all_identities = [hit['identity'] for hit in significant_hits]
        all_coverages = [hit['coverage_query'] for hit in significant_hits]
        
        intra_identities = [hit['identity'] for hit in intraspecies_hits]
        inter_identities = [hit['identity'] for hit in interspecies_hits]
        
        intra_coverages = [hit['coverage_query'] for hit in intraspecies_hits]
        inter_coverages = [hit['coverage_query'] for hit in interspecies_hits]
        
        # 计算平均值
        mean_identity = np.mean(all_identities)
        mean_coverage = np.mean(all_coverages)
        
        mean_intra_identity = np.mean(intra_identities) if intra_identities else 0
        mean_inter_identity = np.mean(inter_identities) if inter_identities else 0
        
        mean_intra_coverage = np.mean(intra_coverages) if intra_coverages else 0
        mean_inter_coverage = np.mean(inter_coverages) if inter_coverages else 0
        
        # 评分逻辑
        # 1. 种内一致性评分（如果有种内比对）
        if intra_identities:
            if mean_intra_identity >= 0.95:
                intra_score = 1.0
            elif mean_intra_identity >= 0.90:
                intra_score = 0.9
            elif mean_intra_identity >= 0.80:
                intra_score = 0.7
            else:
                intra_score = 0.5
        else:
            intra_score = 0.8  # 没有种内比对时的默认分数
        
        # 2. 种间一致性评分
        if inter_identities:
            if mean_inter_identity >= 0.80:
                inter_score = 1.0
            elif mean_inter_identity >= 0.70:
                inter_score = 0.8
            elif mean_inter_identity >= 0.60:
                inter_score = 0.6
            elif mean_inter_identity >= 0.50:
                inter_score = 0.4
            else:
                inter_score = 0.2
        else:
            inter_score = 0.3  # 没有种间比对的惩罚
        
        # 3. 覆盖度评分
        if mean_coverage >= 0.90:
            coverage_score = 1.0
        elif mean_coverage >= 0.80:
            coverage_score = 0.8
        elif mean_coverage >= 0.70:
            coverage_score = 0.6
        elif mean_coverage >= 0.50:
            coverage_score = 0.4
        else:
            coverage_score = 0.2
        
        # 4. 比对广度评分（与多少个不同物种有显著比对）
        unique_species = set(hit['hit_species'] for hit in significant_hits)
        total_species = len(set(seq.species_name for seq in all_sequences))
        species_coverage = len(unique_species) / max(total_species - 1, 1)  # 减1因为不包括自己
        
        if species_coverage >= 0.8:
            breadth_score = 1.0
        elif species_coverage >= 0.6:
            breadth_score = 0.8
        elif species_coverage >= 0.4:
            breadth_score = 0.6
        elif species_coverage >= 0.2:
            breadth_score = 0.4
        else:
            breadth_score = 0.2
        
        # 5. E-value质量评分
        mean_evalue = np.mean([hit['evalue'] for hit in significant_hits])
        if mean_evalue <= 1e-100:
            evalue_score = 1.0
        elif mean_evalue <= 1e-50:
            evalue_score = 0.9
        elif mean_evalue <= 1e-20:
            evalue_score = 0.7
        elif mean_evalue <= 1e-10:
            evalue_score = 0.5
        else:
            evalue_score = 0.3
        
        # 综合全局一致性评分
        global_consistency_score = (
            intra_score * 0.25 +        # 种内一致性
            inter_score * 0.30 +        # 种间一致性
            coverage_score * 0.20 +     # 覆盖度
            breadth_score * 0.15 +      # 比对广度
            evalue_score * 0.10         # E-value质量
        )
        
        # 确定一致性等级
        if global_consistency_score >= 0.8:
            consistency_rank = 'excellent'
        elif global_consistency_score >= 0.6:
            consistency_rank = 'good'
        elif global_consistency_score >= 0.4:
            consistency_rank = 'moderate'
        elif global_consistency_score >= 0.2:
            consistency_rank = 'poor'
        else:
            consistency_rank = 'very_poor'
        
        return {
            'global_consistency_score': global_consistency_score,
            'mean_identity': mean_identity,
            'mean_coverage': mean_coverage,
            'mean_intra_identity': mean_intra_identity,
            'mean_inter_identity': mean_inter_identity,
            'mean_intra_coverage': mean_intra_coverage,
            'mean_inter_coverage': mean_inter_coverage,
            'num_significant_hits': len(significant_hits),
            'intraspecies_hits': len(intraspecies_hits),
            'interspecies_hits': len(interspecies_hits),
            'species_coverage': species_coverage,
            'unique_hit_species': len(unique_species),
            'consistency_rank': consistency_rank,
            'mean_evalue': mean_evalue,
            'detailed_hits': significant_hits[:10]  # 保存前10个最佳比对
        }
    
    def analyze_phylogenetic_position(self, gene_name, blast_results, all_sequences):
        """分析序列在系统发育中的位置"""
        if gene_name not in blast_results:
            return {
                'phylo_position_score': 0.5,
                'is_outlier': False,
                'outlier_evidence': []
            }
        
        hits = blast_results[gene_name]
        significant_hits = [hit for hit in hits if hit['evalue'] <= self.blast_params['evalue']]
        
        if not significant_hits:
            return {
                'phylo_position_score': 0.3,
                'is_outlier': True,
                'outlier_evidence': ['no_significant_hits']
            }
        
        # 分析与其他物种的相似性模式
        species_similarities = defaultdict(list)
        for hit in significant_hits:
            species_similarities[hit['hit_species']].append(hit['identity'])
        
        # 计算每个物种的平均相似性
        species_mean_similarities = {
            species: np.mean(identities) 
            for species, identities in species_similarities.items()
        }
        
        outlier_evidence = []
        is_outlier = False
        
        # 检查异常模式
        # 1. 与某些物种异常高的相似性
        max_similarity = max(species_mean_similarities.values()) if species_mean_similarities else 0
        if max_similarity > 0.98:
            outlier_evidence.append('extremely_high_similarity')
        
        # 2. 与大多数物种的相似性都很低
        low_similarity_species = sum(1 for sim in species_mean_similarities.values() if sim < 0.6)
        if low_similarity_species > len(species_mean_similarities) * 0.7:
            outlier_evidence.append('broadly_low_similarity')
        
        # 3. 相似性分布异常（方差很大）
        if len(species_mean_similarities) > 2:
            similarities = list(species_mean_similarities.values())
            similarity_std = np.std(similarities)
            if similarity_std > 0.2:
                outlier_evidence.append('high_similarity_variance')
        
        # 4. 只与少数物种有显著比对
        total_species = len(set(seq.species_name for seq in all_sequences))
        hit_species_ratio = len(species_similarities) / max(total_species - 1, 1)
        if hit_species_ratio < 0.3:
            outlier_evidence.append('limited_species_hits')
        
        # 判断是否为异常序列
        is_outlier = len(outlier_evidence) >= 2
        
        # 计算系统发育位置评分
        if is_outlier:
            phylo_position_score = max(0.2, 0.8 - len(outlier_evidence) * 0.15)
        else:
            # 基于比对广度和一致性计算
            breadth_score = hit_species_ratio
            consistency_score = 1 - (np.std(list(species_mean_similarities.values())) 
                                   if len(species_mean_similarities) > 1 else 0)
            phylo_position_score = (breadth_score + consistency_score) / 2
        
        return {
            'phylo_position_score': phylo_position_score,
            'is_outlier': is_outlier,
            'outlier_evidence': outlier_evidence,
            'species_similarities': species_mean_similarities,
            'hit_species_ratio': hit_species_ratio,
            'max_similarity': max_similarity
        }
    
    def detect_global_sequence_anomalies(self, gene_name, blast_results, all_sequences):
        """基于全局BLAST结果检测序列异常"""
        anomalies = []
        
        if gene_name not in blast_results:
            anomalies.append("no_global_hits")
            return anomalies
        
        hits = blast_results[gene_name]
        significant_hits = [hit for hit in hits if hit['evalue'] <= self.blast_params['evalue']]
        
        if not significant_hits:
            anomalies.append("no_significant_global_hits")
            return anomalies
        
        # 1. 全局低相似性
        all_identities = [hit['identity'] for hit in significant_hits]
        if np.mean(all_identities) < 0.6:
            anomalies.append("global_low_similarity")
        
        # 2. 种内种间相似性倒挂
        intra_hits = [hit for hit in significant_hits if hit['is_intraspecies']]
        inter_hits = [hit for hit in significant_hits if not hit['is_intraspecies']]
        
        if intra_hits and inter_hits:
            mean_intra = np.mean([hit['identity'] for hit in intra_hits])
            mean_inter = np.mean([hit['identity'] for hit in inter_hits])
            if mean_inter > mean_intra:
                anomalies.append("intraspecies_lower_than_interspecies")
        
        # 3. 覆盖度异常
        coverages = [hit['coverage_query'] for hit in significant_hits]
        if np.mean(coverages) < 0.5:
            anomalies.append("global_low_coverage")
        
        # 4. 比对长度高度不一致
        align_lengths = [hit['align_length'] for hit in significant_hits]
        if len(set(align_lengths)) > len(align_lengths) * 0.8:
            anomalies.append("highly_variable_alignment_lengths")
        
        # 5. 与太少物种有比对
        unique_species = set(hit['hit_species'] for hit in significant_hits)
        total_species = len(set(seq.species_name for seq in all_sequences))
        if len(unique_species) < total_species * 0.3:
            anomalies.append("limited_species_representation")
        
        # 6. E-value分布异常
        evalues = [hit['evalue'] for hit in significant_hits]
        if np.mean(evalues) > 1e-5:
            anomalies.append("poor_evalue_distribution")
        
        return anomalies
    
    def calculate_basic_quality(self, sequence):
        """计算基础质量指标（复用之前的逻辑）"""
        seq_str = str(sequence.seq).upper()
        length = len(seq_str)
        
        # 长度质量
        if self.sequence_type == 'nucleotide':
            optimal_range = (300, 3000)
            min_acceptable = 100
            max_acceptable = 10000
        else:
            optimal_range = (100, 1000)
            min_acceptable = 50
            max_acceptable = 3000
        
        if length < min_acceptable:
            length_score = length / min_acceptable * 0.3
        elif length > max_acceptable:
            length_score = max_acceptable / length * 0.7
        elif optimal_range[0] <= length <= optimal_range[1]:
            length_score = 1.0
        elif length < optimal_range[0]:
            length_score = 0.3 + 0.7 * (length - min_acceptable) / (optimal_range[0] - min_acceptable)
        else:
            length_score = 1.0 - 0.3 * (length - optimal_range[1]) / (max_acceptable - optimal_range[1])
        
        length_score = max(0, min(1, length_score))
        
        # 组成质量
        if self.sequence_type == 'nucleotide':
            gc_content = gc_fraction(sequence.seq)
            optimal_gc_range = (0.35, 0.65)
            if optimal_gc_range[0] <= gc_content <= optimal_gc_range[1]:
                composition_score = 1.0
            else:
                if gc_content < optimal_gc_range[0]:
                    composition_score = max(0, gc_content / optimal_gc_range[0])
                else:
                    composition_score = max(0, (1 - gc_content) / (1 - optimal_gc_range[1]))
        else:
            composition_score = 0.8
        
        # 序列完整性
        if self.sequence_type == 'nucleotide':
            ambiguous_chars = set('NRYSWKMBDHV')
        else:
            ambiguous_chars = set('X')
        
        ambiguous_count = sum(seq_str.count(char) for char in ambiguous_chars)
        gap_count = seq_str.count('-')
        
        ambiguous_ratio = ambiguous_count / length if length > 0 else 0
        gap_ratio = gap_count / length if length > 0 else 0
        
        integrity_score = max(0, 1 - ambiguous_ratio * 10 - gap_ratio * 5)
        
        # 序列复杂度
        if length > 0:
            char_counts = Counter(seq_str)
            entropy = 0
            for count in char_counts.values():
                freq = count / length
                if freq > 0:
                    entropy -= freq * np.log2(freq)
            
            if self.sequence_type == 'nucleotide':
                max_entropy = np.log2(4)
            else:
                max_entropy = np.log2(20)
            
            complexity_score = entropy / max_entropy if max_entropy > 0 else 0
        else:
            complexity_score = 0
        
        # 编码潜力（仅核苷酸）
        if self.sequence_type == 'nucleotide':
            start_codons = ['ATG']
            stop_codons = ['TAA', 'TAG', 'TGA']
            
            max_orf = 0
            for frame in range(3):
                frame_seq = seq_str[frame:]
                current_orf = 0
                in_orf = False
                
                for i in range(0, len(frame_seq) - 2, 3):
                    codon = frame_seq[i:i+3]
                    if len(codon) != 3:
                        break
                    
                    if codon in start_codons and not in_orf:
                        in_orf = True
                        current_orf = 3
                    elif in_orf:
                        if codon in stop_codons:
                            max_orf = max(max_orf, current_orf)
                            in_orf = False
                            current_orf = 0
                        else:
                            current_orf += 3
            
            coding_score = min(1.0, max_orf / 300) if max_orf > 0 else 0.3
        else:
            coding_score = 1.0
        
        return {
            'length_score': length_score,
            'composition_score': composition_score,
            'integrity_score': integrity_score,
            'complexity_score': complexity_score,
            'coding_score': coding_score
        }
    
    def calculate_comprehensive_quality_with_global_blast(self, sequence, global_blast_results, all_sequences):
        """计算包含全局BLAST信息的综合质量评分"""
        # 基础质量分析
        base_quality = self.calculate_basic_quality(sequence)
        
        # 全局BLAST一致性分析
        global_blast_analysis = self.analyze_global_blast_consistency(
            sequence.gene_name, global_blast_results, all_sequences
        )
        
        # 系统发育位置分析
        phylo_analysis = self.analyze_phylogenetic_position(
            sequence.gene_name, global_blast_results, all_sequences
        )
        
        # 检测全局异常
        global_anomalies = self.detect_global_sequence_anomalies(
            sequence.gene_name, global_blast_results, all_sequences
        )
        
        # 计算加权总分
        total_score = (
            base_quality['length_score'] * self.weights['length_quality'] +
            base_quality['composition_score'] * self.weights['composition_quality'] +
            base_quality['integrity_score'] * self.weights['sequence_integrity'] +
            base_quality['complexity_score'] * self.weights['complexity'] +
            global_blast_analysis['global_consistency_score'] * self.weights['global_blast_consistency'] +
            phylo_analysis['phylo_position_score'] * self.weights['phylogenetic_position']
        )
        
        if self.sequence_type == 'nucleotide' and 'coding_potential' in self.weights:
            total_score += base_quality['coding_score'] * self.weights['coding_potential']
        
        # 全局异常惩罚
        global_anomaly_penalty = len(global_anomalies) * 0.08  # 每个异常扣8%
        total_score = max(0, total_score - global_anomaly_penalty)
        
        return {
            'total_score': total_score,
            'sequence_id': sequence.id,
            'gene_name': sequence.gene_name,
            'species_name': sequence.species_name,
            'sequence_length': len(sequence.seq),
            'base_quality': base_quality,
            'global_blast_analysis': global_blast_analysis,
            'phylo_analysis': phylo_analysis,
            'global_anomalies': global_anomalies,
            'global_anomaly_count': len(global_anomalies)
        }
    
    def run_global_blast_enhanced_analysis(self, input_fasta, output_prefix="global_blast"):
        """运行完整的全局BLAST增强质量分析"""
        logger.info(f"Starting global BLAST-enhanced quality analysis: {input_fasta}")
        
        try:
            # 解析所有序列
            all_sequences = []
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
                
                all_sequences.append(record)
                species_groups[species_name].append(record)
            
            self.all_sequences = all_sequences
            self.species_groups = dict(species_groups)
            
            logger.info(f"Found {len(all_sequences)} sequences from {len(self.species_groups)} species")
            logger.info(f"Gene family analysis mode: all sequences will be compared globally")
            
            # 执行全局BLAST分析
            with tempfile.TemporaryDirectory() as temp_dir:
                global_blast_results = self.perform_global_blast_analysis(all_sequences, temp_dir)
                
                if not global_blast_results:
                    logger.warning("Global BLAST analysis failed, using basic quality metrics only")
                    global_blast_results = {}
                
                # 分析每个序列
                all_results = []
                for seq in all_sequences:
                    analysis = self.calculate_comprehensive_quality_with_global_blast(
                        seq, global_blast_results, all_sequences
                    )
                    all_results.append(analysis)
                
                # 按总分排序
                all_results.sort(key=lambda x: x['total_score'], reverse=True)
                
                # 添加全局排名信息
                for i, result in enumerate(all_results):
                    result['global_rank'] = i + 1
                    result['global_percentile'] = (len(all_results) - i) / len(all_results)
                
                # 按物种重新组织结果
                species_results = defaultdict(list)
                for result in all_results:
                    species_results[result['species_name']].append(result)
                
                # 为每个物种内部添加排名
                for species, results in species_results.items():
                    results.sort(key=lambda x: x['total_score'], reverse=True)
                    for i, result in enumerate(results):
                        result['species_rank'] = i + 1
                
                # 生成报告
                self.generate_global_blast_report(species_results, all_results, f"{output_prefix}_global_report.txt")
                
                # 选择每个物种的最佳序列
                selected_sequences = []
                removed_sequences = []
                
                for species, results in species_results.items():
                    best_result = results[0]  # 物种内最佳
                    
                    # 找到对应的序列对象
                    for seq in all_sequences:
                        if (seq.gene_name == best_result['gene_name'] and 
                            seq.species_name == best_result['species_name']):
                            seq.description = f"{seq.gene_name} {seq.species_name}"
                            seq.id = seq.gene_name
                            selected_sequences.append(seq)
                            break
                    
                    # 收集被移除的序列
                    for result in results[1:]:
                        for seq in all_sequences:
                            if (seq.gene_name == result['gene_name'] and 
                                seq.species_name == result['species_name']):
                                seq.description = f"{seq.gene_name} {seq.species_name}"
                                seq.id = seq.gene_name
                                removed_sequences.append(seq)
                                break
                
                # 保存结果
                SeqIO.write(selected_sequences, f"{output_prefix}_best_sequences.fasta", "fasta")
                if removed_sequences:
                    SeqIO.write(removed_sequences, f"{output_prefix}_removed_sequences.fasta", "fasta")
                
                # 统计信息
                high_quality_sequences = sum(1 for r in all_results if r['total_score'] >= 0.7)
                global_outliers = sum(1 for r in all_results if r['phylo_analysis']['is_outlier'])
                
                logger.info("=" * 60)
                logger.info("GLOBAL BLAST-ENHANCED ANALYSIS COMPLETED")
                logger.info("=" * 60)
                logger.info(f"Total sequences analyzed: {len(all_sequences)}")
                logger.info(f"Total species: {len(self.species_groups)}")
                logger.info(f"High-quality sequences: {high_quality_sequences}")
                logger.info(f"Phylogenetic outliers detected: {global_outliers}")
                logger.info(f"Sequences selected (1 per species): {len(selected_sequences)}")
                logger.info(f"Sequences removed: {len(removed_sequences)}")
                logger.info("=" * 60)
                
                return {
                    'total_sequences': len(all_sequences),
                    'selected_count': len(selected_sequences),
                    'removed_count': len(removed_sequences),
                    'high_quality_count': high_quality_sequences,
                    'outlier_count': global_outliers
                }
                
        except Exception as e:
            logger.error(f"Global BLAST-enhanced analysis failed: {e}")
            raise
    
    def generate_global_blast_report(self, species_results, all_results, output_file):
        """生成全局BLAST分析报告"""
        try:
            with open(output_file, 'w') as f:
                f.write("=" * 120 + "\n")
                f.write("GLOBAL BLAST-ENHANCED GENE FAMILY QUALITY ANALYSIS REPORT\n")
                f.write("=" * 120 + "\n")
                f.write(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S UTC')}\n")
                f.write(f"User: Biols9527\n")
                f.write(f"Analysis type: Gene family global comparison\n")
                f.write(f"Sequence type: {self.sequence_type}\n")
                f.write("\n")
                
                # 分析配置
                f.write("ANALYSIS CONFIGURATION:\n")
                f.write("-" * 60 + "\n")
                f.write("Quality weights:\n")
                for metric, weight in self.weights.items():
                    f.write(f"  {metric}: {weight:.3f}\n")
                f.write("\nGlobal BLAST parameters:\n")
                for param, value in self.blast_params.items():
                    f.write(f"  {param}: {value}\n")
                f.write("\n")
                
                # 全局统计
                f.write("GLOBAL GENE FAMILY STATISTICS:\n")
                f.write("-" * 60 + "\n")
                f.write(f"Total sequences in gene family: {len(all_results)}\n")
                f.write(f"Total species represented: {len(species_results)}\n")
                
                global_scores = [r['total_score'] for r in all_results]
                blast_scores = [r['global_blast_analysis']['global_consistency_score'] for r in all_results]
                phylo_scores = [r['phylo_analysis']['phylo_position_score'] for r in all_results]
                
                f.write(f"Mean total quality score: {np.mean(global_scores):.3f} ± {np.std(global_scores):.3f}\n")
                f.write(f"Mean global BLAST consistency: {np.mean(blast_scores):.3f} ± {np.std(blast_scores):.3f}\n")
                f.write(f"Mean phylogenetic position score: {np.mean(phylo_scores):.3f} ± {np.std(phylo_scores):.3f}\n")
                
                # 质量分布
                high_quality = sum(1 for score in global_scores if score >= 0.7)
                medium_quality = sum(1 for score in global_scores if 0.4 <= score < 0.7)
                low_quality = sum(1 for score in global_scores if score < 0.4)
                
                f.write(f"\nQuality distribution:\n")
                f.write(f"  High quality (≥0.7): {high_quality} ({high_quality/len(all_results)*100:.1f}%)\n")
                f.write(f"  Medium quality (0.4-0.7): {medium_quality} ({medium_quality/len(all_results)*100:.1f}%)\n")
                f.write(f"  Low quality (<0.4): {low_quality} ({low_quality/len(all_results)*100:.1f}%)\n")
                
                # 异常检测统计
                outliers = sum(1 for r in all_results if r['phylo_analysis']['is_outlier'])
                anomalous = sum(1 for r in all_results if r['global_anomaly_count'] > 2)
                
                f.write(f"\nAnomaly detection:\n")
                f.write(f"  Phylogenetic outliers: {outliers} ({outliers/len(all_results)*100:.1f}%)\n")
                f.write(f"  Sequences with multiple anomalies: {anomalous} ({anomalous/len(all_results)*100:.1f}%)\n")
                f.write("\n")
                
                # 全局排名（前10和后10）
                f.write("GLOBAL RANKING:\n")
                f.write("-" * 60 + "\n")
                f.write("Top 10 highest quality sequences:\n")
                f.write(f"{'Rank':<4} {'Gene':<15} {'Species':<20} {'Score':<7} {'BLAST':<7} {'Phylo':<7} {'Issues':<10}\n")
                f.write("-" * 80 + "\n")
                
                for i, result in enumerate(all_results[:10]):
                    issues = len(result['global_anomalies'])
                    f.write(f"{result['global_rank']:<4} "
                           f"{result['gene_name']:<15} "
                           f"{result['species_name']:<20} "
                           f"{result['total_score']:<7.3f} "
                           f"{result['global_blast_analysis']['global_consistency_score']:<7.3f} "
                           f"{result['phylo_analysis']['phylo_position_score']:<7.3f} "
                           f"{issues if issues > 0 else '-':<10}\n")
                
                if len(all_results) > 10:
                    f.write("\nBottom 10 lowest quality sequences:\n")
                    f.write(f"{'Rank':<4} {'Gene':<15} {'Species':<20} {'Score':<7} {'BLAST':<7} {'Phylo':<7} {'Issues':<10}\n")
                    f.write("-" * 80 + "\n")
                    
                    for result in all_results[-10:]:
                        issues = len(result['global_anomalies'])
                        f.write(f"{result['global_rank']:<4} "
                               f"{result['gene_name']:<15} "
                               f"{result['species_name']:<20} "
                               f"{result['total_score']:<7.3f} "
                               f"{result['global_blast_analysis']['global_consistency_score']:<7.3f} "
                               f"{result['phylo_analysis']['phylo_position_score']:<7.3f} "
                               f"{issues if issues > 0 else '-':<10}\n")
                
                f.write("\n")
                
                # 各物种详细分析
                f.write("SPECIES-WISE DETAILED ANALYSIS:\n")
                f.write("=" * 120 + "\n")
                
                for species, results in species_results.items():
                    f.write(f"\nSpecies: {species}\n")
                    f.write("-" * 90 + "\n")
                    f.write(f"Number of sequences: {len(results)}\n")
                    
                    if len(results) == 1:
                        result = results[0]
                        f.write(f"Single representative: {result['gene_name']}\n")
                        f.write(f"Global rank: {result['global_rank']}/{len(all_results)}\n")
                        f.write(f"Quality score: {result['total_score']:.3f}\n")
                        f.write(f"Global BLAST consistency: {result['global_blast_analysis']['global_consistency_score']:.3f}\n")
                        f.write(f"Phylogenetic position: {result['phylo_analysis']['phylo_position_score']:.3f}\n")
                        
                        if result['global_anomalies']:
                            f.write(f"Detected issues: {', '.join(result['global_anomalies'])}\n")
                        else:
                            f.write("No issues detected\n")
                    else:
                        # 多序列比较
                        f.write(f"Quality scores range: {results[-1]['total_score']:.3f} - {results[0]['total_score']:.3f}\n")
                        f.write(f"Global ranks range: {results[-1]['global_rank']} - {results[0]['global_rank']}\n")
                        
                        f.write("\nSequence comparison within species:\n")
                        f.write(f"{'SRank':<5} {'Gene':<15} {'GRank':<6} {'Total':<7} {'BLAST':<7} {'Phylo':<7} {'Issues':<15}\n")
                        f.write("-" * 70 + "\n")
                        
                        for result in results:
                            issues_str = ','.join(result['global_anomalies'][:2]) if result['global_anomalies'] else '-'
                            if len(result['global_anomalies']) > 2:
                                issues_str += f"+{len(result['global_anomalies'])-2}"
                            
                            f.write(f"{result['species_rank']:<5} "
                                   f"{result['gene_name']:<15} "
                                   f"{result['global_rank']:<6} "
                                   f"{result['total_score']:<7.3f} "
                                   f"{result['global_blast_analysis']['global_consistency_score']:<7.3f} "
                                   f"{result['phylo_analysis']['phylo_position_score']:<7.3f} "
                                   f"{issues_str:<15}\n")
                        
                        # 推荐的最佳序列
                        best = results[0]
                        f.write(f"\nRECOMMENDED: {best['gene_name']} (Global rank: {best['global_rank']})\n")
                        f.write(f"Reasoning:\n")
                        f.write(f"  Total score: {best['total_score']:.3f} (top within species)\n")
                        f.write(f"  Global BLAST hits: {best['global_blast_analysis']['num_significant_hits']}\n")
                        f.write(f"  Inter-species similarity: {best['global_blast_analysis']['mean_inter_identity']:.3f}\n")
                        f.write(f"  Species coverage: {best['global_blast_analysis']['species_coverage']:.1%}\n")
                        
                        if best['phylo_analysis']['is_outlier']:
                            f.write(f"  ⚠️  Potential outlier: {', '.join(best['phylo_analysis']['outlier_evidence'])}\n")
                        
                        if best['global_anomalies']:
                            f.write(f"  Issues to note: {', '.join(best['global_anomalies'])}\n")
                    
                    f.write("\n" + "=" * 90 + "\n")
                
                # 方法说明
                f.write("\nMETHODOLOGY:\n")
                f.write("-" * 60 + "\n")
                f.write("GLOBAL GENE FAMILY ANALYSIS:\n")
                f.write("This analysis treats all input sequences as members of the same gene family\n")
                f.write("and performs comprehensive global comparisons.\n\n")
                f.write("Key features:\n")
                f.write("1. Global BLAST analysis:\n")
                f.write("   - All-vs-all BLAST comparison within the gene family\n")
                f.write("   - Intra-species vs inter-species similarity assessment\n")
                f.write("   - Species coverage and consistency evaluation\n\n")
                f.write("2. Phylogenetic position analysis:\n")
                f.write("   - Detection of sequences with unusual similarity patterns\n")
                f.write("   - Identification of potential outliers or contaminants\n")
                f.write("   - Assessment of sequence relationships across species\n\n")
                f.write("3. Multi-dimensional quality scoring:\n")
                f.write("   - Traditional sequence quality metrics\n")
                f.write("   - Global consistency within gene family\n")
                f.write("   - Phylogenetic appropriateness\n\n")
                f.write("Selection strategy:\n")
                f.write("- For each species, select the highest-scoring representative\n")
                f.write("- Prioritize sequences with good global consistency\n")
                f.write("- Avoid phylogenetic outliers when possible\n")
                
            logger.info(f"Global BLAST analysis report saved to {output_file}")
            
        except Exception as e:
            logger.error(f"Failed to generate global BLAST report: {e}")

def main():
    """主函数"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Global BLAST-enhanced gene family quality analyzer")
    parser.add_argument("input_fasta", help="Input FASTA file (gene family sequences)")
    parser.add_argument("--sequence_type", choices=['nucleotide', 'protein'], 
                       default='nucleotide', help="Sequence type")
    parser.add_argument("--output_prefix", default="global_blast", help="Output prefix")
    parser.add_argument("--evalue", type=float, default=1e-10, help="BLAST E-value threshold")
    parser.add_argument("--num_threads", type=int, default=8, help="Number of BLAST threads")
    
    args = parser.parse_args()
    
    # 配置BLAST参数
    blast_params = {
        'evalue': args.evalue,
        'word_size': 11 if args.sequence_type == 'nucleotide' else 3,
        'max_target_seqs': 1000,
        'outfmt': 6,
        'num_threads': args.num_threads
    }
    
    # 初始化分析器
    analyzer = GlobalBlastQualityAnalyzer(
        sequence_type=args.sequence_type,
        blast_params=blast_params
    )
    
    # 运行分析
    result = analyzer.run_global_blast_enhanced_analysis(args.input_fasta, args.output_prefix)
    
    print("\n" + "="*80)
    print("✅ GLOBAL GENE FAMILY ANALYSIS COMPLETED!")
    print("="*80)
    print(f"ð Best sequences (1 per species): {args.output_prefix}_best_sequences.fasta")
    if result['removed_count'] > 0:
        print(f"ð️  Alternative sequences: {args.output_prefix}_removed_sequences.fasta")
    print(f"ð Comprehensive report: {args.output_prefix}_global_report.txt")
    print(f"ð Summary: {result['selected_count']} selected from {result['total_sequences']} total")
    print(f"⭐ High quality: {result['high_quality_count']} sequences")
    if result['outlier_count'] > 0:
        print(f"⚠️  Outliers detected: {result['outlier_count']} sequences")
    print("="*80)

if __name__ == "__main__":
    main()
