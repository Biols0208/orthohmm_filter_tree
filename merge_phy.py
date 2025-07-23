#!/dellfsqd2/ST_OCEAN/USER/lishuo11/01_soft/mambaforge/bin/python
import sys
import os
import logging
import glob
import argparse
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
import time

# Global logger - will be configured later
logger = logging.getLogger(__name__)

# ==================================================
# Format Conversion Functions
# ==================================================
def phy_to_fasta(sequences_dict: Dict[str, str]) -> str:
    """
    将PHYLIP格式的序列字典转换为FASTA格式字符串
    
    Args:
        sequences_dict: 物种名到序列的映射字典
        
    Returns:
        str: FASTA格式的字符串
    """
    fasta_lines = []
    for species_name, sequence in sequences_dict.items():
        fasta_lines.append(f">{species_name}")
        fasta_lines.append(sequence)
    return "\n".join(fasta_lines)

def fasta_to_phy_dict(fasta_content: str) -> Tuple[Dict[str, str], List[str]]:
    """
    将FASTA格式内容解析为PHYLIP格式所需的数据结构
    
    Args:
        fasta_content: FASTA格式的文件内容
        
    Returns:
        tuple: (序列字典, 序列ID列表)
    """
    sequence_dict = {}
    sequence_list = []
    current_id = None
    seq_parts = []
    
    for line in fasta_content.split('\n'):
        line = line.rstrip()
        if not line:
            continue
            
        if line.startswith('>'):
            # 处理前一个序列
            if current_id:
                sequence_dict[current_id] = ''.join(seq_parts)
            
            # 提取新序列ID
            header = line[1:].strip()
            current_id = header.split()[0]
            sequence_list.append(current_id)
            seq_parts = []
        else:
            seq_parts.append(line)
    
    # 处理最后一个序列
    if current_id and seq_parts:
        sequence_dict[current_id] = ''.join(seq_parts)
    
    return sequence_dict, sequence_list

def write_fasta_file(sequences_dict: Dict[str, str], output_file: str) -> None:
    """
    将序列字典写入FASTA格式文件
    
    Args:
        sequences_dict: 物种名到序列的映射字典
        output_file: 输出文件路径
    """
    logger.info(f"写入FASTA文件: {output_file}")
    try:
        with open(output_file, 'w') as f:
            fasta_content = phy_to_fasta(sequences_dict)
            f.write(fasta_content)
            if not fasta_content.endswith('\n'):
                f.write('\n')
        logger.info(f"成功写入FASTA文件: {output_file}")
    except IOError as e:
        logger.error(f"写入FASTA文件错误 {output_file}: {e}")
        raise

def write_phylip_from_fasta(fasta_file: str, phylip_file: str) -> None:
    """
    从FASTA文件转换并写入PHYLIP格式文件
    
    Args:
        fasta_file: 输入FASTA文件路径
        phylip_file: 输出PHYLIP文件路径
    """
    logger.info(f"将FASTA文件 {fasta_file} 转换为PHYLIP格式 {phylip_file}")
    
    try:
        with open(fasta_file, 'r') as f:
            fasta_content = f.read()
        
        sequence_dict, sequence_list = fasta_to_phy_dict(fasta_content)
        
        # 验证序列长度一致性
        alignment_length = 0
        for gene in sequence_dict:
            if alignment_length == 0:
                alignment_length = len(sequence_dict[gene])
            elif len(sequence_dict[gene]) != alignment_length:
                raise ValueError(
                    f"序列长度不一致: {gene}序列长度({len(sequence_dict[gene])})与其他序列长度({alignment_length})不匹配"
                )
        
        if not sequence_list:
            raise ValueError("没有找到有效的序列")
        
        # 写入PHYLIP文件
        longest_id_len = max(len(id) for id in sequence_list)
        
        with open(phylip_file, "w") as phyfile:
            phyfile.write(f"{len(sequence_dict)} {alignment_length}\n")
            for gene in sequence_list:
                phyfile.write(f"{gene.ljust(longest_id_len)}   {sequence_dict[gene]}\n")
        
        logger.info(f"成功转换: {len(sequence_dict)}个序列, 长度{alignment_length}")
        
    except Exception as e:
        logger.error(f"FASTA到PHYLIP转换错误: {e}")
        raise

# ==================================================
# Configuration and Validation Functions
# ==================================================
def setup_logging(verbose: int = 0, quiet: bool = False) -> None:
    """
    Configure logging based on verbosity level.
    
    Args:
        verbose: Verbosity level (0=INFO, 1=DEBUG, 2+=DEBUG with more detail)
        quiet: If True, suppress all output except errors
    """
    if quiet:
        level = logging.ERROR
    elif verbose == 0:
        level = logging.INFO
    else:
        level = logging.DEBUG
    
    # Configure the logger
    logging.basicConfig(
        level=level,
        format='%(levelname)s: %(message)s',
        force=True  # Override any existing configuration
    )
    
    # Set the level for our specific logger
    logger.setLevel(level)

def validate_sequence(sequence: str, sequence_type: str = "auto") -> bool:
    """
    Validate if a sequence is a valid DNA/RNA/Protein sequence.
    
    Args:
        sequence: The sequence to validate
        sequence_type: Type of sequence ("dna", "rna", "protein", "auto")
    
    Returns:
        bool: True if sequence is valid, False otherwise
    """
    if not sequence or sequence.isspace():
        return False
    
    # Remove gaps and unknown characters for validation
    clean_seq = re.sub(r'[-\?NnXx]', '', sequence.upper())
    
    if not clean_seq:  # All gaps/unknowns
        return True
    
    if sequence_type == "auto":
        # Auto-detect sequence type
        dna_chars = set('ATCG')
        rna_chars = set('AUCG')
        protein_chars = set('ACDEFGHIKLMNPQRSTVWY')
        
        seq_chars = set(clean_seq)
        
        if seq_chars.issubset(dna_chars):
            sequence_type = "dna"
        elif seq_chars.issubset(rna_chars):
            sequence_type = "rna"
        elif seq_chars.issubset(protein_chars):
            sequence_type = "protein"
        else:
            logger.warning(f"Sequence contains unexpected characters: {seq_chars - protein_chars}")
            return False
    
    # Validate based on detected/specified type
    valid_chars = {
        "dna": set('ATCG'),
        "rna": set('AUCG'), 
        "protein": set('ACDEFGHIKLMNPQRSTVWY')
    }
    
    if sequence_type in valid_chars:
        return set(clean_seq).issubset(valid_chars[sequence_type])
    
    return False

def print_statistics(final_sequences: Dict[str, str], file_count: int, processing_time: float) -> None:
    """
    Print comprehensive statistics about the merge operation.
    
    Args:
        final_sequences: Dictionary of final merged sequences
        file_count: Number of input files processed
        processing_time: Time taken for processing in seconds
    """
    if not final_sequences:
        print("\n=== MERGE STATISTICS ===")
        print("No sequences were merged.")
        return
    
    print("\n=== MERGE STATISTICS ===")
    print(f"Files processed: {file_count}")
    print(f"Total genes/species: {len(final_sequences)}")
    print(f"Final sequence length: {len(next(iter(final_sequences.values())))}")
    print(f"Processing time: {processing_time:.2f} seconds")
    
    # Sequence composition analysis
    if final_sequences:
        sample_seq = next(iter(final_sequences.values()))
        gap_count = sample_seq.count('-')
        total_length = len(sample_seq)
        gap_percentage = (gap_count / total_length) * 100 if total_length > 0 else 0
        
        print(f"Gap percentage: {gap_percentage:.1f}%")
        
        # Gene name analysis
        gene_names = list(final_sequences.keys())
        print(f"Gene/Species names: {', '.join(gene_names[:5])}")
        if len(gene_names) > 5:
            print(f"  ... and {len(gene_names) - 5} more")
    
    print("=" * 25)

# ==================================================
# Function to Read and Validate a Single PHYLIP File
# ==================================================
def read_phy_file_stream(phy_file: str, validate_sequences: bool = False, sequence_type: str = "auto"):
    """
    Generator version of read_phy_file for memory efficiency.
    Yields sequences one by one instead of loading all into memory.
    
    Args:
        phy_file: Path to the PHYLIP file
        validate_sequences: Whether to validate sequence content  
        sequence_type: Type of sequences to validate
        
    Yields:
        tuple: (gene_name, sequence) pairs
        
    Returns:
        int: Expected sequence length from header (returned after StopIteration)
    """
    expected_num_seqs = None
    expected_seq_length = None
    sequences_read = 0

    logger.info(f"Reading file: {phy_file}")

    try:
        with open(phy_file, 'r') as infile:
            # Read and parse header
            try:
                header = infile.readline().strip()
                if not header:
                    raise ValueError(f"File {phy_file} is empty or header is missing.")
                parts = header.split()
                if len(parts) != 2:
                    raise ValueError(f"Invalid header format in {phy_file}: '{header}'. Expected 'num_sequences sequence_length'.")
                expected_num_seqs = int(parts[0])
                expected_seq_length = int(parts[1])
                logger.debug(f"  Header indicates {expected_num_seqs} sequences, {expected_seq_length} length.")
            except ValueError as e:
                raise ValueError(f"Error parsing header in {phy_file}: {e}") from e

            if expected_num_seqs <= 0 or expected_seq_length <= 0:
                raise ValueError(f"Header values (sequences, length) in {phy_file} must be positive integers.")

            # Read sequences one by one
            seen_genes = set()
            for line_num, line in enumerate(infile, 2):
                line = line.strip()
                if not line:  # Skip empty lines
                    continue

                parts = line.split(maxsplit=1)
                if len(parts) != 2:
                    logger.warning(f"  Skipping malformed line {line_num} in {phy_file}: '{line}'. Expected 'name sequence'.")
                    continue

                gene_name = parts[0]
                sequence = parts[1].replace(" ", "")

                # Validate sequence length
                if len(sequence) != expected_seq_length:
                    raise ValueError(
                        f"Sequence length mismatch in {phy_file}, line {line_num} for gene '{gene_name}'. "
                        f"Expected {expected_seq_length}, got {len(sequence)}."
                    )

                # Validate sequence content (if requested)
                if validate_sequences and not validate_sequence(sequence, sequence_type):
                    logger.warning(
                        f"Invalid sequence content for gene '{gene_name}' in {phy_file}, line {line_num}. "
                        f"Sequence may contain unexpected characters."
                    )

                # Handle duplicates
                if gene_name in seen_genes:
                    logger.warning(f"  Duplicate gene name '{gene_name}' found in {phy_file}. Overwriting previous entry.")
                seen_genes.add(gene_name)
                
                sequences_read += 1
                yield gene_name, sequence

            # Final validation: sequence count
            if sequences_read != expected_num_seqs:
                raise ValueError(
                    f"Number of sequences read ({sequences_read}) in {phy_file} "
                    f"does not match header value ({expected_num_seqs})."
                )

    except FileNotFoundError:
        logger.error(f"Error: Input file not found: {phy_file}")
        raise
    except Exception as e:
        logger.error(f"An error occurred while processing file {phy_file}: {e}")
        raise

    logger.debug(f"  Successfully read {sequences_read} sequences from {phy_file}.")
    return expected_seq_length

def read_phy_file(phy_file: str, validate_sequences: bool = False, sequence_type: str = "auto") -> Tuple[Dict[str, str], int]:
    """
    Reads a single PHYLIP (.phy) file, validates its format, and extracts
    gene names and sequences.

    Args:
        phy_file (str): The path to the PHYLIP file.
        validate_sequences (bool): Whether to validate sequence content
        sequence_type (str): Type of sequences to validate ("dna", "rna", "protein", "auto")

    Returns:
        tuple: A tuple containing:
            - dict: A dictionary mapping gene names (str) to sequences (str).
            - int: The expected sequence length for this file based on its header.

    Raises:
        ValueError: If the file format is invalid (header issues, sequence length
                    mismatches, incorrect sequence count).
        FileNotFoundError: If the input file cannot be found.
        Exception: For other potential I/O errors.
    """
    gene_sequences = {}
    expected_num_seqs = None
    expected_seq_length = None
    sequences_read = 0

    logger.info(f"Reading file: {phy_file}")

    try:
        with open(phy_file, 'r') as infile:
            # --- Read and Parse Header ---
            try:
                header = infile.readline().strip()
                if not header:
                    raise ValueError(f"File {phy_file} is empty or header is missing.")
                parts = header.split()
                if len(parts) != 2:
                    raise ValueError(f"Invalid header format in {phy_file}: '{header}'. Expected 'num_sequences sequence_length'.")
                expected_num_seqs = int(parts[0])
                expected_seq_length = int(parts[1])
                logger.debug(f"  Header indicates {expected_num_seqs} sequences, {expected_seq_length} length.")
            except ValueError as e:
                # Catch issues with int() conversion or wrong number of parts
                raise ValueError(f"Error parsing header in {phy_file}: {e}") from e
            except IndexError:
                 # Should be caught by len(parts) check, but included for safety
                 raise ValueError(f"Invalid header format in {phy_file}: '{header}'. Not enough parts.")

            if expected_num_seqs <= 0 or expected_seq_length <= 0:
                 raise ValueError(f"Header values (sequences, length) in {phy_file} must be positive integers.")

            # --- Read Sequence Data ---
            for line_num, line in enumerate(infile, 2): # Start line count from 2 for content
                line = line.strip()
                if not line:  # Skip empty lines
                    continue

                # Split only on the first whitespace to separate name from sequence
                parts = line.split(maxsplit=1)
                if len(parts) != 2:
                    logger.warning(f"  Skipping malformed line {line_num} in {phy_file}: '{line}'. Expected 'name sequence'.")
                    continue

                gene_name = parts[0]
                # Remove potential spaces within the sequence part itself
                sequence = parts[1].replace(" ", "")

                # --- Validate Sequence Length ---
                if len(sequence) != expected_seq_length:
                    raise ValueError(
                        f"Sequence length mismatch in {phy_file}, line {line_num} for gene '{gene_name}'. "
                        f"Expected {expected_seq_length}, got {len(sequence)}."
                    )

                # --- Validate Sequence Content (if requested) ---
                if validate_sequences and not validate_sequence(sequence, sequence_type):
                    logger.warning(
                        f"Invalid sequence content for gene '{gene_name}' in {phy_file}, line {line_num}. "
                        f"Sequence may contain unexpected characters."
                    )

                # --- Store Sequence (Handle Duplicates) ---
                if gene_name in gene_sequences:
                     logger.warning(f"  Duplicate gene name '{gene_name}' found in {phy_file}. Overwriting previous entry.")
                gene_sequences[gene_name] = sequence
                sequences_read += 1

            # --- Final Validation: Sequence Count ---
            if sequences_read != expected_num_seqs:
                raise ValueError(
                    f"Number of sequences read ({sequences_read}) in {phy_file} "
                    f"does not match header value ({expected_num_seqs})."
                )

    except FileNotFoundError:
        logger.error(f"Error: Input file not found: {phy_file}")
        raise # Re-raise the specific error
    except Exception as e:
        logger.error(f"An error occurred while processing file {phy_file}: {e}")
        raise # Re-raise after logging

    logger.debug(f"  Successfully read {sequences_read} sequences from {phy_file}.")
    return gene_sequences, expected_seq_length

# ==================================================
# Function to Merge Multiple PHYLIP Files
# ==================================================
def merge_phy_files(phy_files: List[str], gap_char: str = "-", continue_on_error: bool = False, 
                   validate_sequences: bool = False, sequence_type: str = "auto") -> Tuple[Dict[str, str], int, int]:
    """
    Merges sequence data from multiple PHYLIP files.

    Aligns sequences by gene name across files. If a gene is missing in a
    file, it's padded with the specified gap character repeated for the
    length of sequences in that specific file.

    Args:
        phy_files (list): A list of paths to the input PHYLIP files. The order
                          in this list determines the concatenation order.
        gap_char (str, optional): The character to use for padding missing
                                  sequences. Defaults to "-".

    Returns:
        tuple: A tuple containing:
            - dict: A dictionary mapping gene names (str) to the final merged
                    and padded sequences (str).
            - int: The total number of unique gene names found across all files.
            - int: The total length of the concatenated sequences.

    Raises:
        SystemExit: If any input file fails to read or process.
    """
    # Structure: {gene_name: [seq_file0, seq_file1, None, seq_file3, ...]}
    merged_sequences_fragments = {}
    all_genes = set()
    file_lengths = [] # Stores sequence length for each file index

    # --- Stage 1: Read all files, collect data ---
    logger.info("Starting merge process...")
    # Sort files for consistent concatenation order regardless of glob expansion order
    phy_files.sort()
    logger.info(f"Processing {len(phy_files)} input files")

    for i, phy_file in enumerate(phy_files):
        try:
            gene_sequences, seq_length = read_phy_file(phy_file, validate_sequences, sequence_type)
            file_lengths.append(seq_length)
            current_file_genes = set(gene_sequences.keys())
            all_genes.update(current_file_genes)

            # Add sequences from this file to the fragment dictionary
            for gene_name, sequence in gene_sequences.items():
                if gene_name not in merged_sequences_fragments:
                    # Initialize list with placeholders if gene is new
                    merged_sequences_fragments[gene_name] = [None] * len(phy_files)
                # Place the sequence at the correct index for this file
                merged_sequences_fragments[gene_name][i] = sequence

        except (ValueError, FileNotFoundError, Exception) as e:
            # Handle errors based on continue_on_error flag
            if continue_on_error:
                logger.warning(f"Skipping file {phy_file} due to error: {e}")
                file_lengths.append(0)  # Placeholder for skipped file
                continue
            else:
                logger.error(f"Failed to read or process file {phy_file}. Aborting merge. Error: {e}")
                sys.exit(1) # Exit the script

    num_genes = len(all_genes)
    if not num_genes:
        logging.warning("No genes found across any input files.")
        return {}, 0, 0 # Return empty results if no genes

    total_merged_length = sum(file_lengths)
    logger.info(f"Total unique genes found: {num_genes}")
    logger.debug(f"Individual file lengths: {file_lengths}")
    logger.info(f"Calculated total merged sequence length: {total_merged_length}")

    # --- Stage 2: Pad missing sequences with gaps ---
    logger.debug("Padding missing sequences with gaps...")
    for gene_name in all_genes:
        # Ensure every gene has an entry (should exist from Stage 1 if all_genes is correct)
        if gene_name not in merged_sequences_fragments:
             logger.warning(f"Gene '{gene_name}' in all_genes set but missing from fragments dict. Initializing.")
             merged_sequences_fragments[gene_name] = [None] * len(phy_files)

        # Iterate through the file indices for this gene's fragments
        for file_idx in range(len(phy_files)):
            if merged_sequences_fragments[gene_name][file_idx] is None:
                # If fragment is missing (None), create gap sequence of correct length
                gap_sequence = gap_char * file_lengths[file_idx]
                merged_sequences_fragments[gene_name][file_idx] = gap_sequence
                # logging.debug(f"  Padded '{gene_name}' for file index {file_idx} (length {file_lengths[file_idx]})")

    # --- Stage 3: Concatenate fragments into final sequences ---
    final_sequences = {}
    logger.debug("Concatenating sequence fragments...")
    for gene_name, sequence_parts in merged_sequences_fragments.items():
        final_sequence = "".join(sequence_parts)
        # Final integrity check on length
        if len(final_sequence) != total_merged_length:
             logger.error(f"Internal Error: Final sequence length for '{gene_name}' is {len(final_sequence)}, but expected {total_merged_length}. Aborting.")
             sys.exit(1) # Critical error if lengths don't match
        final_sequences[gene_name] = final_sequence

    logger.info("Sequence padding and concatenation complete.")
    return final_sequences, num_genes, total_merged_length

# ==================================================
# Function to Write Merged Data to PHYLIP File
# ==================================================
def write_merged_phy_stream(gene_iterator, output_file_path: str, num_genes: int, total_merged_length: int) -> None:
    """
    Memory-efficient version that writes sequences as they are generated.
    
    Args:
        gene_iterator: Iterator yielding (gene_name, sequence) tuples  
        output_file_path: Path for the output PHYLIP file
        num_genes: Total number of unique genes
        total_merged_length: Total length of merged sequences
    """
    logger.info(f"Writing merged data to {output_file_path}...")
    try:
        with open(output_file_path, 'w') as outfile:
            # Write header
            logger.debug(f"  Writing header: {num_genes} sequences, {total_merged_length} total length.")
            outfile.write(f" {num_genes}  {total_merged_length}\n")
            
            # Write sequences as they come from the iterator
            for gene_name, sequence in gene_iterator:
                outfile.write(f"{gene_name}   {sequence}\n")
                
        logger.info(f"Successfully wrote merged file: {output_file_path}")
    except IOError as e:
        logger.error(f"Error writing output file {output_file_path}: {e}")
        raise
    except Exception as e:
        logger.error(f"An unexpected error occurred during writing to {output_file_path}: {e}")
        raise

def write_merged_phy(final_sequences: Dict[str, str], output_file_path: str, num_genes: int, total_merged_length: int) -> None:
    """
    Writes the merged gene sequences to a new file in PHYLIP format.

    Args:
        final_sequences (dict): Dictionary mapping gene names to merged sequences.
        output_file_path (str): Path for the output PHYLIP file.
        num_genes (int): Total number of unique genes (sequences).
        total_merged_length (int): The total length of the merged sequences.

    Raises:
        IOError: If there's an error writing to the output file.
        Exception: For other unexpected errors during writing.
    """
    logger.info(f"Writing merged data to {output_file_path}...")
    try:
        with open(output_file_path, 'w') as outfile:
            # --- Write Header ---
            # PHYLIP format often expects space padding in the header
            logger.debug(f"  Writing header: {num_genes} sequences, {total_merged_length} total length.")
            outfile.write(f" {num_genes}  {total_merged_length}\n")

            # --- Write Sequences ---
            # Sort gene names for consistent output order (recommended)
            for gene_name in sorted(final_sequences.keys()):
                sequence = final_sequences[gene_name]
                # Consider PHYLIP name length constraints (often 10 chars) if needed:
                # outfile.write(f"{gene_name:<10} {sequence}\n") # Example: Pad name to 10 chars
                outfile.write(f"{gene_name}   {sequence}\n") # Using simple spacing

        logger.info(f"Successfully wrote merged file: {output_file_path}")
    except IOError as e:
        logger.error(f"Error writing output file {output_file_path}: {e}")
        raise # Re-raise the specific error
    except Exception as e:
        logger.error(f"An unexpected error occurred during writing to {output_file_path}: {e}")
        raise # Re-raise after logging


# ==================================================
# Command Line Argument Parsing
# ==================================================
def create_parser() -> argparse.ArgumentParser:
    """Create and configure the argument parser."""
    parser = argparse.ArgumentParser(
        prog='merge_phy.py',
        description='多功能PHYLIP文件处理工具：合并多个PHYLIP文件，支持格式转换和自动对齐',
        epilog='''
示例用法:
  # 合并PHYLIP文件
  python merge_phy.py merged.phy *.phy
  python merge_phy.py -v output.phy gene1.phy gene2.phy gene3.phy  
  python merge_phy.py --gap-char=N --continue-on-error output.phy data/*.phy
  
  # 格式转换
  python merge_phy.py --to-fasta output.fasta input.phy
  python merge_phy.py --from-fasta --to-phy output.phy input.fasta
  
  # 组合操作：合并后转换为FASTA
  python merge_phy.py --to-fasta merged.fasta *.phy

注意:
  - 输入文件必须是有效的PHYLIP格式 (.phy扩展名)
  - 序列会按基因/物种名进行对齐
  - 缺失的序列会用指定的gap字符填充
  - 支持glob模式匹配 (如 *.phy, gene_?.phy)
        ''',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # 位置参数
    parser.add_argument('output_file', 
                       help='输出合并文件的路径')
    parser.add_argument('input_files', nargs='+',
                       help='输入的PHYLIP文件或glob模式 (如 *.phy)')
    
    # 格式转换参数
    format_group = parser.add_mutually_exclusive_group()
    format_group.add_argument('--to-fasta', action='store_true',
                             help='将输出转换为FASTA格式')
    format_group.add_argument('--to-phy', action='store_true',
                             help='强制输出为PHYLIP格式 (默认行为)')
    parser.add_argument('--from-fasta', action='store_true',
                       help='输入文件为FASTA格式 (需要配合--to-phy使用)')
    
    # 合并参数
    parser.add_argument('-g', '--gap-char', default='-', 
                       help='用于填充缺失序列的字符 (默认: "-")')
    parser.add_argument('-c', '--continue-on-error', action='store_true',
                       help='遇到单个文件错误时继续处理其他文件')
    parser.add_argument('--validate-sequences', action='store_true',
                       help='验证序列内容是否为有效的生物序列')
    parser.add_argument('--sequence-type', choices=['auto', 'dna', 'rna', 'protein'], 
                       default='auto', help='序列类型 (默认: auto-detect)')
    
    # 输出控制
    parser.add_argument('-v', '--verbose', action='count', default=0,
                       help='增加输出详细程度 (-v: DEBUG, -vv: 更详细)')
    parser.add_argument('-q', '--quiet', action='store_true',
                       help='静默模式，只显示错误信息')
    parser.add_argument('--statistics', action='store_true',
                       help='显示详细的合并统计信息')
    parser.add_argument('--dry-run', action='store_true',
                       help='预览模式，不实际写入文件')
    
    return parser


# ==================================================
# Main Execution Block
# ==================================================
def main():
    """Main function with improved argument parsing and error handling."""
    # Parse command line arguments
    parser = create_parser()
    args = parser.parse_args()
    
    # Configure logging
    setup_logging(args.verbose, args.quiet)
    
    start_time = time.time()
    
    # --- File Discovery using Glob ---
    # 根据输入格式确定文件后缀
    if args.from_fasta:
        target_suffixes = [".fasta", ".fa", ".fas"]
    else:
        target_suffixes = [".phy"]
    
    found_files_set = set() 

    logger.info(f"搜索匹配的输入文件: {args.input_files}")
    for pattern in args.input_files:
        # Use glob to find files matching the pattern
        matches = glob.glob(pattern)
        if not matches and not os.path.exists(pattern):
             logger.warning(f"输入模式或文件 '{pattern}' 未匹配到任何现有文件")
        for match in matches:
            # Verify that the match is a file and has the correct suffix
            if os.path.isfile(match) and any(match.lower().endswith(suffix) for suffix in target_suffixes):
                found_files_set.add(os.path.abspath(match))
            else:
                logger.debug(f"忽略 '{match}': 不是文件或不以 {target_suffixes} 结尾")

    # Convert to sorted list
    input_phy_files = sorted(list(found_files_set))

    # --- Validation of Found Files ---
    if not input_phy_files:
        print(f"\n错误: 未找到以 {target_suffixes} 结尾的有效输入文件")
        print("提供的参数/模式:", args.input_files)
        sys.exit(1)

    if len(input_phy_files) == 1:
        logger.warning(f"仅找到一个输入文件 ({input_phy_files[0]})，输出将基本上是该文件的重新格式化版本")

    logger.info(f"找到 {len(input_phy_files)} 个有效的输入文件进行处理")
    if logger.isEnabledFor(logging.DEBUG):
        for f_path in input_phy_files:
            logger.debug(f"  - {f_path}")

    # Dry run mode
    if args.dry_run:
        print(f"\n=== 预览模式 ===")
        print(f"将处理 {len(input_phy_files)} 个文件")
        print(f"输出文件: {args.output_file}")
        print(f"Gap字符: '{args.gap_char}'")
        print(f"序列验证: {'是' if args.validate_sequences else '否'}")
        print(f"错误时继续: {'是' if args.continue_on_error else '否'}")
        return

    # --- Main Execution Logic ---
    try:
        # 处理格式转换的特殊情况
        if args.from_fasta and args.to_phy:
            # FASTA到PHYLIP的单文件转换
            if len(input_phy_files) != 1:
                logger.error("FASTA到PHYLIP转换只支持单个输入文件")
                sys.exit(1)
            write_phylip_from_fasta(input_phy_files[0], args.output_file)
            print(f"\n成功将FASTA文件 {input_phy_files[0]} 转换为PHYLIP格式 {args.output_file}")
            return
        
        # 执行合并操作
        final_sequences, num_genes, total_merged_length = merge_phy_files(
            input_phy_files, 
            gap_char=args.gap_char,
            continue_on_error=args.continue_on_error,
            validate_sequences=args.validate_sequences,
            sequence_type=args.sequence_type
        )

        processing_time = time.time() - start_time

        # 检查是否有序列需要写入
        if num_genes > 0:
            # 根据输出格式选择写入方式
            if args.to_fasta:
                write_fasta_file(final_sequences, args.output_file)
                print(f"\n成功合并 {len(input_phy_files)} 个文件并转换为FASTA格式: {args.output_file}")
            else:
                write_merged_phy(final_sequences, args.output_file, num_genes, total_merged_length)
                print(f"\n成功合并 {len(input_phy_files)} 个文件到 {args.output_file}")
            
            # 显示统计信息
            if args.statistics:
                print_statistics(final_sequences, len(input_phy_files), processing_time)
        else:
             print(f"\n合并过程完成，但在输入文件中未找到序列。输出文件 '{args.output_file}' 未创建。")

    except SystemExit:
        print("\n由于文件处理错误，脚本执行被中止。")
        sys.exit(1)
    except Exception as e:
        print(f"\n合并或写入过程中发生意外错误: {e}")
        logger.exception("完整错误信息:")
        sys.exit(1)


if __name__ == "__main__":
    main()

