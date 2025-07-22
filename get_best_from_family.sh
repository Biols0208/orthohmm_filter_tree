#!/bin/bash

# 基因家族序列处理流水线 - 优化版本
# 完整的序列校正、质量分析和清理工作流
# Author: Biols9527
# Date: 2025-06-24
# Version: 1.3 - 增强脚本路径配置和其他优化

set -e  # 遇到错误立即退出

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# 打印带颜色的信息
print_info() {
    echo -e "${BLUE}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

print_step() {
    echo -e "${PURPLE}[STEP]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

# 显示使用帮助
show_help() {
    cat << EOF
基因家族序列处理流水线 v1.3
Author: Biols9527

用法: $0 [选项] <输入文件>

必需参数:
  输入文件                    FASTA格式的基因家族序列文件

可选参数:
  -h, --help                 显示此帮助信息
  -o, --output PREFIX        输出前缀 (默认: 基于输入文件名_YYYYMMDD_HHMMSS)
  -d, --output-dir DIR       输出目录 (默认: 当前目录)
  -t, --threads NUM          线程数 (默认: 4)
  -e, --evalue VALUE         BLAST E-value阈值 (默认: 1e-5)
  -s, --sequence-type TYPE   序列类型 [nucleotide|protein] (默认: nucleotide)
  -m, --method METHOD        校正方法 [consensus|longest_reference|majority_vote] (默认: consensus)
  --min-overlap NUM          最小重叠长度 (默认: 20)
  
  脚本路径配置:
  --step1-script PATH        step1校正脚本路径 (默认: chek_correction.py)
  --step2-script PATH        step2质量分析脚本路径 (默认: quality_analyzer.py)
  --python-cmd CMD           Python命令 (默认: python3)
  --scripts-dir DIR          脚本目录 (会自动查找脚本)
  
  处理控制:
  --skip-correction          跳过序列校正步骤
  --skip-quality             跳过质量分析步骤
  --keep-intermediate        保留中间文件
  --compress                 压缩输出文件
  --create-subdirs           创建结构化子目录
  --no-timestamp             不在输出文件名中包含时间戳
  --force                    强制覆盖已存在的输出
  
  运行模式:
  -v, --verbose              详细输出模式
  --dry-run                  仅显示将要执行的命令，不实际运行
  --config FILE              从配置文件读取参数
  --save-config FILE         保存当前配置到文件
  
  质量控制:
  --max-memory GB            最大内存使用限制 (默认: 8)
  --temp-dir DIR             临时文件目录 (默认: /tmp)
  --log-level LEVEL          日志级别 [DEBUG|INFO|WARNING|ERROR] (默认: INFO)

示例:
  基本使用:
    $0 OG04298.cds
  
  指定脚本路径:
    $0 --step1-script /path/to/step1.py --step2-script /path/to/step2.py OG04298.cds
    $0 --scripts-dir /opt/gene_tools/ OG04298.cds
  
  完整配置:
    $0 -o my_analysis -d /path/to/output -t 8 --create-subdirs \\
       --step1-script ./scripts/step1_chek_correction.py OG04298.cds
  
  使用配置文件:
    $0 --config analysis.conf OG04298.cds

输出文件结构:
  标准模式:
    \${PREFIX}_best_sequences.fasta          最终最佳序列
    \${PREFIX}_removed_sequences.fasta       被移除的序列
    \${PREFIX}_pipeline_report.txt           完整流水线报告
    \${BASE_NAME}.best                       兼容性输出

  子目录模式 (--create-subdirs):
    results/\${PREFIX}_best_sequences.fasta
    results/\${PREFIX}_removed_sequences.fasta
    intermediate/\${PREFIX}_corrected.fasta
    reports/\${PREFIX}_pipeline_report.txt
    logs/\${PREFIX}_processing.log

EOF
}

# 默认参数
OUTPUT_PREFIX=""
OUTPUT_DIR="."
THREADS=4
EVALUE="1e-5"
SEQUENCE_TYPE="nucleotide"
CORRECTION_METHOD="consensus"
MIN_OVERLAP=20

# 脚本路径配置
STEP1_SCRIPT="step1_chek_correction.py"
STEP2_SCRIPT="step2_quality_analyzer.py"
PYTHON_CMD="/dellfsqd2/ST_OCEAN/USER/lishuo11/01_soft/mambaforge/bin/python"
SCRIPTS_DIR="/dellfsqd2/ST_OCEAN/USER/lishuo1/08_shell"

# 处理控制
SKIP_CORRECTION=false
SKIP_QUALITY=false
KEEP_INTERMEDIATE=false
COMPRESS=false
CREATE_SUBDIRS=false
NO_TIMESTAMP=false
FORCE=false

# 运行模式
VERBOSE=false
DRY_RUN=false
CONFIG_FILE=""
SAVE_CONFIG=""

# 质量控制
MAX_MEMORY=8
TEMP_DIR="/tmp"
LOG_LEVEL="INFO"

# 初始化统计变量
INPUT_TOTAL_SEQS=0
INPUT_TOTAL_SPECIES=0
CORRECTION_SEQS=0
FINAL_SEQS=0
REMOVED_SEQS=0

# 全局路径变量
BASE_NAME=""
TIMESTAMP=""
ACTUAL_OUTPUT_PREFIX=""
RESULTS_DIR=""
INTERMEDIATE_DIR=""
REPORTS_DIR=""
LOGS_DIR=""
LOG_FILE=""

# 日志函数
write_log() {
    local level="$1"
    local message="$2"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    
    if [ -n "$LOG_FILE" ]; then
        echo "[$timestamp] [$level] $message" >> "$LOG_FILE"
    fi
    
    # 根据日志级别决定是否输出到控制台
    case "$LOG_LEVEL" in
        "DEBUG")
            ;;
        "INFO")
            if [ "$level" != "DEBUG" ]; then
                case "$level" in
                    "INFO") print_info "$message" ;;
                    "WARNING") print_warning "$message" ;;
                    "ERROR") print_error "$message" ;;
                esac
            fi
            ;;
        "WARNING")
            if [ "$level" = "WARNING" ] || [ "$level" = "ERROR" ]; then
                case "$level" in
                    "WARNING") print_warning "$message" ;;
                    "ERROR") print_error "$message" ;;
                esac
            fi
            ;;
        "ERROR")
            if [ "$level" = "ERROR" ]; then
                print_error "$message"
            fi
            ;;
    esac
}

# 加载配置文件
load_config() {
    local config_file="$1"
    if [ ! -f "$config_file" ]; then
        print_error "配置文件不存在: $config_file"
        exit 1
    fi
    
    print_info "加载配置文件: $config_file"
    # 使用source加载配置，但只允许特定变量
    while IFS='=' read -r key value; do
        # 跳过注释和空行
        if [[ $key =~ ^[[:space:]]*# ]] || [[ -z "$key" ]]; then
            continue
        fi
        
        # 去除前后空白
        key=$(echo "$key" | xargs)
        value=$(echo "$value" | xargs)
        
        # 只允许设置特定的变量
        case "$key" in
            OUTPUT_PREFIX|OUTPUT_DIR|THREADS|EVALUE|SEQUENCE_TYPE|CORRECTION_METHOD|MIN_OVERLAP|\
            STEP1_SCRIPT|STEP2_SCRIPT|PYTHON_CMD|SCRIPTS_DIR|MAX_MEMORY|TEMP_DIR|LOG_LEVEL)
                eval "$key=\"$value\""
                write_log "DEBUG" "从配置文件设置 $key=$value"
                ;;
            SKIP_CORRECTION|SKIP_QUALITY|KEEP_INTERMEDIATE|COMPRESS|CREATE_SUBDIRS|NO_TIMESTAMP|FORCE|VERBOSE)
                eval "$key=$value"
                write_log "DEBUG" "从配置文件设置 $key=$value"
                ;;
        esac
    done < "$config_file"
}

# 保存配置文件
save_config() {
    local config_file="$1"
    
    cat > "$config_file" << EOF
# 基因家族序列处理流水线配置文件
# 生成时间: $(date '+%Y-%m-%d %H:%M:%S')

# 基本参数
OUTPUT_PREFIX=$OUTPUT_PREFIX
OUTPUT_DIR=$OUTPUT_DIR
THREADS=$THREADS
EVALUE=$EVALUE
SEQUENCE_TYPE=$SEQUENCE_TYPE
CORRECTION_METHOD=$CORRECTION_METHOD
MIN_OVERLAP=$MIN_OVERLAP

# 脚本路径
STEP1_SCRIPT=$STEP1_SCRIPT
STEP2_SCRIPT=$STEP2_SCRIPT
PYTHON_CMD=$PYTHON_CMD
SCRIPTS_DIR=$SCRIPTS_DIR

# 处理控制
SKIP_CORRECTION=$SKIP_CORRECTION
SKIP_QUALITY=$SKIP_QUALITY
KEEP_INTERMEDIATE=$KEEP_INTERMEDIATE
COMPRESS=$COMPRESS
CREATE_SUBDIRS=$CREATE_SUBDIRS
NO_TIMESTAMP=$NO_TIMESTAMP
FORCE=$FORCE

# 质量控制
MAX_MEMORY=$MAX_MEMORY
TEMP_DIR=$TEMP_DIR
LOG_LEVEL=$LOG_LEVEL
VERBOSE=$VERBOSE
EOF
    
    print_success "配置已保存到: $config_file"
}

# 安全的数值比较函数
compare_numbers() {
    local num1="$1"
    local op="$2"
    local num2="$3"

    # 确保数值有效
    if [ -z "$num1" ]; then num1=0; fi
    if [ -z "$num2" ]; then num2=0; fi

    case "$op" in
        "gt"|">")
            [ "$num1" -gt "$num2" ]
            ;;
        "lt"|"<")
            [ "$num1" -lt "$num2" ]
            ;;
        "eq"|"=")
            [ "$num1" -eq "$num2" ]
            ;;
        "ge"|">=")
            [ "$num1" -ge "$num2" ]
            ;;
        "le"|"<=")
            [ "$num1" -le "$num2" ]
            ;;
        *)
            return 1
            ;;
    esac
}

# 解析命令行参数
parse_arguments() {
    while [ $# -gt 0 ]; do
        case $1 in
            -h|--help)
                show_help
                exit 0
                ;;
            -o|--output)
                OUTPUT_PREFIX="$2"
                shift 2
                ;;
            -d|--output-dir)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            -t|--threads)
                THREADS="$2"
                shift 2
                ;;
            -e|--evalue)
                EVALUE="$2"
                shift 2
                ;;
            -s|--sequence-type)
                SEQUENCE_TYPE="$2"
                shift 2
                ;;
            -m|--method)
                CORRECTION_METHOD="$2"
                shift 2
                ;;
            --min-overlap)
                MIN_OVERLAP="$2"
                shift 2
                ;;
            --step1-script)
                STEP1_SCRIPT="$2"
                shift 2
                ;;
            --step2-script)
                STEP2_SCRIPT="$2"
                shift 2
                ;;
            --python-cmd)
                PYTHON_CMD="$2"
                shift 2
                ;;
            --scripts-dir)
                SCRIPTS_DIR="$2"
                shift 2
                ;;
            --skip-correction)
                SKIP_CORRECTION=true
                shift
                ;;
            --skip-quality)
                SKIP_QUALITY=true
                shift
                ;;
            --keep-intermediate)
                KEEP_INTERMEDIATE=true
                shift
                ;;
            --compress)
                COMPRESS=true
                shift
                ;;
            --create-subdirs)
                CREATE_SUBDIRS=true
                shift
                ;;
            --no-timestamp)
                NO_TIMESTAMP=true
                shift
                ;;
            --force)
                FORCE=true
                shift
                ;;
            -v|--verbose)
                VERBOSE=true
                shift
                ;;
            --dry-run)
                DRY_RUN=true
                shift
                ;;
            --config)
                CONFIG_FILE="$2"
                shift 2
                ;;
            --save-config)
                SAVE_CONFIG="$2"
                shift 2
                ;;
            --max-memory)
                MAX_MEMORY="$2"
                shift 2
                ;;
            --temp-dir)
                TEMP_DIR="$2"
                shift 2
                ;;
            --log-level)
                LOG_LEVEL="$2"
                shift 2
                ;;
            -*)
                print_error "未知选项: $1"
                show_help
                exit 1
                ;;
            *)
                if [ -z "$INPUT_FILE" ]; then
                    INPUT_FILE="$1"
                else
                    print_error "多余的参数: $1"
                    exit 1
                fi
                shift
                ;;
        esac
    done
}

# 解析脚本路径
resolve_script_paths() {
    # 如果指定了脚本目录，则在该目录下查找脚本
    if [ -n "$SCRIPTS_DIR" ]; then
        if [ ! -d "$SCRIPTS_DIR" ]; then
            print_error "脚本目录不存在: $SCRIPTS_DIR"
            exit 1
        fi
        
        # 更新脚本路径
        if [ "$STEP1_SCRIPT" = "step1_chek_correction.py" ]; then
            STEP1_SCRIPT="$SCRIPTS_DIR/step1_chek_correction.py"
        fi
        if [ "$STEP2_SCRIPT" = "step2_quality_analyzer.py" ]; then
            STEP2_SCRIPT="$SCRIPTS_DIR/step2_quality_analyzer.py"
        fi
    fi
    
    # 转换为绝对路径
    if [ ! "${STEP1_SCRIPT:0:1}" = "/" ]; then
        STEP1_SCRIPT="$(pwd)/$STEP1_SCRIPT"
    fi
    if [ ! "${STEP2_SCRIPT:0:1}" = "/" ]; then
        STEP2_SCRIPT="$(pwd)/$STEP2_SCRIPT"
    fi
    
    write_log "DEBUG" "Step1 脚本路径: $STEP1_SCRIPT"
    write_log "DEBUG" "Step2 脚本路径: $STEP2_SCRIPT"
}

# 初始化输出路径和文件名
initialize_output_paths() {
    # 获取基础文件名（去除路径和扩展名）
    BASE_NAME=$(basename "$INPUT_FILE")
    BASE_NAME="${BASE_NAME%.*}"  # 移除扩展名
    
    # 生成时间戳
    TIMESTAMP=$(date '+%Y%m%d_%H%M%S')
    
    # 确定输出前缀
    if [ -z "$OUTPUT_PREFIX" ]; then
        if [ "$NO_TIMESTAMP" = true ]; then
            ACTUAL_OUTPUT_PREFIX="${BASE_NAME}"
        else
            ACTUAL_OUTPUT_PREFIX="${BASE_NAME}_${TIMESTAMP}"
        fi
    else
        if [ "$NO_TIMESTAMP" = true ]; then
            ACTUAL_OUTPUT_PREFIX="$OUTPUT_PREFIX"
        else
            ACTUAL_OUTPUT_PREFIX="${OUTPUT_PREFIX}_${TIMESTAMP}"
        fi
    fi
    
    # 创建输出目录
    if [ ! -d "$OUTPUT_DIR" ]; then
        mkdir -p "$OUTPUT_DIR"
    fi
    
    # 检查输出目录是否已存在相同前缀的文件
    if [ "$FORCE" = false ]; then
        local existing_files=$(find "$OUTPUT_DIR" -name "${ACTUAL_OUTPUT_PREFIX}*" -type f 2>/dev/null | head -1)
        if [ -n "$existing_files" ]; then
            print_warning "发现已存在的输出文件: $existing_files"
            print_warning "使用 --force 选项强制覆盖，或更改输出前缀"
            read -p "是否继续？(y/N): " -n 1 -r
            echo
            if [[ ! $REPLY =~ ^[Yy]$ ]]; then
                print_info "操作已取消"
                exit 0
            fi
        fi
    fi
    
    # 设置目录结构
    if [ "$CREATE_SUBDIRS" = true ]; then
        RESULTS_DIR="$OUTPUT_DIR/results"
        INTERMEDIATE_DIR="$OUTPUT_DIR/intermediate"
        REPORTS_DIR="$OUTPUT_DIR/reports"
        LOGS_DIR="$OUTPUT_DIR/logs"
        
        mkdir -p "$RESULTS_DIR" "$INTERMEDIATE_DIR" "$REPORTS_DIR" "$LOGS_DIR"
    else
        RESULTS_DIR="$OUTPUT_DIR"
        INTERMEDIATE_DIR="$OUTPUT_DIR"
        REPORTS_DIR="$OUTPUT_DIR"
        LOGS_DIR="$OUTPUT_DIR"
    fi
    
    # 设置日志文件
    LOG_FILE="${LOGS_DIR}/${ACTUAL_OUTPUT_PREFIX}_processing.log"
    
    write_log "INFO" "输出配置初始化完成"
    write_log "INFO" "基础名称: $BASE_NAME"
    write_log "INFO" "输出前缀: $ACTUAL_OUTPUT_PREFIX"
    write_log "INFO" "输出目录: $OUTPUT_DIR"
    if [ "$CREATE_SUBDIRS" = true ]; then
        write_log "INFO" "使用子目录结构"
    fi
}

# 验证输入参数
validate_parameters() {
    write_log "INFO" "验证输入参数..."
    
    # 检查输入文件
    if [ -z "$INPUT_FILE" ]; then
        print_error "缺少输入文件参数"
        show_help
        exit 1
    fi

    if [ ! -f "$INPUT_FILE" ]; then
        print_error "输入文件不存在: $INPUT_FILE"
        exit 1
    fi

    # 验证序列类型
    if [ "$SEQUENCE_TYPE" != "nucleotide" ] && [ "$SEQUENCE_TYPE" != "protein" ]; then
        print_error "序列类型必须是 'nucleotide' 或 'protein'"
        exit 1
    fi

    # 验证校正方法
    case "$CORRECTION_METHOD" in
        consensus|longest_reference|majority_vote)
            ;;
        *)
            print_error "校正方法必须是 'consensus', 'longest_reference' 或 'majority_vote'"
            exit 1
            ;;
    esac

    # 验证数值参数
    if ! echo "$THREADS" | grep -q '^[0-9][0-9]*$' || compare_numbers "$THREADS" "lt" "1"; then
        print_error "线程数必须是正整数"
        exit 1
    fi

    if ! echo "$MIN_OVERLAP" | grep -q '^[0-9][0-9]*$' || compare_numbers "$MIN_OVERLAP" "lt" "1"; then
        print_error "最小重叠长度必须是正整数"
        exit 1
    fi
    
    if ! echo "$MAX_MEMORY" | grep -q '^[0-9][0-9]*$' || compare_numbers "$MAX_MEMORY" "lt" "1"; then
        print_error "最大内存必须是正整数(GB)"
        exit 1
    fi
    
    # 验证日志级别
    case "$LOG_LEVEL" in
        DEBUG|INFO|WARNING|ERROR)
            ;;
        *)
            print_error "日志级别必须是 DEBUG, INFO, WARNING 或 ERROR"
            exit 1
            ;;
    esac
    
    # 验证输出目录是否可写
    if [ ! -w "$(dirname "$OUTPUT_DIR")" ]; then
        print_error "无法在指定位置创建输出目录: $OUTPUT_DIR"
        exit 1
    fi
    
    # 验证临时目录
    if [ ! -d "$TEMP_DIR" ] || [ ! -w "$TEMP_DIR" ]; then
        print_warning "临时目录不可用: $TEMP_DIR，使用当前目录"
        TEMP_DIR="."
    fi
    
    write_log "INFO" "参数验证通过"
}

# 检查依赖软件
check_dependencies() {
    write_log "INFO" "检查依赖软件..."

    local missing_deps=""

    # 检查Python命令
    if ! command -v "$PYTHON_CMD" >/dev/null 2>&1; then
        missing_deps="$missing_deps $PYTHON_CMD"
    fi

    # 检查Python脚本
    if [ "$SKIP_CORRECTION" = false ]; then
        if [ ! -f "$STEP1_SCRIPT" ]; then
            missing_deps="$missing_deps $STEP1_SCRIPT"
        fi
    fi

    if [ "$SKIP_QUALITY" = false ]; then
        if [ ! -f "$STEP2_SCRIPT" ]; then
            missing_deps="$missing_deps $STEP2_SCRIPT"
        fi
    fi

    # 检查BLAST工具（如果需要质量分析）
    if [ "$SKIP_QUALITY" = false ]; then
        for cmd in makeblastdb blastn blastp; do
            if ! command -v "$cmd" >/dev/null 2>&1; then
                write_log "WARNING" "BLAST工具 $cmd 未找到，质量分析可能受限"
            fi
        done
    fi

    # 检查MAFFT（如果需要序列校正）
    if [ "$SKIP_CORRECTION" = false ]; then
        if ! command -v mafft >/dev/null 2>&1; then
            write_log "WARNING" "MAFFT未找到，序列校正可能受限"
        fi
    fi

    # 检查系统资源
    local available_memory=$(free -g | awk '/^Mem:/{print $7}' 2>/dev/null || echo "unknown")
    if [ "$available_memory" != "unknown" ] && compare_numbers "$available_memory" "lt" "$MAX_MEMORY"; then
        write_log "WARNING" "可用内存($available_memory GB)小于设置的最大内存($MAX_MEMORY GB)"
    fi

    if [ -n "$missing_deps" ]; then
        print_error "缺少以下关键依赖:"
        for dep in $missing_deps; do
            echo "  - $dep"
        done
        exit 1
    fi

    write_log "INFO" "依赖检查通过"
}

# 执行命令（支持干运行模式）
execute_command() {
    local cmd="$1"
    local description="$2"

    write_log "DEBUG" "准备执行: $description"
    write_log "DEBUG" "命令: $cmd"

    if [ "$VERBOSE" = true ]; then
        print_info "执行: $description"
        echo "命令: $cmd"
    fi

    if [ "$DRY_RUN" = true ]; then
        echo "[DRY-RUN] $cmd"
        return 0
    fi

    local start_time=$(date +%s)
    if eval "$cmd"; then
        local end_time=$(date +%s)
        local duration=$((end_time - start_time))
        write_log "INFO" "$description 完成 (耗时: ${duration}s)"
        return 0
    else
        write_log "ERROR" "$description 失败"
        return 1
    fi
}

# 分析输入文件
analyze_input() {
    write_log "INFO" "分析输入文件..."

    INPUT_TOTAL_SEQS=$(grep -c "^>" "$INPUT_FILE" 2>/dev/null || echo "0")
    INPUT_TOTAL_SPECIES=$(grep "^>" "$INPUT_FILE" | cut -d' ' -f2- | sort -u | wc -l 2>/dev/null || echo "0")
    
    local file_size=$(du -h "$INPUT_FILE" | cut -f1)

    write_log "INFO" "输入文件: $INPUT_FILE ($file_size)"
    write_log "INFO" "总序列数: $INPUT_TOTAL_SEQS"
    write_log "INFO" "物种数: $INPUT_TOTAL_SPECIES"

    if compare_numbers "$INPUT_TOTAL_SEQS" "eq" "0"; then
        print_error "输入文件中没有找到序列"
        exit 1
    fi

    # 检查序列头部格式
    local first_header=$(head -1 "$INPUT_FILE")
    if ! echo "$first_header" | grep -q '^>.*[[:space:]].*'; then
        write_log "WARNING" "序列头部格式可能不标准，建议格式: >gene_name species_name"
    fi
    
    # 估算内存需求
    local estimated_memory=$((INPUT_TOTAL_SEQS / 1000 + 1))
    if compare_numbers "$estimated_memory" "gt" "$MAX_MEMORY"; then
        write_log "WARNING" "估算内存需求(${estimated_memory}GB)可能超过设置限制(${MAX_MEMORY}GB)"
    fi
}

# 步骤1: 序列校正
run_sequence_correction() {
    if [ "$SKIP_CORRECTION" = true ]; then
        write_log "INFO" "跳过序列校正步骤"
        CORRECTED_FILE="$INPUT_FILE"
        return 0
    fi

    write_log "INFO" "开始序列校正和完善..."

    local correction_prefix="${INTERMEDIATE_DIR}/${ACTUAL_OUTPUT_PREFIX}_correction"
    local cmd="\"$PYTHON_CMD\" \"$STEP1_SCRIPT\" \"$INPUT_FILE\" \
        --output_prefix \"$correction_prefix\" \
        --correction_method \"$CORRECTION_METHOD\" \
        --min_overlap $MIN_OVERLAP"

    if ! execute_command "$cmd" "序列校正"; then
        print_error "序列校正失败"
        exit 1
    fi

    CORRECTED_FILE="${correction_prefix}_corrected_sequences.fasta"

    if [ ! -f "$CORRECTED_FILE" ]; then
        print_error "校正后的序列文件未生成: $CORRECTED_FILE"
        exit 1
    fi

    CORRECTION_SEQS=$(grep -c "^>" "$CORRECTED_FILE" 2>/dev/null || echo "0")
    write_log "INFO" "序列校正完成，输出 $CORRECTION_SEQS 个序列"
}

# 步骤2: 质量分析
run_quality_analysis() {
    if [ "$SKIP_QUALITY" = true ]; then
        write_log "INFO" "跳过质量分析步骤"
        FINAL_FILE="$CORRECTED_FILE"
        return 0
    fi

    write_log "INFO" "开始全局BLAST质量分析..."

    local quality_prefix="${INTERMEDIATE_DIR}/${ACTUAL_OUTPUT_PREFIX}_quality"
    local cmd="\"$PYTHON_CMD\" \"$STEP2_SCRIPT\" \"$CORRECTED_FILE\" \
        --sequence_type \"$SEQUENCE_TYPE\" \
        --output_prefix \"$quality_prefix\" \
        --evalue $EVALUE \
        --num_threads $THREADS"

    if ! execute_command "$cmd" "质量分析"; then
        print_error "质量分析失败"
        exit 1
    fi

    FINAL_FILE="${quality_prefix}_best_sequences.fasta"
    REMOVED_FILE="${quality_prefix}_removed_sequences.fasta"

    if [ ! -f "$FINAL_FILE" ]; then
        print_error "最终序列文件未生成: $FINAL_FILE"
        exit 1
    fi

    FINAL_SEQS=$(grep -c "^>" "$FINAL_FILE" 2>/dev/null || echo "0")
    REMOVED_SEQS=0
    if [ -f "$REMOVED_FILE" ]; then
        REMOVED_SEQS=$(grep -c "^>" "$REMOVED_FILE" 2>/dev/null || echo "0")
    fi

    write_log "INFO" "质量分析完成，选择 $FINAL_SEQS 个最佳序列，移除 $REMOVED_SEQS 个序列"
}

# 步骤3: 整理输出文件
organize_output() {
    write_log "INFO" "整理输出文件..."

    # 定义最终输出文件路径
    local final_best_sequences="${RESULTS_DIR}/${ACTUAL_OUTPUT_PREFIX}_best_sequences.fasta"
    local final_removed_sequences="${RESULTS_DIR}/${ACTUAL_OUTPUT_PREFIX}_removed_sequences.fasta"
    local compatibility_output="${RESULTS_DIR}/${BASE_NAME}.best"

    # 复制最终序列文件
    if [ -f "$FINAL_FILE" ]; then
        execute_command "cp \"$FINAL_FILE\" \"$final_best_sequences\"" "复制最终序列文件"
        
        # 创建兼容性输出文件
        execute_command "cp \"$FINAL_FILE\" \"$compatibility_output\"" "创建兼容输出文件"
    fi

    # 复制移除的序列文件
    if [ -f "$REMOVED_FILE" ]; then
        execute_command "cp \"$REMOVED_FILE\" \"$final_removed_sequences\"" "复制移除的序列文件"
    fi

    # 更新全局变量以指向最终位置
    FINAL_BEST_SEQUENCES="$final_best_sequences"
    FINAL_REMOVED_SEQUENCES="$final_removed_sequences"
    COMPATIBILITY_OUTPUT="$compatibility_output"
}

# 生成流水线报告
generate_pipeline_report() {
    write_log "INFO" "生成流水线报告..."

    local report_file="${REPORTS_DIR}/${ACTUAL_OUTPUT_PREFIX}_pipeline_report.txt"
    local manifest_file="${REPORTS_DIR}/${ACTUAL_OUTPUT_PREFIX}_file_manifest.txt"
    local config_backup="${REPORTS_DIR}/${ACTUAL_OUTPUT_PREFIX}_config_backup.conf"

    # 备份当前配置
    save_config "$config_backup"

    # 生成详细报告
    cat > "$report_file" << EOF
==========================================
基因家族序列处理流水线报告
==========================================
运行时间: $(date '+%Y-%m-%d %H:%M:%S %Z')
用户: $(whoami)
主机: $(hostname)
工作目录: $(pwd)
流水线版本: 1.3

输入信息:
----------
输入文件: $INPUT_FILE
文件大小: $(du -h "$INPUT_FILE" | cut -f1)
输出目录: $OUTPUT_DIR
输出前缀: $ACTUAL_OUTPUT_PREFIX

运行参数:
----------
序列类型: $SEQUENCE_TYPE
线程数: $THREADS
BLAST E-value: $EVALUE
校正方法: $CORRECTION_METHOD
最小重叠: $MIN_OVERLAP
最大内存限制: ${MAX_MEMORY}GB
Python命令: $PYTHON_CMD

脚本路径:
----------
Step1 校正脚本: $STEP1_SCRIPT
Step2 质量分析脚本: $STEP2_SCRIPT

处理步骤:
----------
EOF

    if [ "$SKIP_CORRECTION" = true ]; then
        echo "1. 序列校正: 跳过" >> "$report_file"
    else
        echo "1. 序列校正: 执行 ($CORRECTION_METHOD 方法)" >> "$report_file"
    fi

    if [ "$SKIP_QUALITY" = true ]; then
        echo "2. 质量分析: 跳过" >> "$report_file"
    else
        echo "2. 质量分析: 执行 (E-value: $EVALUE, 线程: $THREADS)" >> "$report_file"
    fi

    cat >> "$report_file" << EOF

统计结果:
----------
输入序列总数: $INPUT_TOTAL_SEQS
输入物种数: $INPUT_TOTAL_SPECIES
EOF

    if [ "$SKIP_CORRECTION" = false ]; then
        echo "校正后序列数: $CORRECTION_SEQS" >> "$report_file"
    fi

    local efficiency_rate="N/A"
    if compare_numbers "$INPUT_TOTAL_SEQS" "gt" "0"; then
        efficiency_rate=$(echo "scale=2; $FINAL_SEQS * 100 / $INPUT_TOTAL_SEQS" | bc 2>/dev/null || echo "N/A")
        efficiency_rate="${efficiency_rate}%"
    fi

    cat >> "$report_file" << EOF
最终保留序列数: $FINAL_SEQS
移除序列数: $REMOVED_SEQS
处理效率: $efficiency_rate

系统信息:
----------
可用内存: $(free -h | awk '/^Mem:/{print $7}' 2>/dev/null || echo "unknown")
CPU核心数: $(nproc 2>/dev/null || echo "unknown")
磁盘可用空间: $(df -h "$OUTPUT_DIR" | awk 'NR==2{print $4}' 2>/dev/null || echo "unknown")

输出文件位置:
-------------
EOF

    # 生成文件清单
    cat > "$manifest_file" << EOF
文件清单 - $(date '+%Y-%m-%d %H:%M:%S')
========================================

主要结果文件:
EOF

    # 列出实际生成的文件
    if [ -f "$FINAL_BEST_SEQUENCES" ]; then
        local size=$(du -h "$FINAL_BEST_SEQUENCES" | cut -f1)
        echo "✓ $FINAL_BEST_SEQUENCES (${size}) - 最终最佳序列" >> "$report_file"
        echo "- $(basename "$FINAL_BEST_SEQUENCES") (${size}) - 最终最佳序列" >> "$manifest_file"
    fi

    if [ -f "$FINAL_REMOVED_SEQUENCES" ]; then
        local size=$(du -h "$FINAL_REMOVED_SEQUENCES" | cut -f1)
        echo "✓ $FINAL_REMOVED_SEQUENCES (${size}) - 移除的序列" >> "$report_file"
        echo "- $(basename "$FINAL_REMOVED_SEQUENCES") (${size}) - 移除的序列" >> "$manifest_file"
    fi

    if [ -f "$COMPATIBILITY_OUTPUT" ]; then
        local size=$(du -h "$COMPATIBILITY_OUTPUT" | cut -f1)
        echo "✓ $COMPATIBILITY_OUTPUT (${size}) - 兼容性输出" >> "$report_file"
        echo "- $(basename "$COMPATIBILITY_OUTPUT") (${size}) - 兼容性输出" >> "$manifest_file"
    fi

    echo "✓ $report_file - 详细处理报告" >> "$report_file"
    echo "✓ $manifest_file - 文件清单" >> "$report_file"
    echo "✓ $config_backup - 配置备份" >> "$report_file"

    cat >> "$manifest_file" << EOF

报告文件:
- $(basename "$report_file") - 详细处理报告
- $(basename "$manifest_file") - 此文件清单
- $(basename "$config_backup") - 运行时配置备份

日志文件:
- $(basename "$LOG_FILE") - 处理日志

中间文件 (如果保留):
EOF

    if [ "$KEEP_INTERMEDIATE" = true ]; then
        for file in "${INTERMEDIATE_DIR}"/*; do
            if [ -f "$file" ]; then
                local size=$(du -h "$file" | cut -f1)
                echo "- $(basename "$file") (${size})" >> "$manifest_file"
            fi
        done
    else
        echo "- (已清理)" >> "$manifest_file"
    fi

    write_log "INFO" "报告生成完成: $report_file, $manifest_file"
}

# 压缩输出文件
compress_outputs() {
    if [ "$COMPRESS" = false ]; then
        return 0
    fi

    write_log "INFO" "压缩输出文件..."

    local archive_name="${RESULTS_DIR}/${ACTUAL_OUTPUT_PREFIX}_analysis.tar.gz"
    local files_to_compress=""

    # 收集要压缩的文件
    if [ -f "$FINAL_BEST_SEQUENCES" ]; then
        files_to_compress="$files_to_compress $FINAL_BEST_SEQUENCES"
    fi
    if [ -f "$FINAL_REMOVED_SEQUENCES" ]; then
        files_to_compress="$files_to_compress $FINAL_REMOVED_SEQUENCES"
    fi
    if [ -f "$COMPATIBILITY_OUTPUT" ]; then
        files_to_compress="$files_to_compress $COMPATIBILITY_OUTPUT"
    fi

    # 添加报告文件
    for file in "${REPORTS_DIR}"/*; do
        if [ -f "$file" ]; then
            files_to_compress="$files_to_compress $file"
        fi
    done

    if [ -n "$files_to_compress" ]; then
        execute_command "tar zcf \"$archive_name\" $files_to_compress --remove-files" "压缩输出文件"
        write_log "INFO" "压缩文件创建: $archive_name"
    fi
}

# 清理中间文件
cleanup_intermediate() {
    if [ "$KEEP_INTERMEDIATE" = true ]; then
        write_log "INFO" "保留中间文件"
        return 0
    fi

    write_log "INFO" "清理中间文件..."

    local cleaned_count=0
    
    # 删除中间目录中的文件（但保留最终输出）
    if [ "$CREATE_SUBDIRS" = true ]; then
        for file in "${INTERMEDIATE_DIR}"/*; do
            if [ -f "$file" ]; then
                write_log "DEBUG" "删除中间文件: $(basename "$file")"
                rm -f "$file"
                cleaned_count=$((cleaned_count + 1))
            fi
        done
    else
        # 在单目录模式下，删除特定的中间文件
        for pattern in "*_correction_*" "*_quality_*"; do
            for file in ${OUTPUT_DIR}/$pattern; do
                if [ -f "$file" ] && [[ ! "$file" =~ _best_sequences\.fasta$ ]] && [[ ! "$file" =~ _removed_sequences\.fasta$ ]]; then
                    write_log "DEBUG" "删除中间文件: $(basename "$file")"
                    rm -f "$file"
                    cleaned_count=$((cleaned_count + 1))
                fi
            done
        done
    fi

    write_log "INFO" "清理了 $cleaned_count 个中间文件"
}

# 创建便捷访问链接
create_convenience_links() {
    if [ "$CREATE_SUBDIRS" = false ]; then
        return 0
    fi

    write_log "INFO" "创建便捷访问链接..."

    # 在输出目录根部创建到最重要文件的链接
    if [ -f "$FINAL_BEST_SEQUENCES" ]; then
        execute_command "ln -sf \"results/$(basename "$FINAL_BEST_SEQUENCES")\" \"${OUTPUT_DIR}/latest_best_sequences.fasta\"" "创建最佳序列链接"
    fi

    if [ -f "$COMPATIBILITY_OUTPUT" ]; then
        execute_command "ln -sf \"results/$(basename "$COMPATIBILITY_OUTPUT")\" \"${OUTPUT_DIR}/$(basename "$COMPATIBILITY_OUTPUT")\"" "创建兼容性输出链接"
    fi

    local latest_report="${REPORTS_DIR}/${ACTUAL_OUTPUT_PREFIX}_pipeline_report.txt"
    if [ -f "$latest_report" ]; then
        execute_command "ln -sf \"reports/$(basename "$latest_report")\" \"${OUTPUT_DIR}/latest_report.txt\"" "创建最新报告链接"
    fi

    write_log "INFO" "便捷链接创建完成"
}

# 显示最终结果
show_final_results() {
    echo
    print_success "=========================================="
    print_success "基因家族序列处理流水线完成！"
    print_success "=========================================="
    echo
    print_info "处理总结:"
    print_info "  输入序列: $INPUT_TOTAL_SEQS 个 (来自 $INPUT_TOTAL_SPECIES 个物种)"
    if [ "$SKIP_CORRECTION" = false ]; then
        print_info "  校正后序列: $CORRECTION_SEQS 个"
    fi
    print_info "  最终序列: $FINAL_SEQS 个 (每个物种的最佳代表)"
    if compare_numbers "$REMOVED_SEQS" "gt" "0"; then
        print_info "  移除序列: $REMOVED_SEQS 个"
    fi
    echo
    print_info "输出位置: $OUTPUT_DIR"
    print_info "输出前缀: $ACTUAL_OUTPUT_PREFIX"
    echo
    print_info "主要输出文件:"
    
    if [ -f "$FINAL_BEST_SEQUENCES" ]; then
        print_info "  ð§¬ $(basename "$FINAL_BEST_SEQUENCES") - 最终最佳序列"
    fi

    if [ -f "$COMPATIBILITY_OUTPUT" ]; then
        print_info "  ð $(basename "$COMPATIBILITY_OUTPUT") - 兼容格式输出"
    fi

    if [ -f "$FINAL_REMOVED_SEQUENCES" ]; then
        print_info "  ð️ $(basename "$FINAL_REMOVED_SEQUENCES") - 被移除的序列"
    fi

    local report_file="${REPORTS_DIR}/${ACTUAL_OUTPUT_PREFIX}_pipeline_report.txt"
    if [ -f "$report_file" ]; then
        print_info "  ð $(basename "$report_file") - 详细报告"
    fi

    local manifest_file="${REPORTS_DIR}/${ACTUAL_OUTPUT_PREFIX}_file_manifest.txt"
    if [ -f "$manifest_file" ]; then
        print_info "  ð $(basename "$manifest_file") - 文件清单"
    fi
    
    if [ -f "$LOG_FILE" ]; then
        print_info "  ð $(basename "$LOG_FILE") - 处理日志"
    fi

    local archive_file="${RESULTS_DIR}/${ACTUAL_OUTPUT_PREFIX}_analysis.tar.gz"
    if [ -f "$archive_file" ]; then
        print_info "  ð¦ $(basename "$archive_file") - 压缩归档"
    fi

    if [ "$CREATE_SUBDIRS" = true ]; then
        echo
        print_info "目录结构:"
        print_info "  ð $OUTPUT_DIR/"
        print_info "    ├── results/     - 最终结果文件"
        print_info "    ├── reports/     - 报告和配置"
        print_info "    ├── logs/        - 处理日志"
        if [ "$KEEP_INTERMEDIATE" = true ]; then
            print_info "    ├── intermediate/ - 中间处理文件"
        fi
        print_info "    └── *.fasta      - 便捷访问链接"
    fi

    echo
    print_success "流水线执行成功！"
    
    if [ -n "$SAVE_CONFIG" ]; then
        save_config "$SAVE_CONFIG"
    fi
}

# 主函数
main() {
    # 显示开始信息
    echo -e "${CYAN}"
    echo "========================================================"
    echo "    基因家族序列处理流水线 v1.3"
    echo "    Author: Biols9527"
    echo "    Date: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "========================================================"
    echo -e "${NC}"

    # 解析参数
    parse_arguments "$@"
    
    # 加载配置文件（如果指定）
    if [ -n "$CONFIG_FILE" ]; then
        load_config "$CONFIG_FILE"
    fi

    # 验证参数
    validate_parameters
    
    # 解析脚本路径
    resolve_script_paths

    # 初始化输出路径
    initialize_output_paths

    # 显示参数信息
    write_log "INFO" "=== 运行参数总览 ==="
    write_log "INFO" "输入文件: $INPUT_FILE"
    write_log "INFO" "输出目录: $OUTPUT_DIR"
    write_log "INFO" "输出前缀: $ACTUAL_OUTPUT_PREFIX"
    write_log "INFO" "序列类型: $SEQUENCE_TYPE"
    write_log "INFO" "线程数: $THREADS"
    write_log "INFO" "BLAST E-value: $EVALUE"
    write_log "INFO" "校正方法: $CORRECTION_METHOD"
    write_log "INFO" "Python命令: $PYTHON_CMD"
    write_log "INFO" "Step1脚本: $STEP1_SCRIPT"
    write_log "INFO" "Step2脚本: $STEP2_SCRIPT"
    write_log "INFO" "日志级别: $LOG_LEVEL"
    
    if [ "$SKIP_CORRECTION" = true ]; then
        write_log "WARNING" "将跳过序列校正步骤"
    fi
    if [ "$SKIP_QUALITY" = true ]; then
        write_log "WARNING" "将跳过质量分析步骤"
    fi
    if [ "$DRY_RUN" = true ]; then
        write_log "WARNING" "干运行模式 - 不会实际执行命令"
    fi

    # 检查依赖
    check_dependencies

    # 分析输入
    analyze_input

    # 记录开始时间
    START_TIME=$(date +%s)

    # 执行主要步骤
    run_sequence_correction
    run_quality_analysis
    organize_output

    # 生成报告
    generate_pipeline_report

    # 创建便捷链接
    create_convenience_links

    # 压缩文件
    compress_outputs

    # 清理中间文件
    cleanup_intermediate

    # 计算运行时间
    END_TIME=$(date +%s)
    RUNTIME=$((END_TIME - START_TIME))

    write_log "INFO" "总运行时间: ${RUNTIME} 秒"

    # 显示最终结果
    show_final_results
}

# 错误处理
trap 'write_log "ERROR" "脚本异常退出，行号: $LINENO"' ERR

# 运行主函数
main "$@"
