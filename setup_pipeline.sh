#!/bin/bash

#===============================================================================
# 16S PROFILING PIPELINE - SETUP SCRIPT
#===============================================================================
# This script downloads and sets up all required databases and containers
# for the metagenomic profiling pipeline.
#
# Usage:
#   ./setup_pipeline.sh [INSTALLATION_DIRECTORY]
#
# Example:
#   ./setup_pipeline.sh /shared/pipeline_resources
#   ./setup_pipeline.sh  (interactive mode - will prompt for directory)
#
# Author: Petra Polakovicova 
# Version: 1.0.0
#===============================================================================

set -e  # Exit on error
set -u  # Exit on undefined variable

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Script start time
START_TIME=$(date +%s)

#===============================================================================
# FUNCTIONS
#===============================================================================

print_header() {
    echo -e "${BLUE}======================================================================${NC}"
    echo -e "${BLUE}  16S PROFILING PIPELINE - SETUP${NC}"
    echo -e "${BLUE}======================================================================${NC}"
    echo ""
}

log_info() {
    if [ -n "${LOGFILE:-}" ]; then
        echo -e "${GREEN}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "${LOGFILE}"
    else
        echo -e "${GREEN}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
    fi
}

log_warn() {
    if [ -n "${LOGFILE:-}" ]; then
        echo -e "${YELLOW}[WARN]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "${LOGFILE}"
    else
        echo -e "${YELLOW}[WARN]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
    fi
}

log_error() {
    if [ -n "${LOGFILE:-}" ]; then
        echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "${LOGFILE}"
    else
        echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
    fi
}

log_success() {
    if [ -n "${LOGFILE:-}" ]; then
        echo -e "${GREEN}[SUCCESS]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "${LOGFILE}"
    else
        echo -e "${GREEN}[SUCCESS]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
    fi
}

check_command() {
    if ! command -v $1 &> /dev/null; then
        log_error "$1 is not installed or not in PATH"
        exit 1
    fi
}

#===============================================================================
# GET INSTALLATION DIRECTORY
#===============================================================================

print_header

if [ $# -eq 0 ]; then
    # Interactive mode
    echo -e "${YELLOW}No installation directory provided.${NC}"
    echo ""
    echo "Please enter the full path where you want to install pipeline resources:"
    echo "(This will create subdirectories: databases/ and singularity_cache/)"
    echo ""
    read -p "Installation directory: " INSTALL_DIR
    
    # Trim whitespace
    INSTALL_DIR=$(echo "$INSTALL_DIR" | xargs)
    
    if [ -z "$INSTALL_DIR" ]; then
        log_error "No directory provided. Exiting."
        exit 1
    fi
else
    # Command line argument
    INSTALL_DIR="$1"
fi

# Convert to absolute path
INSTALL_DIR=$(realpath -m "$INSTALL_DIR")

log_info "Installation directory: ${INSTALL_DIR}"
echo ""

# Confirm with user
read -p "Continue with this directory? (yes/no): " CONFIRM
if [[ ! "$CONFIRM" =~ ^[Yy][Ee][Ss]$|^[Yy]$ ]]; then
    log_warn "Installation cancelled by user."
    exit 0
fi

#===============================================================================
# SETUP DIRECTORIES
#===============================================================================

echo ""
log_info "Setting up directory structure..."

# Check current directory
PIPELINE_DIR=$(pwd)
PIPELINE_DIR=$(realpath -m "$PIPELINE_DIR")
echo ${PIPELINE_DIR}

# Create main directory
mkdir -p "${INSTALL_DIR}"
cd "${INSTALL_DIR}"

# Create subdirectories
mkdir -p classifiers
mkdir -p singularity_cache
mkdir -p logs
mkdir -p bowtie_phix

# Setup log file
LOGFILE="${INSTALL_DIR}/logs/setup_$(date +%Y%m%d_%H%M%S).log"
touch "${LOGFILE}"

log_success "Directory structure created:"
log_info "  - ${INSTALL_DIR}/classifiers"
log_info "  - ${INSTALL_DIR}/singularity_cache"
log_info "  - ${INSTALL_DIR}/logs"
log_info "Log file: ${LOGFILE}"

#===============================================================================
# CHECK REQUIRED TOOLS
#===============================================================================

echo ""
log_info "Checking required tools..."

check_command wget
check_command tar
check_command singularity

log_success "All required tools are available"

#===============================================================================
# DOWNLOAD CLASSIFIERS
#===============================================================================

echo ""
log_info "========================================="
log_info "STEP 1/3: Downloading pre-built taxonomic classifiers"
log_info "========================================="

CLASSIFIERS_DIR="${INSTALL_DIR}/classifiers"
cd "${CLASSIFIERS_DIR}"

if [ -f "silva-138.2-ssu-nr99-341F-805R-classifier.qza" ]; then
    log_warn "Classifiers already present. Skipping download."
else
    log_info "Downloading classifiers from Zenodo..."
    
    wget "https://filesender.cesnet.cz/download.php?token=a1fbed9f-128d-4353-ab69-62d15298aa2d&files_ids=831961" \
        -O classifiers.tar.gz 2>&1 | tee -a "${LOGFILE}"
    
    if [ $? -eq 0 ]; then
        log_success "Download completed"
        
        log_info "Extracting folder..."
        tar -xzvf classifiers.tar.gz >> "${LOGFILE}" 2>&1
        
        if [ -f "silva-138.2-ssu-nr99-341F-805R-classifier.qza" ]; then
            log_success "Classifiers extracted successfully"
            
            # Cleanup
            log_info "Removing tarball to save space..."
            rm classifiers.tar.gz
        else
            log_error "Database extraction failed"
            exit 1
        fi
    else
        log_error "Failed to download classifiers"
        exit 1
    fi
fi

#===============================================================================
# DOWNLOAD SINGULARITY CONTAINERS
#===============================================================================

echo ""
log_info "========================================="
log_info "STEP 2/4: Downloading Singularity containers"
log_info "========================================="

SING_DIR="${INSTALL_DIR}/singularity_cache"
cd "${SING_DIR}"

declare -A CONTAINERS=(
    ["quay.io-biocontainers-bbmap-39.52--he5f24ec_0.img"]="docker://quay.io/biocontainers/bbmap:39.52--he5f24ec_0"
    ["quay.io-biocontainers-bioconductor-dada2-1.38.0--r45ha27e39d_0.img"]="docker://quay.io/biocontainers/bioconductor-dada2:1.38.0--r45ha27e39d_0"
    ["quay.io-biocontainers-bioconductor-decipher-3.6.0--r45h01b2380_0.img"]="docker://quay.io/biocontainers/bioconductor-decipher:3.6.0--r45h01b2380_0"
    ["quay.io-biocontainers-cutadapt-5.2--py311haab0aaa_0.img"]="docker://quay.io/biocontainers/cutadapt:5.2--py311haab0aaa_0"
    ["quay.io-biocontainers-fastqc-0.12.1--hdfd78af_0.img"]="docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    ["quay.io-biocontainers-multiqc-1.21--pyhdfd78af_0.img"]="docker://quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0"
    ["quay.io-biocontainers-pandas-2.2.1.img"]="docker://quay.io/biocontainers/pandas:2.2.1"
    ["quay.io-qiime2-amplicon-2026.1.img"]="docker://quay.io/qiime2/amplicon:2026.1"
    ["quay.io-biocontainers-entrez-direct-24.0--he881be0_0.img"]="docker://quay.io/biocontainers/entrez-direct:24.0--he881be0_0"
    ["quay.io-biocontainers-bowtie2-2.5.5--ha27dd3b_0.img"]="docker://quay.io/biocontainers/bowtie2:2.5.5--ha27dd3b_0"
)
CONTAINER_COUNT=0
TOTAL_CONTAINERS=${#CONTAINERS[@]}

for img_name in "${!CONTAINERS[@]}"; do
    CONTAINER_COUNT=$((CONTAINER_COUNT + 1))
    
    if [ -f "${img_name}" ]; then
        log_warn "[${CONTAINER_COUNT}/${TOTAL_CONTAINERS}] ${img_name} already exists. Skipping."
    else
        log_info "[${CONTAINER_COUNT}/${TOTAL_CONTAINERS}] Pulling ${img_name}..."
        
        uri="${CONTAINERS[$img_name]}"
        
        singularity pull --name "${img_name}" "${uri}" >> "${LOGFILE}" 2>&1
        
        if [ $? -eq 0 ]; then
            log_success "[${CONTAINER_COUNT}/${TOTAL_CONTAINERS}] ${img_name} downloaded successfully"
        else
            log_error "[${CONTAINER_COUNT}/${TOTAL_CONTAINERS}] Failed to download ${img_name}"
            exit 1
        fi
    fi
done

#===============================================================================
# BUILDING BOWTIE INDEXES
#===============================================================================

echo ""
log_info "========================================="
log_info "STEP 3/4: Building bowtie indexes"
log_info "========================================="

BOWTIE_DIR="${INSTALL_DIR}/bowtie_phix"
cd "${BOWTIE_DIR}"

if [ -f "phiX174.1.bt2" ]; then
    log_warn "Bowtie indexes already present. Skipping."
else
    log_info "Downloading the phiX174 genome."
    singularity exec ${SING_DIR}/quay.io-biocontainers-entrez-direct-24.0--he881be0_0.img efetch -db nucleotide -id NC_001422.1 -format fasta > phiX174.fasta
    
    if [ $? -eq 0 ]; then
        log_success "Download completed"
        
        log_info "Building bowtie indexes for phiX174 genome..."

        singularity exec ${SING_DIR}/quay.io-biocontainers-bowtie2-2.5.5--ha27dd3b_0.img bowtie2-build phiX174.fasta phiX174

        if [ -f "phiX174.1.bt2" ]; then
            log_success "Indexes built succesfully"
        
        else
            log_error "Process failed"
            exit 1
        fi
    else
        log_error "Process failed"
        exit 1
    fi
fi

#===============================================================================
# DOWNLOAD TESTING DATASET
#===============================================================================

echo ""
log_info "========================================="
log_info "STEP 4/4: Downloading test dataset"
log_info "========================================="

TEST_DIR="${PIPELINE_DIR}"
cd "${TEST_DIR}"

if [ -f "test/SRR13005968-test_2.fastq.gz" ]; then
    log_warn "Test data already present. Skipping download."
else
    log_info "Downloading test data from Zenodo..."
    
    wget "https://filesender.cesnet.cz/download.php?token=323bf327-a8a5-4ec2-b82c-0c9026388202&files_ids=836539" \
        -O test.zip 2>&1 | tee -a "${LOGFILE}"
    
    if [ $? -eq 0 ]; then
        log_success "Download completed"
        
        log_info "Extracting folder..."
        unzip test.zip >> "${LOGFILE}" 2>&1
        
        if [ -f "test/SRR13005968-test_2.fastq.gz" ]; then
            log_success "Test data extracted successfully"
            
            # Cleanup
            log_info "Removing zip to save space..."
            rm test.zip
        else
            log_error "Data extraction failed"
            exit 1
        fi
    else
        log_error "Failed to download test data"
        exit 1
    fi
fi


#===============================================================================
# GENERATE CONFIGURATION FILE
#===============================================================================

echo ""
log_info "========================================="
log_info "Generating pipeline configuration"
log_info "========================================="

CONFIG_FILE="${INSTALL_DIR}/pipeline_paths.config"

cat > "${CONFIG_FILE}" << EOF
/*
========================================================================================
    Pipeline Resource Paths Configuration
========================================================================================
    Generated by setup_pipeline.sh on $(date)
    
    Use this configuration with:
    nextflow run main.nf -c ${CONFIG_FILE} [other options]
========================================================================================
*/

params {
    // Cache directories
    singularity_cache_dir = '${INSTALL_DIR}/singularity_cache'
    classifiers_dir    = '${INSTALL_DIR}/classifiers'
}
EOF

log_success "Configuration file created: ${CONFIG_FILE}"

#===============================================================================
# SUMMARY
#===============================================================================

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
MINUTES=$((DURATION / 60))
SECONDS=$((DURATION % 60))

echo ""
echo -e "${GREEN}======================================================================${NC}"
echo -e "${GREEN}  SETUP COMPLETE!${NC}"
echo -e "${GREEN}======================================================================${NC}"
echo ""
log_success "All resources downloaded and installed successfully"
echo ""
echo "Installation Summary:"
echo "  - Installation directory: ${INSTALL_DIR}"
echo "  - Total time: ${MINUTES} minutes ${SECONDS} seconds"
echo ""
echo "Resource Locations:"
echo "  - Classifiers:     ${INSTALL_DIR}/classifiers/"
echo "  - Containers:         ${INSTALL_DIR}/singularity_cache/"
echo "  - Configuration:      ${CONFIG_FILE}"
echo "  - Log file:           ${LOGFILE}"
echo ""
echo "Next Steps:"
echo "  1. Run the pipeline with:"
echo "     nextflow run main.nf -c ${CONFIG_FILE} --input samples.csv --outdir results"
echo ""
echo "  2. Or manually specify paths:"
echo "     nextflow run main.nf \\"
echo "       --singularity_cache_dir ${INSTALL_DIR}/singularity_cache \\"
echo "       --classifiers ${INSTALL_DIR}/classifiers \\"
echo "       --input samples.csv --outdir results"
echo ""
echo -e "${GREEN}======================================================================${NC}"