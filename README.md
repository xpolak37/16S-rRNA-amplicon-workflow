# 16S rRNA amplicon workflow

This repository includes the bioinformatic workflow for short-read 16S rRNA amplicon processing, all parameters are set to 
work properly on specific data from our lab. Before any usage, please check if it fits your data as well. 

## 🚀 Quick usage

Whole workflow can be executed by main_metagenomics.sh, where you have to change paths to conda environments and to your data

```bash
bash main_16s.sh
```

### 📥 Inputs

### 📤 Outputs

### 📦 Requirements
primers.fasta
phix bowtie indexes
fastqc
multiqc
hostile
cutdapt
qiime
vsearch
pre-trained taxonomic classifier

## ℹ️ About
The workflow consists of multiple modules, all will be executed if not set otherwise:

---
### 🧹 **<u>1. Quality control</u>**

Runs **FastQC & MultiQC** to provide quality control report of all samples. Additionally, it creates a **QUICK SUMMARY** 
of adapter content, sequencing depth, overrepresented sequences and their BLAST hits, for overall summary of the run. Lastly, it runs a primer check, generating a report on the number of occurrences of specific primer sequences in each sample. This is particularly useful when using phased staggered primers and you want to verify whether their distribution is relatively uniform.

**Input**: raw fastqs  
**Output**:
- *quality_raw* folder with fastqc and multiqc data
- *primers_check* folder with forward and reverse primers summary
- *run_info* folder with:
     - raw_custom_summary.txt - adapter content, sequencing depth, overrepresented sequences
     - raw_multiqc_report.html - copy of the multiqc_report.html
     - reads_count_summary.txt - read counts
---    


### ⚙️ **<u>2. Preprocessing</u>**

Runs **cutadapt** for preprocessing the raw fastq files. It trims primers and Nextera Transposase Adapters. Generally, it is executed with the default parameters with the minor changes. See the command below.

```bash
cutadapt \
    -g <forward_primer_sequence> \
    -G <reverse_primer_sequence> \
    -a TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -A GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
    -a 'A{10}' -a 'G{10}' -g 'A{10}' -g 'G{10} \
    --discard-untrimmed \
    -o /path/to/output/sample_R1_trimmed.fastq.gz \
    -O /path/to/output/sample_R2_trimmed.fastq.gz

```

Next, it runs **hostile** to remove human and phix contamination. As for human decontamination, it uses hostile's prepared indexes (human-t2t-hla-argos985-mycob140) that consist of T2T-CHM13v2.0 + IPD-IMGT/HLA v3.51 masked with 150mers for 985 FDA-ARGOS bacterial & 140 mycobacterial genomes. As for phix decontamination, it uses custom indexes of [phiX174](https://www.ncbi.nlm.nih.gov/nuccore/9626372), built by:

```bash
bowtie2-build phiX174.fasta phiX174
```

See the executed commands below. 

```bash
# human decontamination
hostile clean \
--fastq1 sample_R1_trimmed.fastq.gz \
--fastq2 sample_R2_trimmed.fastq.gz \
--index human-t2t-hla-argos985-mycob140
--output /path/to/project_dir/decontaminated/human/

# phix decontamination
hostile clean \
--fastq1 sample_R1_trimmed.clean_1.fastq.gz \
--fastq2 sample_R2_trimmed.clean_2.fastq.gz \
--index path/to/prebuilt/phix_indexes \
--output /path/to/project_dir/decontaminated/human_phix/
```

**Input**: raw fastqs  
**Output**: 
- *trimmed* folder - results of trimmomatic
- *decontaminated* folder:
    - *decontaminated/human/* folder - results of hostile's human removal
    - *decontmainated/human_phix/* folder - results of hostile's phix removal
- *run_info* folder with:
     - read_counts_summary.txt - updated read counts track with <u>trimmed</u> and <u>decontaminated</u> columns. 
     - trimmed_multiqc_report.html - copy of the multiqc_report.html containing fastqc reports of trimmed reads
---

### 🧬 **<u> 3. Denoising</u>**

Performs denosing using one of two (or both) tools: **DADA2** and/or **Deblur** through QIIME2 plugins. Before denoising, it runs bbmerge for merging paired-end reads and vsearch for reorienting reads to same orientation.

- **DADA2** uses a model-based approach to correct sequencing errors by learning error rates from the data itself and distinguishing true biological sequences (amplicon sequence variants, ASVs) from errors. It performs quality filtering, error modeling, dereplication, and chimera removal, providing single-nucleotide resolution of microbial variants.

- **Deblur** instead applies a static error model based on known Illumina error profiles and removes reads that are likely to be erroneous by comparing them to more abundant sequences. Unlike DADA2, it does not learn error rates per dataset, making it faster and more consistent across runs but potentially less sensitive to dataset-specific error patterns.

*Which to use*:

This pipeline is designed to offer two separate approaches that are independent of each other. For an objective selection of a suitable tool, see some of the benchmark studies, for example [Nearing et al. 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6087418/). When not specified by user, the pipeline executes DADA2.  

We offer a brief overview and comparison of DADA2 and Deblur (generated by Chatgpt 5) below.


***Summary Comparison***

| Feature                  | **DADA2**                                                  | **Deblur**                                                 |
|---------------------------|-------------------------------------------------------------|-------------------------------------------------------------|
| **Methodology**           | Learns and models dataset-specific error rates to infer exact amplicon sequence variants (ASVs). | Uses a pre-defined static error model to remove likely erroneous reads. |
| **Accuracy**              | Very high; adapts to dataset-specific sequencing errors for precise ASV inference. | High but may miss dataset-specific nuances in error patterns. |
| **Abundance Estimation**  | Sensitive and accurate, distinguishing true biological variants from errors. | Consistent across runs but may underrepresent rare variants. |
| **Speed / Resources**     | Slower and more computationally intensive due to error modeling. | Faster and less resource-demanding, suitable for large datasets. |
| **Reproducibility**       | Varies slightly across datasets due to adaptive modeling. | Highly reproducible since the error model is fixed. |
| **Best Use**              | When maximum accuracy and error modeling are critical (e.g., high-quality or smaller datasets). | When consistency, scalability, and speed across many datasets are priorities. |


### 🦠 **<u> 4. Taxonomic assignment</u>**

Right now, **idtaxa()** from DECIPHER R package is implemented in this pipeline. It can work both with DADA2 or Deblur outputs. 

**IdTaxa()** is a function for taxonomic classification of DNA sequences (typically 16S, 18S, or ITS) using a machine-learning approach based on discriminant analysis of sequence k-mer composition. It assigns taxonomy by comparing query sequences to a trained reference database, estimating confidence scores for each taxonomic rank.

👀 Comming soon 👀 

The assigntaxa(), QIIME classifier will be added in near future. 

### 🔬 **5. Functional profiling**     

👀 Comming soon 👀


