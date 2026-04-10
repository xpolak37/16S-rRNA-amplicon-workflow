#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================
    16S PROFILING PIPELINE
========================================================================================
    FastQC -> Cutadapt -> DADA2 + QIIME- DEBLUR + USEARCH -> Naive bayes + BLAST + IDTAXA + assignTaxonomy + SINTAX -> MultiQC
----------------------------------------------------------------------------------------
*/

// Print pipeline header
log.info """\
    ===================================
    16S PROFILING PIPELINE
    ===================================
    Input samplesheet : ${params.input}
    Output directory  : ${params.outdir}
    ===================================
    """
    .stripIndent()

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { FASTQC as FASTQC_RAW }        from './modules/fastqc'
include { FASTQC as FASTQC_TRIMMED }    from './modules/fastqc'
include { MULTIQC }                     from './modules/multiqc'
include { CUSTOM_SUMMARY_PARSE; CUSTOM_SUMMARY_BLAST; CUSTOM_SUMMARY_RENDER } from './modules/custom_summary'
include { CUTADAPT; CUTADAPT_DADA2_ORIENT } from './modules/cutadapt'
include { HOST_REMOVAL }                    from './modules/hostile.nf'
include { PHIX_REMOVAL }                    from './modules/hostile.nf'
include { DADA2_PAIRED }                       from './modules/dada2'
include { DADA2_SINGLE }                       from './modules/dada2'
include { MERGING_READS }               from './modules/se_preprocessing.nf'
include { ORIENTING_READS }               from './modules/se_preprocessing.nf'
include { QIIME_IMPORT; QIIME_DEBLUR }  from './modules/deblur.nf'
include { VSEARCH_UNOISE3 }             from './modules/unoise.nf'
include { SEQTK_SUBSAMPLE }             from './modules/seqtk_subsample'
include { QIIME_NAIVE_BAYES as QIIME_NAIVE_BAYES_DADA2_PAIRED }  from './modules/tax_classifiers.nf'
include { QIIME_NAIVE_BAYES as QIIME_NAIVE_BAYES_DADA2_SINGLE }  from './modules/tax_classifiers.nf'
include { QIIME_NAIVE_BAYES as QIIME_NAIVE_BAYES_DEBLUR }  from './modules/tax_classifiers.nf'
include { QIIME_NAIVE_BAYES as QIIME_NAIVE_BAYES_UNOISE }  from './modules/tax_classifiers.nf'
include { QIIME_BLAST as QIIME_BLAST_DADA2_PAIRED }  from './modules/tax_classifiers.nf'
include { QIIME_BLAST as QIIME_BLAST_DADA2_SINGLE }  from './modules/tax_classifiers.nf'
include { QIIME_BLAST as QIIME_BLAST_DEBLUR }  from './modules/tax_classifiers.nf'
include { QIIME_BLAST as QIIME_BLAST_UNOISE }  from './modules/tax_classifiers.nf'
include { IDTAXA as IDTAXA_DADA2_PAIRED }  from './modules/tax_classifiers.nf'
include { IDTAXA as IDTAXA_DADA2_SINGLE }  from './modules/tax_classifiers.nf'
include { IDTAXA as IDTAXA_DEBLUR }  from './modules/tax_classifiers.nf'
include { IDTAXA as IDTAXA_UNOISE }  from './modules/tax_classifiers.nf'
include { ASSIGNTAXONOMY as ASSIGNTAXONOMY_DADA2_PAIRED }  from './modules/tax_classifiers.nf'
include { ASSIGNTAXONOMY as ASSIGNTAXONOMY_DADA2_SINGLE }  from './modules/tax_classifiers.nf'
include { ASSIGNTAXONOMY as ASSIGNTAXONOMY_DEBLUR }  from './modules/tax_classifiers.nf'
include { ASSIGNTAXONOMY as ASSIGNTAXONOMY_UNOISE }  from './modules/tax_classifiers.nf'
include { METASTANDARD as NB_DADA2_PAIRED_METASTANDARD }        from './modules/MetaStandard16S'
include { METASTANDARD as NB_DADA2_SINGLE_METASTANDARD }        from './modules/MetaStandard16S'
include { METASTANDARD as NB_DEBLUR_METASTANDARD }       from './modules/MetaStandard16S'
include { METASTANDARD as NB_UNOISE_METASTANDARD }       from './modules/MetaStandard16S'
include { METASTANDARD as BLAST_DADA2_PAIRED_METASTANDARD }        from './modules/MetaStandard16S'
include { METASTANDARD as BLAST_DADA2_SINGLE_METASTANDARD }        from './modules/MetaStandard16S'
include { METASTANDARD as BLAST_DEBLUR_METASTANDARD }       from './modules/MetaStandard16S'
include { METASTANDARD as BLAST_UNOISE_METASTANDARD }       from './modules/MetaStandard16S'
include { METASTANDARD as IDTAXA_DADA2_PAIRED_METASTANDARD }        from './modules/MetaStandard16S'
include { METASTANDARD as IDTAXA_DADA2_SINGLE_METASTANDARD }        from './modules/MetaStandard16S'
include { METASTANDARD as IDTAXA_DEBLUR_METASTANDARD }       from './modules/MetaStandard16S'
include { METASTANDARD as IDTAXA_UNOISE_METASTANDARD }       from './modules/MetaStandard16S'
include { METASTANDARD as AT_DADA2_PAIRED_METASTANDARD }        from './modules/MetaStandard16S'
include { METASTANDARD as AT_DADA2_SINGLE_METASTANDARD }        from './modules/MetaStandard16S'
include { METASTANDARD as AT_DEBLUR_METASTANDARD }       from './modules/MetaStandard16S'
include { METASTANDARD as AT_UNOISE_METASTANDARD }       from './modules/MetaStandard16S'

// Mock community evaluation
include { MOCK_EVALUATION as NB_DADA2_PAIRED_MOCK }        from './modules/mock_evaluation'
include { MOCK_EVALUATION as NB_DADA2_SINGLE_MOCK }        from './modules/mock_evaluation'
include { MOCK_EVALUATION as NB_DEBLUR_MOCK }              from './modules/mock_evaluation'
include { MOCK_EVALUATION as BLAST_DADA2_PAIRED_MOCK }     from './modules/mock_evaluation'
include { MOCK_EVALUATION as BLAST_DADA2_SINGLE_MOCK }     from './modules/mock_evaluation'
include { MOCK_EVALUATION as BLAST_DEBLUR_MOCK }           from './modules/mock_evaluation'
include { MOCK_EVALUATION as IDTAXA_DADA2_PAIRED_MOCK }    from './modules/mock_evaluation'
include { MOCK_EVALUATION as IDTAXA_DADA2_SINGLE_MOCK }    from './modules/mock_evaluation'
include { MOCK_EVALUATION as IDTAXA_DEBLUR_MOCK }          from './modules/mock_evaluation'
include { MOCK_EVALUATION as AT_DADA2_PAIRED_MOCK }        from './modules/mock_evaluation'
include { MOCK_EVALUATION as AT_DADA2_SINGLE_MOCK }        from './modules/mock_evaluation'
include { MOCK_EVALUATION as AT_DEBLUR_MOCK }              from './modules/mock_evaluation'
include { MOCK_EVALUATION as NB_UNOISE_MOCK }              from './modules/mock_evaluation'
include { MOCK_EVALUATION as BLAST_UNOISE_MOCK }           from './modules/mock_evaluation'
include { MOCK_EVALUATION as IDTAXA_UNOISE_MOCK }          from './modules/mock_evaluation'
include { MOCK_EVALUATION as AT_UNOISE_MOCK }              from './modules/mock_evaluation'

// ============================================================
// CONSTANTS — all valid tool names per axis
// ============================================================

def VALID_ASV        = ['dada2_paired','dada2_single', 'deblur', 'unoise']
def VALID_CLASSIFIER = ['qnb', 'qblast',"idtaxa","assigntaxonomy"]

def resolveWorkflows(List validASV, List validClassifier) {
    if (params.all) {
        log.info "Mode: --all  →  running every tool on both axes"
        return [
            asv_tools   : validASV        as LinkedHashSet,
            classifiers : validClassifier as LinkedHashSet
        ]
    }

    def asvTools    = resolveTools(params.asv_inference as String, validASV,        'asv_inference')
    def classifiers = resolveTools(params.classifier    as String, validClassifier,  'classifier')

    log.info "Mode: selected  →  ASV: ${asvTools}  |  Classifiers: ${classifiers}"
    return [
        asv_tools   : asvTools,
        classifiers : classifiers
    ]
}
def workflowsToRun = resolveWorkflows(VALID_ASV, VALID_CLASSIFIER)

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {
     
    // Read and parse samplesheet
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> 
            def sample_id = row.sample
            def read1 = file(row.read1)
            def read2 = file(row.read2)
            
            // Validate files exist
            if (!read1.exists()) exit 1, "ERROR: Read1 file does not exist: ${read1}"
            if (!read2.exists()) exit 1, "ERROR: Read2 file does not exist: ${read2}"
            
            return tuple(sample_id, read1, read2)
        }
        .set { ch_input_reads }
    
    // FastQC on raw reads
    FASTQC_RAW(ch_input_reads, "raw")
    
    // cutadapt for primer and adapter trimming
    CUTADAPT(ch_input_reads)

    // fastqc TRIMMED
    FASTQC_TRIMMED(CUTADAPT.out.reads, "trimmed")
    
    // Host and PhiX removal
    path_bowtie_phix = params.bowtie_dir
    path_hostile_index = params.hostile_index_dir
    HOST_REMOVAL(CUTADAPT.out.reads, path_hostile_index)
    PHIX_REMOVAL(HOST_REMOVAL.out.reads, path_bowtie_phix)

    // Optional subsampling for --quick mode
    if (params.quick) {
        SEQTK_SUBSAMPLE(PHIX_REMOVAL.out.reads)
        ch_reads_for_asv = SEQTK_SUBSAMPLE.out.reads
    } else {
        ch_reads_for_asv = PHIX_REMOVAL.out.reads
    }

    // DADA2 PAIRED
    if ('dada2_paired' in workflowsToRun.asv_tools) {
        // Split each sample's reads by orientation so DADA2 learns separate
        // error models per orientation, then merges the ASV tables afterwards.
        CUTADAPT_DADA2_ORIENT(ch_reads_for_asv)

        dada2_input_paired = CUTADAPT_DADA2_ORIENT.out
            .map { sample_id, fwd_r1, fwd_r2, rev_r1, rev_r2 ->
                   [fwd_r1, fwd_r2, rev_r1, rev_r2] }
            .flatten()
            .collect()

        DADA2_PAIRED(dada2_input_paired)
    }

    // Merging (only for tools that need merged reads) ────
    def needs_merging = ['deblur', 'unoise','dada2_single']

    if (workflowsToRun.asv_tools.any { it in needs_merging }) {
        MERGING_READS(ch_reads_for_asv)
    }

    // Orienting (only for tools that can have oriented reads)
    def needs_orienting = ['deblur', 'unoise']

    if (workflowsToRun.asv_tools.any { it in needs_merging }) {
        // TO DO - new orienting db
        orienting_db=file(params.classifiers_dir + "/silva_nrr99_v138.1_orienting.udb")
        ORIENTING_READS(MERGING_READS.out.reads,orienting_db)
    }

    // DADA2 SINGLE
    if ('dada2_single' in workflowsToRun.asv_tools) {
        dada2_input_single = MERGING_READS.out.reads
        .map { sample_id, reads -> reads }
        .collect()
        
        DADA2_SINGLE(dada2_input_single)
    }

    // DEBLUR 
    if ('deblur' in workflowsToRun.asv_tools) {

        qiime_input = ORIENTING_READS.out.reads
        .map { sample_id, reads -> reads }
        .collect()
        
        QIIME_IMPORT(qiime_input)
        QIIME_DEBLUR(QIIME_IMPORT.out)
    }

    // UNOISE3 (vsearch)
    if ('unoise' in workflowsToRun.asv_tools) {
        unoise_input = ORIENTING_READS.out.reads
        .map { sample_id, reads -> reads }
        .collect()

        VSEARCH_UNOISE3(unoise_input)
    }

    // TAXONOMIC ASSIGNMENT + METASTANDARD + MOCK EVALUATION
    
    // Reference files for mock evaluation
    mock_abundance_file = params.mock_evaluation ? file(params.mock_abundance) : null
    mock_taxa_file = params.mock_evaluation ? file(params.mock_taxa) : null
    mock_synonyms_file = params.mock_evaluation ? file(params.mock_synonyms) : null

    // QIIME NAIVE BAYES
    if ('qnb' in workflowsToRun.classifiers) {

        qiime_NB_classifier = file(params.classifiers_dir + "/silva-138.2-ssu-nr99-341F-805R-classifier.qza")

        // DEBLUR
        if ('deblur' in workflowsToRun.asv_tools) {
            QIIME_NAIVE_BAYES_DEBLUR(QIIME_DEBLUR.out.fasta,"deblur",qiime_NB_classifier)
            NB_DEBLUR_METASTANDARD(QIIME_DEBLUR.out.asv_table,QIIME_NAIVE_BAYES_DEBLUR.out,"NaiveBayes")
            if (params.mock_evaluation) {
                NB_DEBLUR_MOCK(NB_DEBLUR_METASTANDARD.out, mock_abundance_file, mock_taxa_file, mock_synonyms_file, "deblur_NaiveBayes")
            }
        }

        // DADA2
        if ('dada2_paired' in workflowsToRun.asv_tools) {
            QIIME_NAIVE_BAYES_DADA2_PAIRED(DADA2_PAIRED.out.fasta,"dada2_paired",qiime_NB_classifier)
            NB_DADA2_PAIRED_METASTANDARD(DADA2_PAIRED.out.asv_table,QIIME_NAIVE_BAYES_DADA2_PAIRED.out,"NaiveBayes")
            if (params.mock_evaluation) {
                NB_DADA2_PAIRED_MOCK(NB_DADA2_PAIRED_METASTANDARD.out, mock_abundance_file, mock_taxa_file, mock_synonyms_file, "dada2PE_NaiveBayes")
            }
        } 

        if ('dada2_single' in workflowsToRun.asv_tools) {
            QIIME_NAIVE_BAYES_DADA2_SINGLE(DADA2_SINGLE.out.fasta,"dada2_single",qiime_NB_classifier)
            NB_DADA2_SINGLE_METASTANDARD(DADA2_SINGLE.out.asv_table,QIIME_NAIVE_BAYES_DADA2_SINGLE.out,"NaiveBayes")
            if (params.mock_evaluation) {
                NB_DADA2_SINGLE_MOCK(NB_DADA2_SINGLE_METASTANDARD.out, mock_abundance_file, mock_taxa_file, mock_synonyms_file, "dada2SE_NaiveBayes")
            }
        }

        // UNOISE
        if ('unoise' in workflowsToRun.asv_tools) {
            QIIME_NAIVE_BAYES_UNOISE(VSEARCH_UNOISE3.out.fasta,"unoise",qiime_NB_classifier)
            NB_UNOISE_METASTANDARD(VSEARCH_UNOISE3.out.asv_table,QIIME_NAIVE_BAYES_UNOISE.out,"NaiveBayes")
            if (params.mock_evaluation) {
                NB_UNOISE_MOCK(NB_UNOISE_METASTANDARD.out, mock_abundance_file, mock_taxa_file, mock_synonyms_file, "unoise_NaiveBayes")
            }
        }

    }

    // QIIME BLAST
    if ('qblast' in workflowsToRun.classifiers) {

        qiime_BLAST_classifier_reads = file(params.classifiers_dir + "/silva-138.2-ssu-nr99-seqs-filt.qza")
        qiime_BLAST_classifier_tax = file(params.classifiers_dir + "/silva-138.2-ssu-nr99-tax.qza")

        // DEBLUR
        if ('deblur' in workflowsToRun.asv_tools) {
            QIIME_BLAST_DEBLUR(QIIME_DEBLUR.out.fasta,"deblur",qiime_BLAST_classifier_reads,qiime_BLAST_classifier_tax)
            BLAST_DEBLUR_METASTANDARD(QIIME_DEBLUR.out.asv_table,QIIME_BLAST_DEBLUR.out,"BLAST")
            if (params.mock_evaluation) {
                BLAST_DEBLUR_MOCK(BLAST_DEBLUR_METASTANDARD.out, mock_abundance_file, mock_taxa_file, mock_synonyms_file, "deblur_BLAST")
            }
        }

        // DADA2
        if ('dada2_paired' in workflowsToRun.asv_tools) {
            QIIME_BLAST_DADA2_PAIRED(DADA2_PAIRED.out.fasta,"dada2_paired", qiime_BLAST_classifier_reads, qiime_BLAST_classifier_tax)
            BLAST_DADA2_PAIRED_METASTANDARD(DADA2_PAIRED.out.asv_table,QIIME_BLAST_DADA2_PAIRED.out,"BLAST")
            if (params.mock_evaluation) {
                BLAST_DADA2_PAIRED_MOCK(BLAST_DADA2_PAIRED_METASTANDARD.out, mock_abundance_file, mock_taxa_file, mock_synonyms_file, "dada2PE_BLAST")
            }
        }

        if ('dada2_single' in workflowsToRun.asv_tools) {
            QIIME_BLAST_DADA2_SINGLE(DADA2_SINGLE.out.fasta,"dada2_single", qiime_BLAST_classifier_reads, qiime_BLAST_classifier_tax)
            BLAST_DADA2_SINGLE_METASTANDARD(DADA2_SINGLE.out.asv_table,QIIME_BLAST_DADA2_SINGLE.out,"BLAST")
            if (params.mock_evaluation) {
                BLAST_DADA2_SINGLE_MOCK(BLAST_DADA2_SINGLE_METASTANDARD.out, mock_abundance_file, mock_taxa_file, mock_synonyms_file, "dada2SE_BLAST")
            }
        }
        // UNOISE
        if ('unoise' in workflowsToRun.asv_tools) {
            QIIME_BLAST_UNOISE(VSEARCH_UNOISE3.out.fasta,"unoise",qiime_BLAST_classifier_reads,qiime_BLAST_classifier_tax)
            BLAST_UNOISE_METASTANDARD(VSEARCH_UNOISE3.out.asv_table,QIIME_BLAST_UNOISE.out,"BLAST")
            if (params.mock_evaluation) {
                BLAST_UNOISE_MOCK(BLAST_UNOISE_METASTANDARD.out, mock_abundance_file, mock_taxa_file, mock_synonyms_file, "unoise_BLAST")
            }
        }

    }

    // IDTAXA
    if ('idtaxa' in workflowsToRun.classifiers) {

        idtaxa_classifier = file(params.classifiers_dir + "/idtaxa_trainingSet_V3V4_silva_138_2.RData")

        // DEBLUR
        if ('deblur' in workflowsToRun.asv_tools) {
            IDTAXA_DEBLUR(QIIME_DEBLUR.out.fasta,"deblur",idtaxa_classifier)
            IDTAXA_DEBLUR_METASTANDARD(QIIME_DEBLUR.out.asv_table,IDTAXA_DEBLUR.out,"IDTAXA")
            if (params.mock_evaluation) {
                IDTAXA_DEBLUR_MOCK(IDTAXA_DEBLUR_METASTANDARD.out, mock_abundance_file, mock_taxa_file, mock_synonyms_file, "deblur_IDTAXA")
            }
        }

        // DADA2
        if ('dada2_paired' in workflowsToRun.asv_tools) {
            IDTAXA_DADA2_PAIRED(DADA2_PAIRED.out.fasta,"dada2_paired",idtaxa_classifier)
            IDTAXA_DADA2_PAIRED_METASTANDARD(DADA2_PAIRED.out.asv_table,IDTAXA_DADA2_PAIRED.out,"IDTAXA")
            if (params.mock_evaluation) {
                IDTAXA_DADA2_PAIRED_MOCK(IDTAXA_DADA2_PAIRED_METASTANDARD.out, mock_abundance_file, mock_taxa_file, mock_synonyms_file, "dada2PE_IDTAXA")
            }
        }

        if ('dada2_single' in workflowsToRun.asv_tools) {
            IDTAXA_DADA2_SINGLE(DADA2_SINGLE.out.fasta,"dada2_single",idtaxa_classifier)
            IDTAXA_DADA2_SINGLE_METASTANDARD(DADA2_SINGLE.out.asv_table,IDTAXA_DADA2_SINGLE.out,"IDTAXA")
            if (params.mock_evaluation) {
                IDTAXA_DADA2_SINGLE_MOCK(IDTAXA_DADA2_SINGLE_METASTANDARD.out, mock_abundance_file, mock_taxa_file, mock_synonyms_file, "dada2SE_IDTAXA")
            }
        }

        // UNOISE
        if ('unoise' in workflowsToRun.asv_tools) {
            IDTAXA_UNOISE(VSEARCH_UNOISE3.out.fasta,"unoise",idtaxa_classifier)
            IDTAXA_UNOISE_METASTANDARD(VSEARCH_UNOISE3.out.asv_table,IDTAXA_UNOISE.out,"IDTAXA")
            if (params.mock_evaluation) {
                IDTAXA_UNOISE_MOCK(IDTAXA_UNOISE_METASTANDARD.out, mock_abundance_file, mock_taxa_file, mock_synonyms_file, "unoise_IDTAXA")
            }
        }

    }

    // ASSIGN TAXONOMY
    if ('assigntaxonomy' in workflowsToRun.classifiers) {

        assigntaxonomy_classifier = file(params.classifiers_dir + "/dada2_trainset_v3v4_silva_nr99_v138_2.fa.gz")

        // DEBLUR
        if ('deblur' in workflowsToRun.asv_tools) {
            ASSIGNTAXONOMY_DEBLUR(QIIME_DEBLUR.out.fasta,"deblur",assigntaxonomy_classifier)
            AT_DEBLUR_METASTANDARD(QIIME_DEBLUR.out.asv_table,ASSIGNTAXONOMY_DEBLUR.out,"AssignTaxonomy")
            if (params.mock_evaluation) {
                AT_DEBLUR_MOCK(AT_DEBLUR_METASTANDARD.out, mock_abundance_file, mock_taxa_file, mock_synonyms_file, "deblur_AssignTaxonomy")
            }
        }

        // DADA2
        if ('dada2_paired' in workflowsToRun.asv_tools) {
            ASSIGNTAXONOMY_DADA2_PAIRED(DADA2_PAIRED.out.fasta,"dada2_paired",assigntaxonomy_classifier)
            AT_DADA2_PAIRED_METASTANDARD(DADA2_PAIRED.out.asv_table,ASSIGNTAXONOMY_DADA2_PAIRED.out,"AssignTaxonomy")
            if (params.mock_evaluation) {
                AT_DADA2_PAIRED_MOCK(AT_DADA2_PAIRED_METASTANDARD.out, mock_abundance_file, mock_taxa_file, mock_synonyms_file, "dada2PE_AssignTaxonomy")
            }
        }

        if ('dada2_single' in workflowsToRun.asv_tools) {
            ASSIGNTAXONOMY_DADA2_SINGLE(DADA2_SINGLE.out.fasta,"dada2_single",assigntaxonomy_classifier)
            AT_DADA2_SINGLE_METASTANDARD(DADA2_SINGLE.out.asv_table,ASSIGNTAXONOMY_DADA2_SINGLE.out,"AssignTaxonomy")
            if (params.mock_evaluation) {
                AT_DADA2_SINGLE_MOCK(AT_DADA2_SINGLE_METASTANDARD.out, mock_abundance_file, mock_taxa_file, mock_synonyms_file, "dada2SE_AssignTaxonomy")
            }
        }

        // UNOISE
        if ('unoise' in workflowsToRun.asv_tools) {
            ASSIGNTAXONOMY_UNOISE(VSEARCH_UNOISE3.out.fasta,"unoise",assigntaxonomy_classifier)
            AT_UNOISE_METASTANDARD(VSEARCH_UNOISE3.out.asv_table,ASSIGNTAXONOMY_UNOISE.out,"AssignTaxonomy")
            if (params.mock_evaluation) {
                AT_UNOISE_MOCK(AT_UNOISE_METASTANDARD.out, mock_abundance_file, mock_taxa_file, mock_synonyms_file, "unoise_AssignTaxonomy")
            }
        }

    }

    // Collect all QC files for MultiQC
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMMED.out.zip.collect().ifEmpty([]))
    
    // MultiQC aggregation
    MULTIQC(ch_multiqc_files.collect())

    // Custom per-run summary report (raw FastQC + optional remote BLAST)
    if (params.custom_summary) {
        CUSTOM_SUMMARY_PARSE(FASTQC_RAW.out.zip.collect())
        if (params.custom_summary_blast) {
            CUSTOM_SUMMARY_BLAST(CUSTOM_SUMMARY_PARSE.out.fasta)
            blast_tsv_ch = CUSTOM_SUMMARY_BLAST.out.tsv
        } else {
            // Provide an empty placeholder so RENDER's input is satisfied
            blast_tsv_ch = Channel.fromPath("${projectDir}/assets/empty_blast_hits.tsv")
        }
        CUSTOM_SUMMARY_RENDER(CUSTOM_SUMMARY_PARSE.out.json, blast_tsv_ch)
    }
}

/*
========================================================================================
    WORKFLOW COMPLETION
========================================================================================
*/

workflow.onComplete {
    log.info """\
        Pipeline completed at: ${workflow.complete}
        Execution status: ${workflow.success ? 'SUCCESS' : 'FAILED'}
        Duration: ${workflow.duration}
        Output directory: ${params.outdir}
        """
        .stripIndent()
}