#!/usr/bin/env nextflow
/*
========================================================================================
    neonatal-amr-nf
========================================================================================
    Training pipeline for neonatal AMR genomic surveillance of the Klebsiella
    pneumoniae species complex and Escherichia coli (Illumina + Oxford Nanopore).

    Layout follows nf-core conventions:
      modules/nf-core/   official nf-core modules (see modules.json for versions)
      modules/local/     pipeline-specific modules (main.nf, meta.yml, environment.yml)
      subworkflows/local one subworkflow per course module group
      workflows/         the main analysis workflow
      conf/              resources (base.config), tool arguments (modules.config),
                         site configs (slurm.config, local.config)

    Help:  nextflow run main.nf --help
----------------------------------------------------------------------------------------
*/

include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_neonatal_amr_pipeline/main'
include { NEONATAL_AMR            } from './workflows/neonatal_amr'

workflow {
    main:
    if (params.help) {
        // --help: print usage and stop; no other option is needed
        log.info(helpMessage())
    }
    else if (!params.input) {
        // Missing --input: show how to run the pipeline, then stop with a short error
        log.info(helpMessage())
        error("Missing required option --input (a samplesheet CSV). See the usage above.")
    }
    else {
        // Shared read-only images: check the course-built images exist before starting
        if (workflow.containerEngine == 'singularity' && params.local_images) {
            def missing = params.local_images
                .collect { process, image -> "${params.shared_container_dir}/local/${image}" }
                .findAll { path -> !file(path).exists() }
            if (missing) {
                error(
                    "Container image(s) missing from the shared read-only image folder: ${missing.join(', ')}. " +
                    "Participants cannot add images: ask a trainer to run " +
                    "'bin/prepare_shared_containers.sh --dir ${params.shared_container_dir} --check'."
                )
            }
        }

        PIPELINE_INITIALISATION(params.input)
        NEONATAL_AMR(PIPELINE_INITIALISATION.out.samplesheet)

        // With the shared read-only images (-profile singularity -c conf/slurm.config or
        // conf/local.config) a missing image shows up as "Permission denied" on a '.lock'
        // file. Tell the participant what it means.
        def run = workflow
        def imageDir = params.shared_container_dir
        workflow.onError {
            def message = (run.errorMessage ?: '') + ' ' + (run.errorReport ?: '')
            if (message.contains('.lock') && message.contains('Permission denied')) {
                log.error(
                    "A container image is missing from the shared read-only image folder (${imageDir}). " +
                    "Participants cannot add images: ask a trainer to run " +
                    "'bin/prepare_shared_containers.sh --dir ${imageDir} --check' and add the missing image."
                )
            }
            else if (message.contains('Missing Singularity library directory')) {
                log.error(
                    "The shared image folder ${imageDir}/nf-core does not exist. Check shared_container_dir " +
                    "in your site config (-c conf/slurm.config or conf/local.config), or pass --shared_container_dir."
                )
            }
        }
    }
}

//
// Usage message (nextflow run main.nf --help)
//
def helpMessage() {
    return """
    ======================================================================================
     neonatal-amr-nf  ${workflow.manifest.version ?: ''}
     Neonatal AMR genomic surveillance: K. pneumoniae species complex and E. coli
    ======================================================================================

    Usage:
      nextflow run main.nf --input samplesheet.csv --outdir results [options]

    Typical runs:
      # ILRI HPC (Slurm) with the course's shared, read-only Singularity images
      nextflow run main.nf -profile singularity -c conf/slurm.config \\
          --input samplesheet.csv --outdir results --kraken2_db /path/to/kraken2_db

      # A single training server with the same images
      nextflow run main.nf -profile singularity -c conf/local.config \\
          --input samplesheet.csv --outdir results

      # Check the wiring without running any tools
      nextflow run main.nf -profile test,singularity -c conf/local.config -stub --outdir test_results

    Required:
      --input                  Samplesheet CSV, one of:
                               reads:    sample,species,fastq_1,fastq_2,long_fastq,phylo_group
                               assembly: sample,species,fasta,phylo_group
                                 species     kpneumoniae | ecoli (or 'K. pneumoniae', 'E. coli')
                                 long_fastq  optional Nanopore reads (hybrid assembly)
                                 fasta       contigs; read QC, Kraken2 and assembly are skipped
                                 phylo_group optional label (e.g. ST39); used with
                                             --phylo_grouping samplesheet
      --input_type             auto (from the columns), reads or assembly [${params.input_type}]
                               Build an assembly samplesheet from metadata and a contigs folder:
                               bin/make_assembly_samplesheet.py --help

    Output:
      --outdir                 Results folder                          [${params.outdir}]

    Databases:
      --kraken2_db             Kraken2 database folder (with Bracken files);
                               omit to skip the contamination check    [${params.kraken2_db ?: 'not set'}]
      --amrfinderplus_db       AMRFinderPlus database (.tar.gz or folder);
                               omit to download the latest             [${params.amrfinderplus_db ?: 'not set'}]

    Phylogeny / transmission:
      --skip_phylogeny         Skip the phylogeny step                 [${params.skip_phylogeny}]
      --phylo_grouping         kleborate: group by Kleborate species + ST, reference
                                          per species/ST from --phylo_reference_sheet
                               samplesheet: group by the phylo_group column
                                                                       [${params.phylo_grouping}]
      --phylo_reference_sheet  CSV species,st,name,accession            [${params.phylo_reference_sheet}]
      --phylo_reference_dir    Folder with <accession>.fasta (bin/fetch_references.py)
                                                                       [${params.phylo_reference_dir}]
      --phylo_references       samplesheet mode: CSV 'group,fasta';
                               default: each group's best own assembly [${params.phylo_references ?: 'not set'}]
      --phylo_min_group_size   Smallest group analysed (tree: 3+)      [${params.phylo_min_group_size}]

    Typing and annotation:
      --skip_kaptive           No standalone Kaptive report (Kleborate
                               still reports K/O loci)                 [${params.skip_kaptive}]
      --run_prokka             Annotate every assembly with Prokka     [${params.run_prokka}]

    Containers:
      --shared_container_dir   Read-only folder of prepared Singularity images
                               (set in conf/slurm.config or conf/local.config)
                                                                       [${params.shared_container_dir}]

    Other:
      --multiqc_config         Custom MultiQC config file              [${params.multiqc_config ?: 'not set'}]
      --help                   Show this message and exit

    Profiles (-profile):
      singularity              Run tools in Singularity/Apptainer images; add a site config with -c
      docker, apptainer        Other container engines
      conda, mamba             Install tools with conda instead of containers
      wave                     Build containers for the local modules with Wave
      laptop                   Cap resources for an 8-core / 32 GB machine
      slurm                    Submit jobs with Slurm (without a site config)
      test                     Small built-in test data (use with -stub)

    Site configs (-c):
      conf/slurm.config        ILRI HPC: Slurm queues and nodes, shared images
      conf/local.config        One shared training server, shared images
    """.stripIndent()
}
