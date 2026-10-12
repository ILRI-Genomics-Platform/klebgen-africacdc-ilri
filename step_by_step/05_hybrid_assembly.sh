#!/usr/bin/env bash
#
# Step 05 - hybrid assembly: long reads first, then polishing with short reads
#           (pipeline: ASSEMBLY -> AUTOCYCLER_ASSEMBLE, MEDAKA, POLYPOLISH_SHORTREAD, PYPOLCA_RUN)
#
#   bash step_by_step/05_hybrid_assembly.sh [sample ... | all]
#
# Only for samples with long reads (others are skipped). Needs: steps 01 and 03.
# Slow: Autocycler makes 4 read subsets x 2 assemblers = 8 assemblies (an hour or more per isolate).
#
#   5a Autocycler  : subsample the ONT reads, assemble each subset with Flye and Raven, and build one
#                    consensus from the 8 assemblies (circular chromosome and plasmids when they agree)
#   5b Medaka      : correct the consensus with the ONT reads (long-read polishing)
#   5c Polypolish  : correct remaining errors with the Illumina reads (repeat-aware short-read polishing)
#   5d Pypolca     : a second short-read polish (--careful: only changes well-supported positions)
#
# Output: 02_assembly/autocycler/<sample>.autocycler.fasta/.gfa,  02_assembly/hybrid/<sample>.hybrid.fasta
#         (steps 06 onwards use the hybrid assembly for this sample instead of Shovill's contigs)
# Look at: number of contigs and whether they are circular (Bandage on the .gfa); Pypolca's report of
#   how many errors it corrected (few = the assembly was already accurate).

source "$(dirname "$0")/config.sh"

aout="${OUTDIR}/02_assembly/autocycler"; hout="${OUTDIR}/02_assembly/hybrid"; mkdir -p "$aout" "$hout"
for s in $(samples_from_args "$@"); do
    if [ -z "$(field "$s" long_fastq)" ]; then info "$s: no long reads; skipped"; continue; fi
    long="$(filtered_ont "$s")"
    need_file "$long" "run step 03 (Filtlong) first"; need_file "$(trimmed_r1 "$s")" "run step 01 (fastp) first"
    work="${hout}/${s}_work"; mkdir -p "$work"

    step "5a Autocycler (Flye + Raven): $s"
    mkdir -p "${work}/autocycler"; cd "${work}/autocycler"
    run_sh AUTOCYCLER_ASSEMBLE "
        genome_size=\$(autocycler helper genome_size --reads $long --threads ${THREADS})
        echo \"estimated genome size: \$genome_size\"

        autocycler subsample --reads $long --out_dir subsampled_reads --genome_size \$genome_size

        mkdir -p assemblies
        for assembler in flye raven; do
            for read_set in subsampled_reads/sample_*.fastq; do
                i=\$(basename \$read_set .fastq | sed 's/sample_//')
                autocycler helper \$assembler --reads \$read_set --out_prefix assemblies/\${assembler}_\${i} \\
                    --threads ${THREADS} --genome_size \$genome_size
            done
        done

        autocycler compress -i assemblies -a autocycler_out
        autocycler cluster -a autocycler_out
        for cluster in autocycler_out/clustering/qc_pass/cluster_*; do
            autocycler trim -c \$cluster
            autocycler resolve -c \$cluster
        done
        autocycler combine -a autocycler_out -i autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa"
    cp autocycler_out/consensus_assembly.fasta "${aout}/${s}.autocycler.fasta"
    cp autocycler_out/consensus_assembly.gfa   "${aout}/${s}.autocycler.gfa"

    step "5b Medaka (long-read polishing): $s"
    mkdir -p "${work}/medaka"; cd "${work}/medaka"
    run MEDAKA medaka_consensus \
        -t "${THREADS}" \
        -i "$long" \
        -d "${aout}/${s}.autocycler.fasta" \
        -o ./
    cp consensus.fasta "${work}/${s}.medaka.fa"

    step "5c Polypolish (short-read polishing): $s"
    mkdir -p "${work}/polypolish"; cd "${work}/polypolish"
    cp "${work}/${s}.medaka.fa" draft.fasta
    run POLYPOLISH_SHORTREAD bwa index draft.fasta
    run_sh POLYPOLISH_SHORTREAD "bwa mem -t ${THREADS} -a draft.fasta $(trimmed_r1 "$s") > alignments_1.sam"
    run_sh POLYPOLISH_SHORTREAD "bwa mem -t ${THREADS} -a draft.fasta $(trimmed_r2 "$s") > alignments_2.sam"
    run POLYPOLISH_SHORTREAD polypolish filter \
        --in1 alignments_1.sam --in2 alignments_2.sam \
        --out1 filtered_1.sam  --out2 filtered_2.sam
    run_sh POLYPOLISH_SHORTREAD "polypolish polish draft.fasta filtered_1.sam filtered_2.sam > ${s}.polypolish.fasta"
    rm -f alignments_?.sam filtered_?.sam             # large; not needed any more

    step "5d Pypolca (short-read polishing): $s"
    cd "$work"
    run PYPOLCA_RUN pypolca run \
        -a "polypolish/${s}.polypolish.fasta" \
        -1 "$(trimmed_r1 "$s")" -2 "$(trimmed_r2 "$s")" \
        -t "${THREADS}" \
        -o "${s}.hybrid" \
        --prefix "${s}.hybrid" \
        --careful
    cp "${s}.hybrid/${s}.hybrid_corrected.fasta" "${hout}/${s}.hybrid.fasta"
    info "hybrid assembly: ${hout}/${s}.hybrid.fasta ($(grep -c '>' "${hout}/${s}.hybrid.fasta") contigs)"
done
