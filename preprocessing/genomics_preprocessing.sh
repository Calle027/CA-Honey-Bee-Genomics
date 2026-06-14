#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Preprocessing workflow for California honey bee genomics
#
# This script documents the command-line workflow used to produce:
#   final_cleaned.vcf.gz
#   final_cleaned.gds
#   genome_coverage.tsv
#
# Edit the paths and sample lists below before running on a new system.
###############################################################################

PROJECT_DIR="${PROJECT_DIR:-/path/to/project}"
REFERENCE="${REFERENCE:-/path/to/Amel_HAv3.1.fa}"
GATK="${GATK:-/path/to/gatk}"
TRIMMOMATIC_JAR="${TRIMMOMATIC_JAR:-/path/to/trimmomatic-0.39.jar}"
ADAPTERS="${ADAPTERS:-/path/to/TruSeq3-PE-2.fa}"
THREADS="${THREADS:-16}"

SAMPLES_FILE="${SAMPLES_FILE:-samples.txt}"
SRA_MAP="${SRA_MAP:-sra_samples.tsv}"
INTERVALS="${INTERVALS:-Amel_HAv3.1_interval.list}"

RAW_DIR="${PROJECT_DIR}/01_raw_fastq"
TRIM_DIR="${PROJECT_DIR}/02_trimmed_fastq"
MAP_DIR="${PROJECT_DIR}/03_mapped"
BAM_DIR="${PROJECT_DIR}/04_bam"
GVCF_DIR="${PROJECT_DIR}/05_gvcf"
VCF_DIR="${PROJECT_DIR}/06_vcf"
LOG_DIR="${PROJECT_DIR}/logs"

mkdir -p "${RAW_DIR}" "${TRIM_DIR}" "${MAP_DIR}" "${BAM_DIR}" "${GVCF_DIR}" "${VCF_DIR}" "${LOG_DIR}"

fetch_sra() {
  local srr="$1"
  local sample="$2"

  fastq-dump "${srr}" --split-3 --skip-technical --outdir "${RAW_DIR}"

  [[ -f "${RAW_DIR}/${srr}_1.fastq" ]] && pigz "${RAW_DIR}/${srr}_1.fastq"
  [[ -f "${RAW_DIR}/${srr}_2.fastq" ]] && pigz "${RAW_DIR}/${srr}_2.fastq"
  [[ -f "${RAW_DIR}/${srr}.fastq" ]] && pigz "${RAW_DIR}/${srr}.fastq"

  if [[ -f "${RAW_DIR}/${srr}_1.fastq.gz" && -f "${RAW_DIR}/${srr}_2.fastq.gz" ]]; then
    mv "${RAW_DIR}/${srr}_1.fastq.gz" "${RAW_DIR}/${sample}_R1.fastq.gz"
    mv "${RAW_DIR}/${srr}_2.fastq.gz" "${RAW_DIR}/${sample}_R2.fastq.gz"
  elif [[ -f "${RAW_DIR}/${srr}.fastq.gz" ]]; then
    mv "${RAW_DIR}/${srr}.fastq.gz" "${RAW_DIR}/${sample}_SE.fastq.gz"
  fi
}

trim_reads() {
  local sample="$1"

  java -jar "${TRIMMOMATIC_JAR}" PE -threads "${THREADS}" \
    "${RAW_DIR}/${sample}_R1.fastq.gz" \
    "${RAW_DIR}/${sample}_R2.fastq.gz" \
    "${TRIM_DIR}/${sample}_R1_trimmed.fastq.gz" \
    "${TRIM_DIR}/${sample}_R1_unpaired.fastq.gz" \
    "${TRIM_DIR}/${sample}_R2_trimmed.fastq.gz" \
    "${TRIM_DIR}/${sample}_R2_unpaired.fastq.gz" \
    "ILLUMINACLIP:${ADAPTERS}:2:30:10" \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
}

map_reads() {
  local sample="$1"

  bwa mem -t "${THREADS}" "${REFERENCE}" \
    "${TRIM_DIR}/${sample}_R1_trimmed.fastq.gz" \
    "${TRIM_DIR}/${sample}_R2_trimmed.fastq.gz" \
    > "${MAP_DIR}/${sample}.sam"
}

make_bam() {
  local sample="$1"

  samtools view -@ "${THREADS}" -q 20 -bS "${MAP_DIR}/${sample}.sam" \
    > "${BAM_DIR}/${sample}.filtered.bam"

  samtools sort -@ "${THREADS}" \
    -o "${BAM_DIR}/${sample}.sort.bam" \
    "${BAM_DIR}/${sample}.filtered.bam"

  samtools index "${BAM_DIR}/${sample}.sort.bam"
}

mark_duplicates() {
  local sample="$1"

  "${GATK}" MarkDuplicates \
    --REMOVE_DUPLICATES true \
    --METRICS_FILE "${LOG_DIR}/${sample}.dedup.metrics.txt" \
    --INPUT "${BAM_DIR}/${sample}.sort.bam" \
    --OUTPUT "${BAM_DIR}/${sample}.dedup.bam"

  samtools index "${BAM_DIR}/${sample}.dedup.bam"
}

add_read_groups() {
  local sample="$1"

  "${GATK}" AddOrReplaceReadGroups \
    --INPUT "${BAM_DIR}/${sample}.dedup.bam" \
    --OUTPUT "${BAM_DIR}/${sample}.RG.bam" \
    --RGLB lib1 \
    --RGPL illumina \
    --RGPU unit1 \
    --RGSM "${sample}" \
    --CREATE_INDEX true
}

calculate_coverage() {
  echo -e "sample\tbam\tavg_coverage" > "${PROJECT_DIR}/genome_coverage.tsv"

  while read -r sample; do
    local bam="${BAM_DIR}/${sample}.RG.bam"
    local avg_coverage
    avg_coverage=$(samtools depth -a "${bam}" | awk '{sum += $3} END {if (NR > 0) print sum / NR; else print "NA"}')
    echo -e "${sample}\t$(basename "${bam}")\t${avg_coverage}" >> "${PROJECT_DIR}/genome_coverage.tsv"
  done < "${SAMPLES_FILE}"
}

call_gvcf() {
  local sample="$1"

  "${GATK}" HaplotypeCaller \
    -R "${REFERENCE}" \
    -I "${BAM_DIR}/${sample}.RG.bam" \
    -O "${GVCF_DIR}/${sample}.g.vcf.gz" \
    -ERC GVCF \
    -G StandardAnnotation \
    -G AS_StandardAnnotation
}

make_sample_map() {
  while read -r sample; do
    echo -e "${sample}\t${GVCF_DIR}/${sample}.g.vcf.gz"
  done < "${SAMPLES_FILE}" > "${GVCF_DIR}/sample.cohort.map"
}

import_genomicsdb() {
  "${GATK}" GenomicsDBImport \
    --genomicsdb-workspace-path "${GVCF_DIR}/genomicsdb" \
    --batch-size 100 \
    --sample-name-map "${GVCF_DIR}/sample.cohort.map" \
    --reader-threads 12 \
    --max-num-intervals-to-import-in-parallel 17 \
    --intervals "${INTERVALS}"
}

joint_genotype() {
  "${GATK}" GenotypeGVCFs \
    -R "${REFERENCE}" \
    -V "gendb://${GVCF_DIR}/genomicsdb" \
    -G StandardAnnotation \
    -G AS_StandardAnnotation \
    --annotation FisherStrand \
    --annotation MappingQualityRankSumTest \
    --annotation QualByDepth \
    --annotation ReadPosRankSumTest \
    --annotation MappingQuality \
    -O "${VCF_DIR}/joint.raw.vcf.gz"

  bcftools view -m2 -M2 -v snps -Oz --threads "${THREADS}" \
    -o "${VCF_DIR}/joint.biallelic.snps.vcf.gz" \
    "${VCF_DIR}/joint.raw.vcf.gz"

  bcftools index "${VCF_DIR}/joint.biallelic.snps.vcf.gz"
}

filter_variants() {
  "${GATK}" IndexFeatureFile \
    -F "${VCF_DIR}/joint.biallelic.snps.vcf.gz"

  "${GATK}" VariantAnnotator \
    -R "${REFERENCE}" \
    -V "${VCF_DIR}/joint.biallelic.snps.vcf.gz" \
    -O "${VCF_DIR}/annotated.snps.vcf.gz" \
    --annotation Coverage \
    --annotation QualByDepth \
    --annotation FisherStrand \
    --annotation RMSMappingQuality \
    --annotation MappingQualityRankSumTest \
    --annotation ReadPosRankSumTest

  "${GATK}" VariantFiltration \
    -R "${REFERENCE}" \
    -V "${VCF_DIR}/annotated.snps.vcf.gz" \
    -O "${VCF_DIR}/filt.final.snps.vcf.gz" \
    --filter-name "Low_QD" -filter "QD < 5.0" \
    --filter-name "High_FS" -filter "FS > 11.0" \
    --filter-name "Low_MQ" -filter "MQ < 40.0" \
    --filter-name "Low_MQRankSum" -filter "MQRankSum < -2.0" \
    --filter-name "High_MQRankSum" -filter "MQRankSum > 2.0" \
    --filter-name "Low_ReadPosRankSum" -filter "ReadPosRankSum < -2.0" \
    --filter-name "High_ReadPosRankSum" -filter "ReadPosRankSum > 2.0" \
    --filter-name "Low_DP" -filter "DP < 20"

  bcftools view -h "${VCF_DIR}/filt.final.snps.vcf.gz" > "${VCF_DIR}/tmp.header"
  bcftools view -H "${VCF_DIR}/filt.final.snps.vcf.gz" \
    | awk -F '\t' '($7 == "Low_DP" || index($7, ";") > 0) { next } 1' \
    > "${VCF_DIR}/tmp.variants"

  cat "${VCF_DIR}/tmp.header" "${VCF_DIR}/tmp.variants" \
    | bgzip > "${PROJECT_DIR}/final_cleaned.vcf.gz"

  bcftools index -t "${PROJECT_DIR}/final_cleaned.vcf.gz"
  rm "${VCF_DIR}/tmp.header" "${VCF_DIR}/tmp.variants"
}

make_gds() {
  Rscript -e 'args <- commandArgs(TRUE); SeqArray::seqVCF2GDS(args[1], args[2], storage.option = "ZIP_RA")' \
    "${PROJECT_DIR}/final_cleaned.vcf.gz" \
    "${PROJECT_DIR}/final_cleaned.gds"
}

run_fastq_to_bam() {
  cut -f1 "${SRA_MAP}" | tail -n +2 > "${LOG_DIR}/srr.ids"
  cut -f2 "${SRA_MAP}" | tail -n +2 > "${LOG_DIR}/sample.ids"

  parallel --link --joblog "${LOG_DIR}/fetch_sra.log" fetch_sra :::: "${LOG_DIR}/srr.ids" :::: "${LOG_DIR}/sample.ids"
  parallel --joblog "${LOG_DIR}/trim_reads.log" trim_reads :::: "${SAMPLES_FILE}"
  parallel --joblog "${LOG_DIR}/map_reads.log" map_reads :::: "${SAMPLES_FILE}"
  parallel --joblog "${LOG_DIR}/make_bam.log" make_bam :::: "${SAMPLES_FILE}"
  parallel --joblog "${LOG_DIR}/mark_duplicates.log" mark_duplicates :::: "${SAMPLES_FILE}"
  parallel --joblog "${LOG_DIR}/add_read_groups.log" add_read_groups :::: "${SAMPLES_FILE}"
}

run_variant_calling() {
  parallel --joblog "${LOG_DIR}/call_gvcf.log" call_gvcf :::: "${SAMPLES_FILE}"
  make_sample_map
  import_genomicsdb
  joint_genotype
  filter_variants
  make_gds
}

export -f fetch_sra trim_reads map_reads make_bam mark_duplicates add_read_groups call_gvcf
export PROJECT_DIR REFERENCE GATK TRIMMOMATIC_JAR ADAPTERS THREADS
export RAW_DIR TRIM_DIR MAP_DIR BAM_DIR GVCF_DIR VCF_DIR LOG_DIR

case "${1:-}" in
  fastq_to_bam)
    run_fastq_to_bam
    ;;
  coverage)
    calculate_coverage
    ;;
  variant_calling)
    run_variant_calling
    ;;
  all)
    run_fastq_to_bam
    calculate_coverage
    run_variant_calling
    ;;
  *)
    echo "Usage: $0 {fastq_to_bam|coverage|variant_calling|all}"
    exit 1
    ;;
esac
