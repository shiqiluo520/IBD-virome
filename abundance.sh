#!/bin/bash

prefix=$1
BAMS=$2
out_dir="out_abundance_${prefix}"
mkdir -p $out_dir
coverm contig --methods covered_fraction --bam-files $BAMS --min-read-percent-identity 0.95 1> ${out_dir}/abundance_contigs_covered_fraction.tsv 2> ${out_dir}/log_contig_covered_fraction.txt
sed -i '1 s/ Covered Fraction//g' ${out_dir}/abundance_contigs_covered_fraction.tsv
coverm contig --methods count --bam-files $BAMS --min-read-percent-identity 0.95 1> ${out_dir}/abundance_contigs_count.tsv 2> ${out_dir}/log_contig_count.txt
sed -i '1 s/ Read Count//g' ${out_dir}/abundance_contigs_count.tsv
coverm contig --methods tpm --bam-files $BAMS --min-read-percent-identity 0.95 1> ${out_dir}/abundance_contigs_tpm.tsv 2> ${out_dir}/log_contig_tpm.txt
sed -i '1 s/ TPM//g' ${out_dir}/abundance_contigs_tpm.tsv
coverm contig --methods rpkm --bam-files $BAMS --min-read-percent-identity 0.95 1> ${out_dir}/abundance_contigs_rpkm.tsv 2> ${out_dir}/log_contig_rpkm.txt
sed -i '1 s/ RPKM//g' ${out_dir}/abundance_contigs_rpkm.tsv