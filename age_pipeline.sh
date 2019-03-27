#! /bin/bash

name=$1
intmodel=$2
ltrmodel=$3

echo "Name:           ${name}"
echo "Internal model: ${intmodel}"
echo "LTR model(s):   ${ltrmodel}"

# Make LTRs array
IFS=',' read -ra ltrA <<< "$ltrmodel"

mkdir -p ${name}

echo -e "[--- Building ERV loci ---]"
if [[ ! -e ${name}/${name}.gtf ]]; then
    buildERV --auto --no_igv ${name} ${intmodel} ${ltrmodel} 2>&1 | tee ${name}/build.log
fi

# Create TSV
gtftools tsv < ${name}/${name}.gtf > ${name}/${name}.tsv

echo -e "[--- Extracting ERV sequences ---]"
gtftools extract --gtfout ${name}/${name}_extracted.gtf ${name}/${name}.gtf > ${name}/${name}.fna
./geneious_gtf < ${name}/${name}_extracted.gtf > ${name}/${name}_extracted.geneious.gtf

# Estimate ages
echo -e "\n\n[--- Estimating ages ---]\n\n"
python estimate_ages.py ${name}/${name}_extracted.gtf ${name}/${name}.fna ERV_human.consensus.fasta > ${name}/distances.tsv

# Plot
Rscript plotages.R ${name}/distances.tsv ${name}/age_distribution.pdf
