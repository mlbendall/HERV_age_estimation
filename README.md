# HERV age estimation

To estimate the age of each locus, a sequence from the locus is aligned to a 
consensus sequence for the subfamily, and the genetic distance is calculated.

## Setup environment

Clone this repository

```bash
git clone git@github.com:mlbendall/HERV_age_estimation.git HERV_age_estimation.git && cd HERV_age_estimation.git
```

Create a conda environment with the necessary packages

```bash
conda create -n hervages python=2.7 requests future biopython intervaltree mafft
```

Activate the environment

```bash
conda activate hervages
```

Install `telebuilder` into this environment

```
git clone git@github.com:mlbendall/telebuilder.git && \
cd telebuilder && \
pip install . && \
cd ..
```

You must also have `R` installed in your current environment with the packages `tidyverse`, `ape`, and `cowplot`.

## 1. Create a database of consensus sequences

The repeat identifiers used in the Telescope annotation are the same as the names in the 
UCSC RepeatMasker database and follow the RepBase convention. However, we have decided to
use Dfam consensus sequences instead, so we need a way to map the Dfam sequences to the
RepBase names. The mapping is provided in [`rmsk_to_dfam.tsv`](./rmsk_to_dfam.tsv). The columns in this file are
RepBase identifier, internal or LTR, Dfam accession, and Dfam description.

> Note 1: There is no Dfam model corresponding to LTR7A in RepBase (We use the RepBase sequence for this subfamily)
> 
> Note 2: The description for [DF0000516](http://dfam.org/family/DF0000516/summary) has a typo. It says LTR40b in the description but the model is for LTR40c.

### 1A. Download consensus sequences from Dfam

Use the Dfam API to download consensus sequences. Query selects human (ancestors and descendants) with repeat type "LTR".

```
python dfam_consensus.py
```

##### Output: 

| | description |
| --- | --- |
| `ERV_human.dfam.fasta` | Consensus sequences from Dfam |
| `ERV_human.dfam.gb` | Consensus sequences and annotations from Dfam|


### 1B. Download consensus sequences from RepBase

Use a `curl` query to download consensus sequences. Query selects human (ancestors and descendants) with repeat type "Endogenous Retrovirus". Registration is required with [RepBase](https://www.girinst.org/accountservices/register.php) and the script will prompt for your username and password.

```
python repbase_consensus.py
```

##### Input:

| | description |
| --- | --- |
| `username` | RepBase user name |
| `password` | RepBase password |


##### Output: 

| | description |
| --- | --- |
| `ERV_human.repbase.fasta` | Consensus sequences from RepBase |
| `ERV_human.repbase.gb` | Consensus sequences and annotations from RepBase|

### 1C. Create the consensus sequence database

Creates a final consensus database that uses RepBase identifiers and Dfam consensus sequences.

```
python final_consensus.py
```

##### Input:

| | description |
| --- | --- |
| `ERV_human.dfam.gb` | Consensus sequences and annotations from Dfam|
| `ERV_human.repbase.fasta` | Consensus sequences from RepBase |
| [`rmsk_to_dfam.tsv`](./rmsk_to_dfam.tsv) | Mapping RepBase identifiers to Dfam accessions|

##### Output: 

| | description |
| --- | --- |
| `ERV_human.consensus.fasta` | Consensus database |

## 2. `age_pipeline.sh`

The bash script [`age_pipeline.sh`](./age_pipeline.sh) encapsulates all the steps needed to perform ages estimates
for an ERV subfamily. The main steps are to 2A) build the ERV annotation, 
2B) extract the sequences for each locus, and 2C) estimate ages for each locus.
Below each step is described in detail.

The following loop will estimate distances for all families in [`families.tsv`](./families.tsv):

```bash
IFS=$'\t' grep -v '^#' families.tsv | while read n im lm ; do
    echo -e "\n\nAge pipeline $n"
    /bin/bash age_pipeline.sh $n $im $lm
done
```

### 2A. Build ERV loci

A GTF file containing the locations of ERV loci can be obtained using the `buildERV` 
method in `telebuilder`. This GTF file contains separate annotations for each "hit", and 
these annotations are grouped into loci with the "locus" tag. 

**NOTE:** If estimating ages of an 
existing annotation, it is best to skip this step and instead use the GTFs that were 
generated while building the annotation. 

Example command:

```bash
buildERV --auto --no_igv ${name} ${intmodel} ${ltrmodel}
```

where `$name` is the name for the subfamily, `$intmodel` is the RepBase identifier
for the internal subfamily, and `$ltrmodel` is a comma-separated list of RepBase 
identifiers for LTR subfamilies. A table containing information about each locus is 
generated using `gtftools`:

```bash
gtftools tsv < ${name}/${name}.gtf > ${name}/${name}.tsv
```

##### Input:

| | description |
| --- | --- |
| [`build.txt`](./build.txt) | genome build, i.e "hg38" |
| [`chrom.sizes`](./chrom.sizes) | chromosome names and sizes in desired sort order |
| [`cytoband.gtf`](./cytoband.gtf) | GTF containing cytogenetic band names and ranges |
| `$name` | subfamily name |
| `$intmodel` | RepBase id for internal region |
| `$ltrmodel` | comma-separated list of RepBase ids for LTRs |

##### Output: 

| | description |
| --- | --- |
| `$name/$name.gtf` | Locations of ERV loci and constituent hits |
| `$name/$name.tsv` | Table with data about ERV loci |


### 2B. Extract ERV sequences

The nucleotide sequence for each locus is downloaded from UCSC and saved as a FASTA file,
one record per locus. Sequences are oriented in the sense direction (negative stranded loci are
reverse complemented). At the same time, the original annotations are adjusted for the extracted
sequence and a GTF annotating the extracted regions is created. Example command:

```bash
gtftools extract --gtfout ${name}/${name}_extracted.gtf ${name}/${name}.gtf > ${name}/${name}.fna
```

Geneious provides a nice way to visualize the extracted loci, but a small change needs to be
made to the GTF file using the `geneious_gtf` script: 


```bash
./geneious_gtf < ${name}/${name}_extracted.gtf > ${name}/${name}_extracted.geneious.gtf
```

To load in Geneious, select both `$name/$name.fna` and `$name/${name}_extracted.geneious.gtf`
and drag them into the window at the same time.

##### Input:

| | description |
| --- | --- |
| [`build.txt`](./build.txt) | genome build, i.e "hg38" |
| `$name/$name.gtf` | Locations of ERV loci and constituent hits |


##### Output: 

| | description |
| --- | --- |
| `$name/$name.fna` | FASTA with sequence for each locus |
| `$name/${name}_extracted.gtf` | GTF with hits relative to extracted FASTA |
| `$name/${name}_extracted.geneious.gtf` | As above, but suitable for loading in Geneious |

### 2C. Estimate ages

The age of each locus is calculated by selecting the "best" LTR that represents the locus,
aligning the LTR to the LTR consensus sequence, and calculating the genetic distance using
several measures of evolutionary distance. The "best" LTR is selected by identifying
flanking 5' and 3' LTRs; the longer of the two is chosen for loci that have both LTRs present.
The locus LTR and consensus LTR are aligned using `MAFFT`. The evolutionary distance
between the aligned sequences is calculated using the `ape` package in R. 

The distance metrics implemented are:

+ **p.dist** (raw, N): This is simply the proportion or the number of sites that differ between each pair of sequences. 
+ **jc.dist** (JC69): This model was developed by Jukes and Cantor (1969). It assumes that all substitutions (i.e. a change of a base by another one) have the same probability. This probability is the same for all sites along the DNA sequence. Another assumption is that the base frequencies are balanced and thus equal to 0.25.
+ **tn.dist** (TN93): Tamura and Nei (1993) developed a model which assumes distinct rates for both kinds of transition (A <-> G versus C <-> T), and transversions. The base frequencies are not assumed to be equal and are estimated from the data.
+ **k2p.dist** (K80): The distance derived by Kimura (1980), sometimes referred to as “Kimura's 2-parameters distance”, has the same underlying assumptions than the Jukes–Cantor distance except that two kinds of substitutions are considered: transitions (A <-> G, C <-> T), and transversions (A <-> C, A <-> T, C <-> G, G <-> T). They are assumed to have different probabilities. A transition is the substitution of a purine (C, T) by another one, or the substitution of a pyrimidine (A, G) by another one. A transversion is the substitution of a purine by a pyrimidine, or vice-versa. Both transition and transversion rates are the same for all sites along the DNA sequence.
+ **t3p.dist** (T92): Tamura (1992) generalized the Kimura model by relaxing the assumption of equal base frequencies. This is done by taking into account the bias in G+C content in the sequences. The substitution rates are assumed to be the same for all sites along the DNA sequence.
+ **cof.dist**: Like p-distance, but initial gaps are counted as a mismatch. No penalty for gap continuations. Described by Subramanian et al. (2011).


```bash
python estimate_ages.py \
    ${name}/${name}_extracted.gtf \
    ${name}/${name}.fna \
    ERV_human.consensus.fasta > ${name}/distances.tsv
```

Finally, we create a plot showing the age distribution by LTR subfamily using the
Kimura distance. 

```bash
Rscript plotages.R ${name}/distances.tsv ${name}/age_distribution.pdf
```

##### Input:

| | description |
| --- | --- |
| `$name/${name}_extracted.gtf` | GTF with hits relative to extracted FASTA |
| `$name/$name.fna` | FASTA with sequence for each locus |
| ERV_human.consensus.fasta | Database of consensus sequences |

##### Output: 

| | description |
| --- | --- |
| `$name/distances.tsv` | Table with distances for each locus |
| `$name/age_distribution.pdf` | Plot with age distribution by LTR model |


## 3. Analysis

This parses all the `distance.tsv` files and plots the age distribution by subfamily, ordered from youngest to oldest. Internal loci and subfamilies with 2 or fewer loci are removed.

```bash
Rscript plotages_all.R
```
