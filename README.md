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
conda create -n hervages python=2.7 future biopython intervaltree mafft
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

## 1. Create a database of consensus sequences

The repeat identifiers used in the Telescope annotation are the same as the names in the 
UCSC RepeatMasker database and follow the RepBase convention. However, we have decided to
use Dfam consensus sequences instead, so we need a way to map the Dfam sequences to the
RepBase names. The mapping is provided in [`rmsk_to_dfam.tsv`](./rmsk_to_dfam.tsv). The columns in this file are
RepBase identifier, internal or LTR, Dfam accession, and Dfam description.

> Note 1: There is no Dfam model corresponding to LTR7A in RepBase (We use the RepBase sequence for this subfamily)
> 
> Note 2: The description for [DF0000516](http://dfam.org/family/DF0000516/summary) has a typo. It says LTR40b in the description but the model is for LTR40c.

#### 1A. Download consensus sequences from Dfam

Use the Dfam API to download consensus sequences. Query selects human (ancestors and descendants) with repeat type "LTR".

```
python dfam_consensus.py
```

Creates two files, `ERV_human.dfam.gb` and `ERV_human.dfam.fasta`.

#### 1B. Download consensus sequences from RepBase

Use a `curl` query to download consensus sequences. Query selects human (ancestors and descendants) with repeat type "Endogenous Retrovirus". Registration is required with [RepBase](https://www.girinst.org/accountservices/register.php) and the script will prompt for your username and password.

```
python repbase_consensus.py
```

Creates two files, `ERV_human.dfam.gb` and `ERV_human.dfam.fasta`.

#### 1C. Create the consensus sequence database

Creates a final consensus database that uses RepBase identifiers and Dfam consensus sequences.

```
python final_consensus.py
```

Creates `ERV_human.consensus.fasta`.

## 2. 