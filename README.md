# Hubble.2D6

## Introduction

Hubble.2D6 is a tool to predict *CYP2D6* haplotype function.  Hubble.2D6 predicts function on an oridal scale: Normal Function, Decreased Function, and No Function, as used in CPIC guidelines\*.  You can use Hubble.2D6 to predict the function of *CYP2D6* star alleles that comprise any arbitrary set of variants within the *CYP2D6* locus.

Read more about Hubble.2D6 in the [manuscript](https://www.biorxiv.org/content/10.1101/684357v2.abstract).

\*Increased function is not predicted because, currently, the only known mechanism the leads to an increased function allele is through increased copy number, which is not evaluated by Hubble.2D6.  



## Installation

### Quick Installation

--- This has not yet been implemented.  Coming very soon. Proceed to Full Installation. ---

This method only works to evaluate star alleles that comprise only SNVs.  All annotation embeddings have been precomputed for SNVs, so installation of annovar and VEP is not necessary.  If you require analysis of INDELs, proceed to "Full Installation"

#### Intall Python libraries

```pip install -r requirements.txt```


### Full Installation
Hubble.2D6 uses annotation embeddings for input variants to predict star allele function.  A number of steps are required to generate the annotation embeddings, which requires installation of several tools.  If you only need to evaluate star alleles with SNVs, you can skip this step.  Annotation embeddings for SNVs have been precomputed.  When running hubble use the flag `--no-indels` to skip the annotation step.

#### Intall Python libraries
* pip install -r requirements.txt

#### Intall Annovar
* Install [Annovar](https://annovar.openbioinformatics.org/en/latest/)
* Download required annovar databases.  Some of these are large databases that will take a while to download.  (ex. `annotate_variation.pl -downdb -buildver hg19 -webfrom annovar DATABASE_NAME humandb/`)
    * wgEncodeAwgDnaseMasterSites
    * wgEncodeRegTfbsClusteredV3
    * tfbsConsSites
    * wgEncodeHaibMethyl450Gm12878SitesRep1
    * wgEncodeHaibMethylRrbsGm12878HaibSitesRep1
    * dann
    * cadd
    * dbnsfp33a
    * fathmm
    * gnomad_genome
    * refGene
* Copy custom database files to annovar directory (ex. `cp data/DATABASE_NAME $ANNOVAR/humandb/`)
    * CYP2D6_eQTLs
    * cyp2d6_active_site

#### Intall VEP
* Install [VEP](https://uswest.ensembl.org/info/docs/tools/vep/index.html)
* Install VEP plugin [LOFTEE](https://github.com/konradjk/loftee)

#### Intall BCFtools
* Install [BCFtools](http://samtools.github.io/bcftools/bcftools.html)

#### Get the reference genome

```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz
```

Set global variables.  These do not need to be set for the quick installation.

```
HUBBLE_PATH=path_to_install_dir
ANNOVAR_PATH=path_to_annovar_dir
VEP_PATH=path_to_vep_dir
BCFtools=path_to_bcftools
REF_GENOME=path_to_ref_genomie
```

## Run

The input to Hubble.2D6 is a VCF and can have as many samples as needed.  To evaluate star alleles, all sites with an alternate allele
should be marked as having a homozygous alternate allele (1/1).  


### Data preprocessing
First the VCF needs to be converted to a sequence of variant embeddings
```
sh bin/vcf2seq.sh INPUT_VCF PREFIX
```

Example:
```
sh bin/vcf2seq.sh data/test.vcf test
```

This will generate a file with a .seq extension that can be used to predict function

### Functional prediction
Once you have the seq file predictions can be made for star allele function.  Samples and their 
functional predictions will be written to file specified with `-o`.
```
python bin/predict.py -f INPUT.seq -o OUT_FILE
```

Example:
```
python bin/predict.py -f data/test.seq -o test_predictions.txt
```






