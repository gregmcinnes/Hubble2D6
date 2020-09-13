# Hubble.2D6

## Introduction

Hubble.2D6 is a tool to predict *CYP2D6* haplotype function.  Hubble.2D6 predicts function on an oridal scale: Normal Function, Decreased Function, and No Function, as used in CPIC guidelines\*.  You can use Hubble.2D6 to predict the function of *CYP2D6* star alleles that comprise any arbitrary set of variants within the *CYP2D6* locus.

Read more about Hubble.2D6 in the [manuscript](https://www.biorxiv.org/content/10.1101/684357v2.abstract).

\*Increased function is not predicted because currently the only known mechanism the leads to an increased function allele is through increased copy number, which is not evaluated by Hubble.2D6.  



## Installation

We recommend using a [conda](https://docs.conda.io/en/latest/) environment for the python libraries. Python version 3.6 required.

### Quick Installation

If you only need to analyze haplotypes containing SNVs and common INDELs (INDELs found in existing star alleles), you can use the quick installation.  We have precomputed embeddings for all possible SNVs and common INDELs so the lengthy annotation step is not needed. If you require analysis of INDELs not in existing star alleles, proceed to "Full Installation".

#### Intall Python libraries

```pip install -r requirements.txt```


### Full Installation
Hubble.2D6 uses annotation embeddings for input variants to predict star allele function.  A number of steps are required to generate the annotation embeddings, which requires installation of several tools.  If you only need to evaluate star alleles with SNVs, you can skip this step.  Annotation embeddings for SNVs have been precomputed.

_Fair warning, installation of VEP can be frustrating and the Annovar libraries are very large._

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

#### Annotations

**_If you are only processing SNVs and common INDELs you can skip this step._**

First the VCF needs to be converted to a sequence of variant embeddings
```
sh bin/vcf2seq.sh INPUT_VCF PREFIX
```

Example:
```
sh bin/vcf2seq.sh data/test.vcf test
```

This will generate a file with a .seq extension that can be used to predict function

#### VCF
When using SNVs and common INDELs, you should split multiallelic sites and realign INDELs. This ensures INDELs can be mapped consistently to annotation embeddings. These functions can be performed using [BCFtools](http://samtools.github.io/bcftools/bcftools.html).

Split multialleleic sites:
```
bcftools norm -m-any -o OUTPUT.vcf  INPUT.vcf
```

Normalize INDELs:
```
bcftools norm -f $REF_GENOME -c s -o OUTPUT.vcf INPUT.vcf
```

### Functional prediction

SNVs and common INDELs only can be mapped from VCF.  Generate star allele predictions using the following command:

```
python bin/hubble.py -s data/sample.vcf
```

If you require INDELs that have not been precomputed, after preparing the seq file as described above use this function to predict the haplotype function.

Example:
```
python bin/hubble.py -s data/test.seq
```






