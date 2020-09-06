#!/bin/bash
#SBATCH --job-name=vcf2seq
#SBATCH --output=vcf2seq.%j.out
#SBATCH --error=vcf2seq.%j.err
#SBATCH --time=48:00:00
#SBATCH -p rbaltman,owners,normal
#SBATCH --nodes=1
#SBATCH --mem 20G
#SBATCH --mail-type=FAIL

# Input: VCF
input_vcf=$1
prefix=$2

# Install path
HUBBLE_PATH=/oak/stanford/groups/rbaltman/gmcinnes/projects/cyp2d6/Hubble2D6
BCFTOOLS_PATH=/share/software/user/open/bcftools/1.8/bin
ANNOVAR_PATH=/oak/stanford/groups/rbaltman/gmcinnes/bin/annovar
VEP_PATH=/home/groups/rbaltman/gmcinnes/bin/vep/ensembl-vep
REF_GENOME=/oak/stanford/groups/rbaltman/gmcinnes/data/human_genome/hg19/hg19/hg19.fa

source ~/.bash_profile
conda activate hubble
module load perl

echo "Normalizing chromosome names"
awk '{if($0 !~ /^#/ && $0 !~ /^chr/) print "chr"$0; else print $0}' $input_vcf > $prefix.chr.vcf

# Normalize indels and flip sites that don't match thhe e asdf;lkj;
echo  "Normalizing INDELs"
$BCFTOOLS_PATH/bcftools norm -f $REF_GENOME -c s -o $prefix.norm.vcf $prefix.chr.vcf

# Merge multiallelic sites that may have been formed from the normalization
echo "Merging multiallelic sites"
python $HUBBLE_PATH/bin/merge_multiallelic.py --vcf $prefix.norm.vcf > $prefix.norm_m.vcf

# vcf2gvcf
echo  "Converting to gVCF"
python $HUBBLE_PATH/bin/vcf2gvcf.py --vcf $prefix.norm_m.vcf --fasta $REF_GENOME --bed $HUBBLE_PATH/data/cyp2d6.variant_region.bed --format GT > $prefix.gvcf

# Split all multi allelic sites and flip sites that don’t match the reference
echo "Splitting multi-allelic sites"
bcftools norm -m-any -o $prefix.split_snps.gvcf  $prefix.gvcf

# Annotate
## VEP
echo  "Running VEP"
/home/groups/rbaltman/gmcinnes/bin/vep/ensembl-vep/vep -i $prefix.split_snps.gvcf --merged --dir_cache /oak/stanford/groups/rbaltman/gmcinnes/bin/vep --species homo_sapiens --offline --fasta /oak/stanford/groups/rbaltman/gmcinnes/bin/vep/homo_sapiens/75_GRCh37/Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz  --force_overwrite --regulatory  --plugin LoF,loftee_path:/home/groups/rbaltman/gmcinnes/bin/loftee/loftee,human_ancestor_fa:/oak/stanford/groups/rbaltman/gmcinnes/bin/vep/loftee/human_ancestor.fa.gz,conservation_file:/oak/stanford/groups/rbaltman/gmcinnes/bin/vep/loftee/phylocsf_gerp.sql --dir_plugins /home/groups/rbaltman/gmcinnes/bin/loftee/loftee --allow_non_variant --vcf -o $prefix.vep.gvcf

## Annovar
echo "Running Annovar"
perl $ANNOVAR_PATH/annovar_latest/table_annovar.pl $prefix.vep.gvcf  $ANNOVAR_PATH/annovar_latest/humandb -buildver hg19 -out $prefix -remove -protocol wgEncodeAwgDnaseMasterSites,wgEncodeRegTfbsClusteredV3,tfbsConsSites,CYP2D6_eQTLs,wgEncodeHaibMethyl450Gm12878SitesRep1,wgEncodeHaibMethylRrbsGm12878HaibSitesRep1,cyp2d6_active_site,dann,cadd,dbnsfp33a,fathmm,gnomad_genome,refGene -operation r,r,r,f,r,r,r,f,f,f,f,f,g -nastring . -polish -vcfinput

## Custom
echo "Adding custom annotations"
python $HUBBLE_PATH/bin/add_custom_annotations.py --vcf $prefix.hg19_multianno.vcf --annotations rare,deleterious,coding,indel,active_site,tfbs,eqtl,dnase,methyl > $prefix.custom.gvcf

# Merge multiallelic
echo "Merging multiallelic sites"
python $HUBBLE_PATH/bin/merge_multiallelic.py --vcf $prefix.custom.gvcf > $prefix.final.gvcf

# vcf2seq
echo "Running vcf2seq"
python $HUBBLE_PATH/bin/vcf2seq.py --vcf $prefix.final.gvcf --annotations rare,deleterious,coding,indel,active_site,tfbs,eqtl,dnase,methyl --zero_pad_length 0 --no_indels > $prefix.seq

# Clean up intermediate files
echo "Cleaning up"
#rm -f $prefix.chr.vcf
#rm -f $prefix.norm.vcf
#rm -f $prefix.gvcf
#rm -f $prefix.split_snps.gvcf
#rm -f $prefix.vep.gvcf
#rm -f $prefix.hg19_multianno.vcf
#rm -f $prefix.custom.gvcf
#rm -f $prefix.hg19_multianno.txt
#rm -f $prefix.vep.vcf_summary.html
#rm -f $prefix.vep.vcf_warnings.txt
#rm -f $prefix.avinput
#rm -f $prefix.header.vcf
