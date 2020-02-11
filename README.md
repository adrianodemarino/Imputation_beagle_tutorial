# Imputation beagle tutorial
Throughout the protocol we assume Bash shell.

This is a tutorial forked from: https://www.protocols.io/run/genotype-imputation-workflow-v3-0-xbgfijw

Split tutorial step by step.

1. Installation anaconda
2. Installation software
3. 

## Create new env that we call imputation

`conda create -n imputation python=3.6 anaconda`

Activate your new env:

`conda activate imputation`

Installation of the required pack/software:
```
conda install -c bioconda eagle
conda install -c bioconda beagle
conda install -c r r-base
conda install -c bioconda bcftools
conda install -c conda-forge r-data.table
conda install -c conda-forge r-sm
```

## or Download and install the software packages

Required software packages are listed below with the versions used in this protocol. However, using the latest versions is
recommended.
* BCFtools v1.7 (or later version) http://www.htslib.org/download/
* R v3.4.1 (or later version) https://www.r-project.org/
* R package data.table https://github.com/Rdatatable/data.table/wiki/Installation
* R package sm https://cran.r-project.org/package=sm
* Eagle v2.3.5 https://data.broadinstitute.org/alkesgroup/Eagle/
* Beagle v4.1 beagle.27Jan18.7e1.jar https://faculty.washington.edu/browning/beagle/b4_1.html
* Beagle bref bref.27Jan18.7e1.jar https://faculty.washington.edu/browning/beagle/bref.27Jan18.7e1.jar

## Reference genome and genetic map files
 
Fasta files
Homo Sapiens assembly hg38 version 0 is used and the required files are:
Homo_sapiens_assembly38.fasta
Homo_sapiens_assembly38.fasta.fai
The files are available for downloading at Broad Insitute storage in Google cloud at: [hg38](https://console.cloud.google.com/storage/browser/broad-references/hg38/v0/?pli=1)


1.2.2. Genetic map files for phasing with Eagle
Genetic map file (all chromosomes in a single file) with recombination frequencies for GRCh38/hg38 are available for downloading at Eagle download page at:
https://data.broadinstitute.org/alkesgroup/Eagle/downloads/tables/
genetic_map_hg38_withX.txt.gz

We have processed the file according to the command below in order to split it per chromosome with correct headers.
The resulting files are saved as:
eagle_chr#_b38.txt (where # is the chromosome number)

```
for CHR in {1..23}; do
    zcat genetic_map_hg38_withX.txt.gz | \
    grep ^${CHR} | \
    sed '1ichr position COMBINED_rate(cM/Mb) Genetic_Map(cM)' \
    > eagle_chr${CHR}_b38.map
done
```
Note: Currently the chromosome notation in the Eagle genetic map files is only the chromosome number without 'chr' and chrX is '23'. Starting from Eagle v2.4, also chromosome notation with 'chr' tag is supported.


1.2.3 Genetic map files for imputation with Beagle
Genetic map files fro Beagle are available for downloading at Beagle download page at
http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/
plink.GRCh38.map.zip

Unzip the files and change the chromosome notation from PLINK format (only number or X) to GRCh38/hg38 standard notation with 'chr' tag as follows:

```
# Unzip the files
unzip plink.GRCh38.map.zip

# Rename chromosome 23
mv plink.chrX.GRCh38.map plink.chr23.GRCh38.map

# Add 'chr' tag to the beginning of the line and 
# store the output with suitable filename
for CHR in {1..23}; do 
    cat plink.chr${CHR}.GRCh38.map | \
    sed 's/^/chr/' \
    > beagle_chr${CHR}_b38.map
done
```

Note: For GRCh38/hg38, the chromosome notation in the Beagle genetic map files is 'chr#' and chromosome 23 is 'chrX'.
1.3. Imputation reference panel files
 
1.3.1 Obtain the reference panel files
For increased imputation accuracy, we recommend using a population-specific imputation reference panel, if available.
 
If population-specific reference data is not available, for instance 1000 Genomes Project (1000 GP) (www.nature.com/articles/nature15393) data can be used instead.
 
GRCh38/hg38 files are available at EBI 1000 genomes ftp site: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/

```
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr{{1..22},X}_GRCh38.genotypes.20170504.vcf.gz{,.tbi}
```
NOTE: 
The reference panel files should contain: 
* phased genotypes,
* chromosome names with 'chr' and chromosome X as 'chrX', 
* all variants as biallelic records,
* only SNPs and INDELs,
* only unique variants, 
* non-missing data, 
* chrX as diploid genotypes, and
* only unique IDs

1.3.2 Minimum quality control

Here, we have piped most of the processing steps together in order to save significant amount of time by avoiding writing out multiple intermediate files. If your imputation reference panel does not require all the steps, modify the command accordingly.

For 1000GP data:
Generate a text file containing space-spearated old and new chromosome names. This is required to rename the numerical chromosome names with 'chr' tag. Apply the new chromosome names with 'bcftools annotate'.
Remove the rare variants, here singletons and doubletons by setting AC threshold with 'bcftools view'.
Split multiallelic sites to biallelic records with 'bcftools norm'.
Keep only SNPs and INDELs with 'bcftools view'. Here, the 1000GP data included a tag VT in the INFO field and data contain also structural variants which should be excluded.
Align the variants to reference genome with 'bcftools norm' in order to have the REF and ALT alleles in the shortest possible representation and to confirm that the REF allele matches the reference genome, additionally remove duplicate variants (-d none).
After alignment, remove multiallelic records with 'bcftools view', since these are formed during the alignment if the REF does not match with the reference genome.
Finally, remove sites containing missing data with 'bcftools view'.


```
# Generate a chromosome renaming file
for CHR in {1..23} X ; do 
    echo ${CHR} chr${CHR}
done >> chr_names.txt

# Multiple processing commands piped together
for CHR in {1..22} X; do
    bcftools annotate --rename-chrs chr_names.txt \
        ALL.chr${CHR}_GRCh38.genotypes.20170504.vcf.gz -Ou | \
    bcftools view -e 'INFO/AC<3 | INFO/AN-INFO/AC<3' -Ou | \
    bcftools norm -m -any -Ou | \
    bcftools view -i 'INFO/VT="SNP" | INFO/VT="INDEL"' -Ou | \
    bcftools norm -f Homo_sapiens_assembly38.fasta -d none -Ou | \
    bcftools view -m 2 -M 2 -Ou | \
    bcftools view -g ^miss -Oz -o 1000GP_chr${CHR}.vcf.gz"
done
```
 
If multiallelic sites are present in your data, in order to preserve them throughout the protocol, set ID field with unique IDs e.g. in format CHR_POS_REF_ALT. (RSIDs might contain duplicates, when the multiallelic sites are decomposed.)

```
for CHR in {1..23}; do
    bcftools annotate \
    --set-id '%CHROM\_%POS\_%REF\_%ALT' \
    panel_chr${CHR}.vcf.gz  \
    -Oz -o panel_SNPID_chr${CHR}.vcf.gz
done
```

1.3.3 Convert haploid genotypes to homozygous diploids
Often chrX is represented as haploid genotypes for males, however, Beagle can only handle diploid genotypes. 
The command here produces unphased diploid genotypes. But since the haploid genotypes are in diploid format as REF/REF or ALT/ALT, we can simply set the phase for those alleles with a simple sed replacement. 

The chrX ploidy can be corrected as follow

```
# Fix the chromosome X ploidy to phased diploid
# Requires a ploidy.txt file containing 
# space-separated CHROM,FROM,TO,SEX,PLOIDY 
echo "chrX 1 156040895 M 2" > ploidy.txt
bcftools +fixploidy \
    1000GP_chrX.vcf.gz -Ov -- -p ploidy.txt | \
    sed 's#0/0#0\|0#g;s#1/1#1\|1#g' | \
bcftools view -Oz -o 1000GP_chr23.vcf.gz"
```

1.3.4 Duplicate ID removal
Remove duplicate IDs. If you wish to preserve all multiallelic sites, replace the ID column with a unique ID e.g. CHR_POS_REF_ALT (as indicated in Step 1.3.2). 

Here, 1000GP did not contain multiallelic sites after AC filtering, and thus, RSIDs were preserved in the ID column. And since RSIDs are not always unique, duplicates should be removed.

```
for CHR in {1..23}; do
    bcftools query -f '%ID\n' 1000GP_chr${CHR}.vcf.gz | \
    sort | uniq -d > 1000GP_chr${CHR}.dup_id

    if [[ -s 1000GP_chr${CHR}.dup_id ]]; then
    	bcftools view -e ID=@1000GP_chr${CHR}.dup_id \
    	1000GP_chr${CHR}.vcf.gz \
        -Oz -o 1000GP_filtered_chr${CHR}.vcf.gz
    else 
    	mv 1000GP_chr${CHR}.vcf.gz \
        1000GP_filtered_chr${CHR}.vcf.gz
    fi
done
```
1.3.5 Reference panel allele frequencies
Generate a tab-delimited file of the reference panel allele frequencies, one variant per line, with columns CHR, SNP (in format CHR_POS_REF_ALT), REF, ALT, AF (including the header line).

First, update (or add) AF values in the INFO field, calculate it with BCFtools plugin +fill-tags:
```
# Check if the VCF does NOT contain AF in the INFO field,
# and calculate it with bcftools +fill-tags plugin

for CHR in {1..23}; do
    bcftools +fill-tags \
        1000GP_filtered_chr${CHR}.vcf.gz \
        -Oz -o 1000GP_AF_chr${CHR}.vcf.gz -- -t AF
done
```
Extract the wanted fields from each VCF file and combine as a single output file with the header:

```
# Generate a tab-delimited header
echo -e 'CHR\tSNP\tREF\tALT\tAF' \
    > 1000GP_imputation_all.frq

# Query the required fields from the VCF file 
# and append to the allele frequency file 
for CHR in {1..23}; do
    bcftools query \
    -f '%CHROM\t%CHROM\_%POS\_%REF\_%ALT\t%REF\t%ALT\t%INFO/AF\n' \
    1000GP_AF_chr${CHR}.vcf.gz \
    >> 1000GP_imputation_all.frq
done
```
Note: Chromosome notation in the panel.frq file should follow the GRCh38/hg38 notations ('chr#' for autosomal chromosomes and 'chrX' for chromosome 23).
 


1.3.6 Create binary reference panel files
The phased reference panel files per chromosome are required in bref format (bref = binary reference). For more information, see Beagle documentation at Beagle site:
https://faculty.washington.edu/browning/beagle/bref.16Dec15.pdf.
 
The required bref.*.jar is downloaded from Beagle site:
```
wget https://faculty.washington.edu/browning/beagle/bref.27Jan18.7e1.jar
```
Use the processed imputation reference panel VCFs as inputs for the example command below. 
The output files have the suffix '.bref' instead of '.vcf.gz'.
```
# Convert each file to bref format
for CHR in {1..23}; do
    java -jar /path/to/bref.27Jan18.7e1.jar \
        1000GP_AF_chr${CHR}.vcf.gz
done
```

1.3.7 Generate a list of the reference panel sample IDs
List of sample IDs present in the reference panel, one line per sample ID can be generated from any of the VCF files as in the example below (assuming that all chromosomes contain the same set of samples):
```
bcftools query -l 1000GP_AF_chr22.vcf.gz \
    > 1000GP_sample_IDs.txt
```
1.4. You are ready to start! 
As the last prepatory step, let's go over the required input data file(s) and also expected final output files!
 
1.4.1 Input file:
 
Post-QC chip genotype data in VCFv4.2 format and chrX genotypes as diploid genotypes:
DATASET.vcf.gz
 
Note: Chromosome notation should follow the GRCh38/hg38 notations (e.g. 'chr#' for autosomal chromosomes, 'chrX', 'chrY' and 'chrM').

Note: If the input data was lifted over from an older genome build to build version 38, cautious inspection of the data is highly recommended before proceeding with the protocol.

Note: If chrX is represented as haploid genotypes, follow step 1.3.3 first 'bcftools +fixploidy' command (not the other two piped commands) to convert to diploid genotypes.
 
 
1.4.2 Final output files:
DATASET_imputed_info_chr#.vcf.gz (where # is chromosome number)
DATASET_postimputation_summary_plots.pdf

Note: Several intermediate files are created during the protocol. Those files can be used for troubleshooting and deleted once the successful imputation is confirmed.



