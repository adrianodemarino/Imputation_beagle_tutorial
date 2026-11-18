# Genotype imputation with Beagle 5.5

> [!TIP]
> There is a new tutorial for **[Selphi 2](https://github.com/omicsedge/selphi2)**, the phasing and imputation engine I work on:
> **[Imputation_selphi2_tutorial](https://github.com/adrianodemarino/Imputation_selphi2_tutorial)**. Same data, same steps,
> with Beagle and Selphi 2 side by side.

This guide imputes chromosome 22 for 100 people from a simulated genotyping array, using the 1000 Genomes
high-coverage reference panel and Beagle 5.5. Every command below was run as written; the numbers in the
"what you should see" boxes come from that run (16 threads, Linux).

1. [Install](#1-install)
2. [Download the panel and the genetic map](#2-download-the-panel-and-the-genetic-map)
3. [Prepare the panel, the test people and their "array"](#3-prepare-the-panel-the-test-people-and-their-array)
4. [Convert the panel to bref3](#4-convert-the-panel-to-bref3)
5. [Phase and impute](#5-phase-and-impute)
6. [Read the output](#6-read-the-output)
7. [Check the accuracy](#7-check-the-accuracy)
8. [Your own data](#8-your-own-data)

## 1. Install

You need Java 8 or later, bcftools 1.15 or later and the two Beagle jars.

```bash
conda create -n imputation -c conda-forge -c bioconda bcftools openjdk
conda activate imputation

mkdir -p tools
wget -P tools https://faculty.washington.edu/browning/beagle/beagle.27Feb25.75f.jar
wget -P tools https://faculty.washington.edu/browning/beagle/bref3.27Feb25.75f.jar
```

The newest jar is always listed on the [Beagle page](https://faculty.washington.edu/browning/beagle/beagle.html).

## 2. Download the panel and the genetic map

The panel is the 1000 Genomes high-coverage release (3,202 people sequenced at 30x by the New York Genome
Center, phased with SHAPEIT), on GRCh38. The map is Beagle's PLINK-format GRCh38 map.

```bash
KG=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage
PANEL=$KG/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz

mkdir -p data maps
wget -P data "$PANEL" "$PANEL.tbi"
wget -P data "$KG/20130606_g1k_3202_samples_ped_population.txt"

wget -P maps https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip
unzip -q maps/plink.GRCh38.map.zip -d maps
MAP=maps/chr_in_chrom_field/plink.chrchr22.GRCh38.map    # yes, "chrchr22": the file uses chr22 inside
```

The zip has two copies of the map: `chr_in_chrom_field/` (chromosome written `chr22`, as in GRCh38 VCFs) and
`no_chr_in_chrom_field/` (`22`). The chromosome name has to match your VCF.

## 3. Prepare the panel, the test people and their "array"

To measure accuracy we need people whose full genome we know but whom the panel has never seen. We take 100 of
them out of the panel (20 per continental group, none with a parent or child in 1000 Genomes), keep their 30x
genotypes as the truth, and give the imputation only what a genotyping chip would read.

[`test_samples.txt`](test_samples.txt) lists the 100 people.

```bash
mkdir -p work
IN=data/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz

# Panel: everyone else, biallelic SNPs and indels, sites still polymorphic without the 100
bcftools view -S ^test_samples.txt -v snps,indels -m2 -M2 "$IN" -Ou |
  bcftools view -c1 -Oz -o work/panel.chr22.vcf.gz
bcftools index -t work/panel.chr22.vcf.gz

# Truth: the 100 people at every site
bcftools view -S test_samples.txt -v snps,indels -m2 -M2 "$IN" -Oz -o work/truth.chr22.vcf.gz
bcftools index -t work/truth.chr22.vcf.gz

# Array: common SNPs (MAF >= 1%), at most one every 5 kb, genotypes unphased as a chip reports them
bcftools view -v snps -q 0.01:minor work/panel.chr22.vcf.gz -Ou |
  bcftools query -f '%CHROM\t%POS\n' |
  awk 'BEGIN{OFS="\t"} $2 >= last + 5000 {print; last = $2}' > work/array_sites.chr22.tsv
bcftools view -T work/array_sites.chr22.tsv work/truth.chr22.vcf.gz -Ou |
  bcftools +setGT -Ou - -- -t a -n u |
  bcftools annotate -x INFO -Oz -o work/array.chr22.vcf.gz
bcftools index -t work/array.chr22.vcf.gz
```

What you should see (`bcftools index -n`, `bcftools query -l | wc -l`):

```text
panel    3102 people    1064451 sites
truth     100 people    1065095 sites
array     100 people       7216 sites
```

## 4. Convert the panel to bref3

bref3 is Beagle's compressed panel format. It is built once and makes every run after it faster.

```bash
java -jar tools/bref3.27Feb25.75f.jar work/panel.chr22.vcf.gz > work/panel.chr22.bref3
```

The 445 MB panel VCF becomes a 196 MB `panel.chr22.bref3`.

## 5. Phase and impute

Beagle 5 phases and imputes in one run. `gt=` is the target, `ref=` the panel, `map=` the genetic map.

```bash
mkdir -p results
java -Xmx16g -jar tools/beagle.27Feb25.75f.jar \
  gt=work/array.chr22.vcf.gz ref=work/panel.chr22.bref3 map="$MAP" \
  out=results/beagle.chr22 nthreads=16
bcftools index -t results/beagle.chr22.vcf.gz
```

What you should see at the end of `results/beagle.chr22.log`:

```text
Reference samples:                3,102
Study     samples:                  100
...
Imputation time:               3 seconds
```

The whole command took 13.4 seconds and 4.9 GB of memory (the Java process, with `-Xmx16g`) on 16 threads,
median of three runs measured with `/usr/bin/time`.

Useful options: `gp=true` adds genotype probabilities, `ap=true` the per-haplotype ALT probabilities,
`impute=false` only phases. `-Xmx` must be large enough for the panel; with too little Beagle stops with an
out-of-memory error.

## 6. Read the output

```bash
# one imputed site and one chip site, first person only
bcftools view -H -i 'INFO/IMP=1 && INFO/AF>0.2' results/beagle.chr22.vcf.gz | head -1 | cut -f1-10
bcftools view -H -e 'INFO/IMP=1' results/beagle.chr22.vcf.gz | head -1 | cut -f1-10
```

```text
chr22  10562747  22:10562747:AGTTTT:A  AGTTTT  A  .  PASS  DR2=0.36;AF=0.2799;IMP  GT:DS  0|0:0.21
chr22  10519276  22:10519276:G:C       G       C  .  PASS  DR2=1;AF=0.025         GT:DS  0|0:0
```

- `GT` is phased (`0|1`): the first allele is on one parental haplotype, the second on the other.
- `DS` is the dosage, the expected number of ALT alleles (0 to 2).
- `INFO/DR2` is the estimated squared correlation between the imputed dosage and the true genotype.
  `DR2 < 0.3` is the usual filter for poorly imputed sites.
- `INFO/IMP` marks imputed sites; sites without it were on the array.

Keep only well-imputed sites:

```bash
bcftools view -i 'INFO/DR2>=0.3' results/beagle.chr22.vcf.gz -Oz -o results/beagle.chr22.dr2.vcf.gz
```

## 7. Check the accuracy

We compare the imputed dosages with the 30x truth, on imputed sites only (the array sites are left out), and
report the squared correlation r² by minor allele frequency. The scorer is in Selphi 2 and works on any
imputed VCF (see [its tutorial](https://github.com/adrianodemarino/Imputation_selphi2_tutorial#install)):

```bash
selphi --evaluate results/beagle.chr22.vcf.gz --truth work/truth.chr22.vcf.gz \
  --exclude-sites work/array.chr22.vcf.gz --out results/beagle.eval
```

```text
MAF bin                N variants    Mean R²  Concordance
------------------------------------------------------------
0.5-1%                     133440     0.4849       0.9927
1-2%                        75037     0.5559       0.9862
2-5%                        59848     0.6604       0.9753
5-10%                       33375     0.7302       0.9548
10-20%                      32447     0.7539       0.9231
20-50%                      59375     0.7997       0.8928

Per-sample (n=100):
  R²:          mean=0.893349, min=0.838404, max=0.953355
```

r² by frequency bin is the squared correlation, over all sites in the bin, between imputed dosage and true
genotype. The bins start at 0.5% because 100 people carry 200 copies of chr22. Selphi 2 scores higher in every
bin on the same data; the [Selphi 2 tutorial](https://github.com/adrianodemarino/Imputation_selphi2_tutorial#step-4-compare-with-beagle)
has the side-by-side.

Rare variants are where imputation struggles: a site carried by a handful of people in the panel can only
be imputed when the target shares a long stretch of haplotype with one of them.

## 8. Your own data

Before imputing your own genotypes:

- **Same build and chromosome names** in target, panel and map (`chr22` vs `22`). Rename with
  `bcftools annotate --rename-chrs`. Data lifted over from GRCh37 needs a careful check.
- **Same strand and alleles.** Arrays often report the other strand. Normalise against the reference genome:
  `bcftools norm -f GRCh38.fa -c ws` fixes swapped REF/ALT, and drop palindromic A/T and C/G SNPs if you are
  not sure of the strand.
- **Biallelic records.** Split multiallelic sites with `bcftools norm -m -any`.
- **Sample QC first.** Remove samples and sites with a high missing rate, and sites far from Hardy-Weinberg
  equilibrium, before imputing.
- **One chromosome per run**, or pass `chrom=chr22` to pick one.
- **A panel close to your people.** Accuracy depends more on the panel than on the tool: the larger and the
  closer to your population, the better.

## Citation

Browning BL, Tian X, Zhou Y, Browning SR. Fast two-stage phasing of large-scale sequence data.
*Am J Hum Genet* 2021;108:1880-1890.

Byrska-Bishop M et al. High-coverage whole-genome sequencing of the expanded 1000 Genomes Project cohort
including 602 trios. *Cell* 2022;185:3426-3440.
