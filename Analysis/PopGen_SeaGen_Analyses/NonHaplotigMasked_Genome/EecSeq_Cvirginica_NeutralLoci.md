# Filtering the neutral loci

Author: A. Zyck

Date: April 14, 2020

**Here, I am taking the `neutralloci.recode.vcf` file and filtering it for different regions of the genome:**
- Exon regions
- Genic regions
- Coding sequences
- Intron regions
- Intergenic regions

**First, I am creating a new working directory.**

In PATH: `/home/azyck/NB_capture/NB_ddocent/NB_neutralloci`

```
$ mkdir NB_neutralloci
$ cd NB_neutralloci
```

**Copy the `neutralloci.recode.vcf` file from the `NB_OutlierDetect` directory in which it was created.**

```
$ cp ../NB_OutlierDetect/neutralloci.recode.vcf .
```

**Link the genome region bed files to `NB_neutralloci`.**

The various oyster genome regions were identified in the [OysterGenomeProject](https://github.com/jpuritz/OysterGenomeProject). [J.Puritz](https://github.com/jpuritz) saved each region as a separate BED file in the `RAID_STORAGE2` directory in KITT:

- Exons: `sorted.ref3.0.exon.sc.bed`
- Genes: `sorted.ref3.0.gene.sc.bed`
- Coding sequences: `sorted.ref3.0.CDS.sc.bed`
- Introns: `cv.ref3.intron.sc.bed`
- Intergenic regions: `cv.ref3.intergenic.sc.bed`

```
# Linking BED files to working directory
$ ln -s /RAID_STORAGE2/Shared_Data/Oyster_Genome/sorted.ref3.0.exon.sc.bed .
$ ln -s /RAID_STORAGE2/Shared_Data/Oyster_Genome/sorted.ref3.0.gene.sc.bed .
$ ln -s /RAID_STORAGE2/Shared_Data/Oyster_Genome/sorted.ref3.0.CDS.sc.bed .
$ ln -s /RAID_STORAGE2/Shared_Data/Oyster_Genome/cv.ref3.intron.sc.bed .
$ ln -s /RAID_STORAGE2/Shared_Data/Oyster_Genome/cv.ref3.intergenic.sc.bed .
```

**Create new VCF files for each genomic region**

Filtering using `vcftools` with flag `--bed` so that VCF file is filtered to include only entries contained with the BED file.

**_Exons_**

```
# VCF file containing all exon regions
$ vcftools --vcf neutralloci.recode.vcf --recode --recode-INFO-all --bed sorted.ref3.0.exon.sc.bed --out neutralexon

output:
VCFtools - 0.1.17
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf neutralloci.recode.vcf
	--recode-INFO-all
	--out neutralexon
	--recode
	--bed sorted.ref3.0.exon.sc.bed

After filtering, kept 40 out of 40 Individuals
Outputting VCF file...
  Read 337295 BED file entries.
After filtering, kept 50238 out of a possible 74520 Sites
Run Time = 22.00 seconds
```

**_Genes_**

```
# VCF file containing all genic regions
$ vcftools --vcf neutralloci.recode.vcf --recode --recode-INFO-all --bed sorted.ref3.0.gene.sc.bed --out neutralgene

output:
VCFtools - 0.1.17
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf neutralloci.recode.vcf
	--recode-INFO-all
	--out neutralgene
  --recode
  --bed sorted.ref3.0.gene.sc.bed

After filtering, kept 40 out of 40 Individuals
Outputting VCF file...
	Read 35016 BED file entries.
After filtering, kept 66939 out of a possible 74520 Sites
Run Time = 12.00 seconds
```

**_Coding Sequences_**

```
# VCF file containing all coding sequences
$ vcftools --vcf neutralloci.recode.vcf --recode --recode-INFO-all --bed sorted.ref3.0.CDS.sc.bed --out neutralCDS

output:
VCFtools - 0.1.17
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf neutralloci.recode.vcf
	--recode-INFO-all
	--out neutralCDS
	--recode
	--bed sorted.ref3.0.CDS.sc.bed

After filtering, kept 40 out of 40 Individuals
Outputting VCF file...
  Read 295525 BED file entries.
After filtering, kept 30107 out of a possible 74520 Sites
Run Time = 20.00 seconds
```

**_Introns_**

```
# VCF file containing all intron regions
$ vcftools --vcf neutralloci.recode.vcf --recode --recode-INFO-all --bed cv.ref3.intron.sc.bed --out neutralintron

output:
VCFtools - 0.1.17
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf neutralloci.recode.vcf
	--recode-INFO-all
	--out neutralintron
	--recode
	--bed cv.ref3.intron.sc.bed

After filtering, kept 40 out of 40 Individuals
Outputting VCF file...
	Read 302280 BED file entries.
After filtering, kept 16701 out of a possible 74520 Sites
Run Time = 21.00 seconds
```

**_Intergenic regions_**

```
# VCF file containing all intergenic regions
$ vcftools --vcf neutralloci.recode.vcf --recode --recode-INFO-all --bed cv.ref3.intergenic.sc.bed --out neutralintergenic

output:
VCFtools - 0.1.17
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf neutralloci.recode.vcf
	--recode-INFO-all
	--out neutralintergenic
	--recode
	--bed cv.ref3.intergenic.sc.bed

After filtering, kept 40 out of 40 Individuals
Outputting VCF file...
	Read 35025 BED file entries.
After filtering, kept 7581 out of a possible 74520 Sites
Run Time = 3.00 seconds
```

Population genomic analyses are performed on each VCF file and documented [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/tree/master/Analysis/PopGen_SeaGen_Analyses/NonHaplotigMasked_Genome/Neutral_SNPs/NB_PopGen_NeutralSNPs.Rmd)
