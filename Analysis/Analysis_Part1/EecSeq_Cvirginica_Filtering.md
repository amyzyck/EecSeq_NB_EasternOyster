# *Crassostrea virginica* EecSeq Seascape Genomics Analysis - SNP Filtering

Author: Amy Zyck

Last Edited: March 12, 2021

EecSeq reads were demultiplexed and underwent quality control, alignment to the Eastern Oyster genome (+ Haplotig Masked Genome), and variant calling for Single Nucleotide Polymorphisms (SNPs). The bioinformatic steps can be accessed [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Analysis/Analysis_Part1/EecSeq_Cvirginica_dDocent.md)

# Variant Filtering

#### For the following filtering steps, each step will be completed twice, once on the original VCF file and again on the masked haplotig VCF file. The outputs for each filtering step on each VCF file will be presented in that order.

First, I created a new directory and linked the VCF file to this new directory for organization purposes.

```
# Within NB_ddocent Directory
$ mkdir NB_SNPFiltering
$ cd NB_SNPFiltering

# Link Raw VCF file to new directory
$ ln -s ../TotalRawSNPs.vcf.gz .
```

Steps taken for filtering are based on steps completed by [J. Puritz](https://github.com/jpuritz) for the [OysterGenomeProject](https://github.com/jpuritz/OysterGenomeProject/blob/master/Bioinformatics/OysterGenome.md) and steps I completed with [RADSeq data for the fiddler crab, _Uca rapax_](https://github.com/amyzyck/RADseq_Uca-rapax_2016/blob/master/Scripts/Filtering/Filtering_UcaRapax.md).

Scripts used for filtering steps are located in the [dDocent repository](https://github.com/jpuritz/dDocent/tree/master/scripts) and the [OysterGenomeProject repository](https://github.com/jpuritz/OysterGenomeProject/tree/master/Bioinformatics/Scripts) on Github.

**First step for filtering is to remove sites with a QUALITY score below 20 and those sites that have more than 80% missing data.**

```
$ vcftools --gzvcf TotalRawSNPs.vcf.gz --minQ 20 --recode --recode-INFO-all --out TRS --max-missing 0.20

output:
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--gzvcf TotalRawSNPs.vcf.gz
	--recode-INFO-all
	--minQ 20
	--max-missing 0.2
	--out TRS
	--recode

Using zlib version: 1.2.11
After filtering, kept 50 out of 50 Individuals
Outputting VCF file...
After filtering, kept 37412226 out of a possible 52332265 Sites
Run Time = 8143.00 seconds
```

```
Masked haplotig output:
After filtering, kept 50 out of 50 Individuals
Outputting VCF file...
After filtering, kept 37409463 out of a possible 48598546 Sites
Run Time = 8075.00 seconds
```

**Next step is to mark any genotypes with less than 5 reads as missing.**

```
$ vcftools --vcf TRS.recode.vcf --minDP 5 --recode --recode-INFO-all --out TRSdp5

output:
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf TRS.recode.vcf
	--recode-INFO-all
	--minDP 5
	--out TRSdp5
	--recode

After filtering, kept 50 out of 50 Individuals
Outputting VCF file...
After filtering, kept 37412226 out of a possible 37412226 Sites
Run Time = 12796.00 seconds
```

```
Masked haplotig output:
After filtering, kept 50 out of 50 Individuals
Outputting VCF file...
After filtering, kept 37409463 out of a possible 37409463 Sites
Run Time = 12603.00 seconds
```

**Next, remove sites with more than 50% missing data.**

```
$ vcftools --vcf TRSdp5.recode.vcf --max-missing 0.5 --recode --recode-INFO-all --out TRSdp5g5

output:
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf TRSdp5.recode.vcf
	--recode-INFO-all
	--maf 0.001
	--max-missing 0.5
	--out TRSdp5g5
	--recode

After filtering, kept 50 out of 50 Individuals
Outputting VCF file...
After filtering, kept 2859613 out of a possible 37412226 Sites
Run Time = 1714.00 seconds
```

```
Masked haplotig output:
After filtering, kept 50 out of 50 Individuals
Outputting VCF file...
After filtering, kept 3712479 out of a possible 37409463 Sites
Run Time = 1016.00 seconds
```

**MAF filtering for FreeBayes Output.**

```
$ curl -L -O https://raw.githubusercontent.com/jpuritz/dDocent/master/scripts/untested/multi.maf.sh
$ chmod +x multi.maf.sh
$ multi.maf.sh TRSdp5g5.recode.vcf 0.001 TRSdp5g5maf

output:
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf TRSdp5g5.recode.vcf
	--recode-INFO-all
	--out TRSdp5g5maf
	--positions maf.loci.to.keep
	--recode

After filtering, kept 50 out of 50 Individuals
Outputting VCF file...
After filtering, kept 3122748 out of a possible 3148322 Sites
Run Time = 526.00 seconds
```

```
Masked haplotig output:
After filtering, kept 50 out of 50 Individuals
Outputting VCF file...
After filtering, kept 3686076 out of a possible 3712479 Sites
Run Time = 782.00 seconds
```

**Use a custom script called `filter_missing_ind.sh` to filter out bad individuals.**

Here I mostly just want to see % of missing data within individuals. We are working with a small number of individuals, so I don't want to remove any.

```
$ curl -L -O https://raw.githubusercontent.com/jpuritz/dDocent/master/scripts/filter_missing_ind.sh
$ chmod +x filter_missing_ind.sh
$ ./filter_missing_ind.sh TRSdp5g5maf.recode.vcf TRSdp5g5mafMI

output:
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf TRSdp5g5maf.recode.vcf
	--missing-indv
	--out TRSdp5g5mafMI

After filtering, kept 50 out of 50 Individuals
Outputting Individual Missingness
After filtering, kept 3122748 out of a possible 3122748 Sites
Run Time = 31.00 seconds

Histogram of % missing data per individual

        5 +---------------------------------------------------------------------------------------------------------+
          |        * * +  **         +            +            +            +             +            +            |
          |        * *    **                                  'totalmissing' using (bin($1,binwidth)):(1.0) ******* |
          |        * *    **                                                                                        |
          |        * *    **                                                                                        |
        4 |-+     ** *    **                                                                                      +-|
          |       ** *    **                                                                                        |
          |       ** *    **                                                                                        |
          |       ** *    **                                                                                        |
        3 |-+  ***** *    **       ******                                           *********************         +-|
          |    *  ** *    **       *    *                                           *                   *           |
          |    *  ** *    **       *    *                                           *                   *           |
          |    *  ** *    **       *    *                                           *                   *           |
          |    *  ** *    **       *    *                                           *                   *           |
        2 |-+  *  ** * ******* *** *    ***************                             *                   *         +-|
          |    *  ** * * *** * * * *    *             *                             *                   *           |
          |    *  ** * * *** * * * *    *             *                             *                   *           |
          |    *  ** * * *** * * * *    *             *                             *                   *           |
        1 |*****  ** *** *** *** ***    *             *******************************                   ***       +-|
          |*   *  ** *** *** * * * *    *             *                             *                   ***         |
          |*   *  ** *** *** * * * *    *             *                             *                   ***         |
          |*   *  ** *** *** * * * *    *             *                             *                   ***         |
          |*   *  ** *** *** * * * * +  *         +   *        +            +       *     +            +***         |
        0 +---------------------------------------------------------------------------------------------------------+
         0.1          0.2           0.3          0.4          0.5          0.6           0.7          0.8          0.9
                                                       % of missing data

The 85% cutoff would be 0.506825
Would you like to set a different cutoff, yes or no
$ yes #a 85% cutoff would remove 5 individuals
Please enter new cutoff
$ 0.95
All individuals with more than 95.0% missing data will be removed.

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf TRSdp5g5maf.recode.vcf
	--remove lowDP.indv
	--recode-INFO-all
	--out TRSdp5g5mafMI
	--recode

Excluding individuals in 'exclude' list
After filtering, kept 50 out of 50 Individuals
Outputting VCF file...
After filtering, kept 3122748 out of a possible 3122748 Sites
Run Time = 523.00 seconds
```

```
Masked haplotig output:
After filtering, kept 50 out of 50 Individuals
Outputting VCF file...
After filtering, kept 3686076 out of a possible 3686076 Sites
Run Time = 772.00 seconds
```

Here is a list of the worst five samples:

```
$ cat <(head -1 TRSdp5g5mafMI.imiss ) <(mawk '!/F_MI/' TRSdp5g5mafMI.imiss | sort -k5 -r -n ) | head -6 | column -t

INDV   N_DATA   N_GENOTYPES_FILTERED  N_MISS   F_MISS
PVD_4  3122748  0                     2579761  0.826119
PVD_5  3122748  0                     2545375  0.815107
PVD_3  3122748  0                     2526932  0.809201
PVD_1  3122748  0                     2520712  0.80721
PVD_2  3122748  0                     2513993  0.805058
```

These 5 samples had very few reads compared to the other samples from the beginning.

**Use a second custom `script pop_missing_filter.sh` to filter loci that have high missing data values in a single population.**

This step needs a file that maps individuals to populations `popmap`.

```
$ head popmap

PVD_1	PVD
PVD_2	PVD
PVD_3	PVD
PVD_4	PVD
PVD_5	PVD
PVD_6	PVD
PVD_7	PVD
PVD_8	PVD
PVD_9	PVD
PVD_10	PVD
```

```
$ curl -L -O https://raw.githubusercontent.com/jpuritz/dDocent/master/scripts/pop_missing_filter.sh
$ chmod +x pop_missing_filter.sh
$ ./pop_missing_filter.sh TRSdp5g5mafMI.recode.vcf popmap 0.1 1 TRSdp5g5mafMIp9

sample output:
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf TRSdp5g5mafMI.recode.vcf
	--keep keep.BIS
	--out BIS
	--missing-site

Keeping individuals in 'keep' list
After filtering, kept 10 out of 50 Individuals
Outputting Site Missingness
After filtering, kept 3122748 out of a possible 3122748 Sites
Run Time = 41.00 seconds

# Repeated for other 4 populations

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf TRSdp5g5mafMI.recode.vcf
	--exclude-positions loci.to.remove
	--recode-INFO-all
	--out TRSdp5g5mafMIp9
	--recode

After filtering, kept 50 out of 50 Individuals
Outputting VCF file...
After filtering, kept 213068 out of a possible 3122748 Sites
Run Time = 51.00 seconds
```

```
Masked haplotig output:
After filtering, kept 50 out of 50 Individuals
Outputting VCF file...
After filtering, kept 259970 out of a possible 3686076 Sites
Run Time = 78.00 seconds
```

**Filter out any sites with less than 95% overall call rate and MAF of 0.001.**

```
$ vcftools --vcf TRSdp5g5mafMIp9.recode.vcf --recode-INFO-all --max-missing 0.95 --maf 0.001 --out TRSdp5g5mafMIp9g9 --recode

output:
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf TRSdp5g5mafMIp9.recode.vcf
	--recode-INFO-all
	--maf 0.001
	--max-missing 0.95
	--out TRSdp5g5mafMIp9g9
	--recode

After filtering, kept 207269 out of a possible 213068 Sites
Run Time = 35.00 seconds
```

```
Masked haplotig output:
After filtering, kept 50 out of 50 Individuals
Outputting VCF file...
After filtering, kept 254116 out of a possible 259970 Sites
Run Time = 58.00 seconds
```

**Next, split the VCF files into nDNA and mtDNA for further filtering**

```
$ cat TRSdp5g5mafMIp9g9.recode.vcf | mawk '!/NC_007175.2/' > TRSdp5g5mafMIp9g9.nDNA.vcf

$ cat TRSdp5g5mafMIp9g9.recode.vcf | head -10000 | mawk '/#/' > header
$ cat TRSdp5g5mafMIp9g9.recode.vcf | head -10000 | mawk '/#/' > header
$ cat header <(cat TRSdp5g5mafMIp9g9.recode.vcf |head -20000| mawk '/NC_007175.2/') > TRSdp5g5mafMIp9g9.mtDNA.vcf
```

The remaining filtering steps will be split between the mtDNA and nDNA

#### mtDNA

**Break complex mutational events (combinations of SNPs and INDELs) into separate SNP and INDEL calls, and then remove INDELs.**

```
$ vcffilter -s -f "AB < 0.001" TRSdp5g5mafMIp9g9.mtDNA.vcf | vcffilter -s -f "QUAL / DP > 0.25" > TRSdp5g5mafMIp9g9.mtDNA.F.vcf
$ vcfallelicprimitives -k -g TRSdp5g5mafMIp9g9.mtDNA.F.vcf | sed 's:\.|\.:\.\/\.:g' > TRSdp5g5mafMIp9g9mtDNAF.prim
$ vcftools --vcf TRSdp5g5mafMIp9g9mtDNAF.prim --remove-indels --recode --recode-INFO-all --out SNP.TRSdp5g5mafMIp9g9mtDNA

output:
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf TRSdp5g5mafMIp9g9mtDNAF.prim
	--recode-INFO-all
	--out SNP.TRSdp5g5mafMIp9g9mtDNA
	--recode
	--remove-indels

After filtering, kept 50 out of 50 Individuals
Outputting VCF file...
After filtering, kept 469 out of a possible 486 Sites
Run Time = 0.00 seconds
```

```
Masked haplotig output:
After filtering, kept 50 out of 50 Individuals
Outputting VCF file...
After filtering, kept 469 out of a possible 486 Sites
Run Time = 0.00 seconds
```

#### nDNA

**Use a custom filter script `dDocent_ngs_filters`.**


```
$ curl -L -O https://raw.githubusercontent.com/jpuritz/OysterGenomeProject/master/Bioinformatics/Scripts/dDocent_ngs_filters
$ chmod +x dDocent_ngs_filters
$ ./dDocent_ngs_filters TRSdp5g5mafMIp9g9.nDNA.vcf TRSdp5g5mafMIp9g9dnDNA

output:
This script will automatically filter a FreeBayes generated VCF file using criteria related to site depth,
quality versus depth, allelic balance at heterzygous individuals, and paired read representation.
The script assumes that loci and individuals with low call rates (or depth) have already been removed.

Contact Jon Puritz (jpuritz@gmail.com) for questions and see script comments for more details on particular filters

Is this from a mixture of SE and PE libraries? Enter yes or no.
$ no
Number of sites filtered based on allele balance at heterozygous loci, locus quality, and mapping quality / Depth
 32795 of 206750

Number of additional sites filtered based on properly paired status
 90 of 173955

Number of sites filtered based on high depth and lower than 2*DEPTH quality score
 27245 of 173955



                                               Histogram of mean depth per site

     7000 +---------------------------------------------------------------------------------------------------------+
          |  +     +     +     +     +    +     +     +     +     +     +     +     +     +    +     +     +     +  |
          |               ****                            'meandepthpersite' using (bin($1,binwidth)):(1.0) ******* |
     6000 |-+            ** *****                                                                                 +-|
          |            **** ** *****                                                                                |
          |            * ** ** ** *****                                                                             |
          |           ** ** ** ** ** **                                                                             |
     5000 |-+       **** ** ** ** ** ***                                                                          +-|
          |         * ** ** ** ** ** *****                                                                          |
          |         * ** ** ** ** ** *** ****                                                                       |
     4000 |-+      ** ** ** ** ** ** *** ** **                                                                    +-|
          |        ** ** ** ** ** ** *** ** ****                                                                    |
          |      **** ** ** ** ** ** *** ** ** **                                                                   |
     3000 |-+    * ** ** ** ** ** ** *** ** ** *****                                                              +-|
          |     ** ** ** ** ** ** ** *** ** ** ** *****                                                             |
          |     ** ** ** ** ** ** ** *** ** ** ** ** ****                                                           |
     2000 |-+   ** ** ** ** ** ** ** *** ** ** ** ** ** ****                                                      +-|
          |   **** ** ** ** ** ** ** *** ** ** ** ** ** ** *                                                        |
          |   * ** ** ** ** ** ** ** *** ** ** ** ** ** ** *                                                        |
          |   * ** ** ** ** ** ** ** *** ** ** ** ** ** ** *                                                        |
     1000 |-+** ** ** ** ** ** ** ** *** ** ** ** ** ** ** **                                                     +-|
          |**** ** ** ** ** ** ** ** *** ** ** ** ** ** ** ************                                             |
          |* ** ** ** ** ** ** ** ** *** ** ** ** ** ** ** ** *** ** ***********************************************|
        0 +---------------------------------------------------------------------------------------------------------+
             12    16    20    24    28   32    36    40    44    48    52    56    60    64   68    72    76    80
                                                          Mean Depth

The 95% cutoff would be 66
Would you like to use a different maximum mean depth cutoff than 66, yes or no
$ no
Maximum mean depth cutoff is 66
Number of sites filtered based on maximum mean depth
 7488 of 146627

Total number of sites filtered
 67611 of 206750

Remaining sites
 139139

Filtered VCF file is called Output_prefix.FIL.recode.vcf

Filter stats stored in TRSdp5g5mafMIp9g9dnDNA.filterstats
```

```
Masked haplotig output:
Number of sites filtered based on maximum mean depth
 9305 of 180550

Total number of sites filtered
 82353 of 253598

Remaining sites
 171245

Filtered VCF file is called Output_prefix.FIL.recode.vcf

Filter stats stored in TRSdp5g5mafMIp9g9dnDNA.filterstats
```

**Break complex mutational events (combinations of SNPs and INDELs) into separate SNP and INDEL calls, and then remove INDELs.**

```
$ vcfallelicprimitives TRSdp5g5mafMIp9g9dnDNA.FIL.recode.vcf --keep-info --keep-geno > TRSdp5g5mafMIp9g9dnDNA.prim.vcf
$ vcftools --vcf TRSdp5g5mafMIp9g9dnDNA.prim.vcf --remove-indels --recode --recode-INFO-all --out SNP.TRSdp5g5mafMIp9g9dnDNA

output:
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf TRSdp5g5mafMIp9g9dnDNA.prim.vcf
	--recode-INFO-all
	--out SNP.TRSdp5g5mafMIp9g9dnDNA
	--recode
	--remove-indels

After filtering, kept 50 out of 50 Individuals
Outputting VCF file...
After filtering, kept 172463 out of a possible 196124 Sites
Run Time = 40.00 seconds
```

```
Masked haplotig output:
After filtering, kept 50 out of 50 Individuals
Outputting VCF file...
After filtering, kept 210181 out of a possible 238353 Sites
Run Time = 36.00 seconds
```

***
**Filtered by a few different Minor Allele Frequencies**

Minor Allele Frequency greater than 1%

```
$ vcftools --vcf SNP.TRSdp5g5mafMIp9g9dnDNA.recode.vcf --maf 0.01 --recode --recode-INFO-all --out SNP.TRSdp5g5mafMIp9g9dnDNAmaf01

output:
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf SNP.TRSdp5g5mafMIp9g9dnDNA.recode.vcf
	--recode-INFO-all
	--maf 0.01
	--out SNP.TRSdp5g5mafMIp9g9dnDNAmaf01
	--recode
After filtering, kept 50 out of 50 Individuals
Outputting VCF file...
After filtering, kept 172463 out of a possible 172463 Sites
Run Time = 40.00 second
```

```
Masked haplotig output:
After filtering, kept 50 out of 50 Individuals
Outputting VCF file...
After filtering, kept 210181 out of a possible 210181 Sites
Run Time = 39.00 seconds
```

Minor Allele Frequency greater than 2.5%

```
$ vcftools --vcf SNP.TRSdp5g5mafMIp9g9dnDNA.recode.vcf --maf 0.025 --recode --recode-INFO-all --out SNP.TRSdp5g5mafMIp9g9dnDNAmaf025

output:
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf SNP.TRSdp5g5mafMIp9g9dnDNA.recode.vcf
	--recode-INFO-all
	--maf 0.025
	--out SNP.TRSdp5g5mafMIp9g9dnDNAmaf025
	--recode
After filtering, kept 50 out of 50 Individuals
Outputting VCF file...
After filtering, kept 101484 out of a possible 172463 Sites
Run Time = 25.00 seconds
```

```
Masked haplotig output:
After filtering, kept 50 out of 50 Individuals
Outputting VCF file...
After filtering, kept 123524 out of a possible 210181 Sites
Run Time = 23.00 seconds
```

Minor Allele Frequency greater than 5%

```
$ vcftools --vcf SNP.TRSdp5g5mafMIp9g9dnDNA.recode.vcf --maf 0.05 --recode --recode-INFO-all --out SNP.TRSdp5g5mafMIp9g9dnDNAmaf05

output:
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf SNP.TRSdp5g5mafMIp9g9dnDNA.recode.vcf
	--recode-INFO-all
	--maf 0.05
	--out SNP.TRSdp5g5mafMIp9g9dnDNAmaf05
	--recode
After filtering, kept 50 out of 50 Individuals
Outputting VCF file...
After filtering, kept 77419 out of a possible 172463 Sites
Run Time = 20.00 seconds
```

```
Masked haplotig output:
After filtering, kept 50 out of 50 Individuals
Outputting VCF file...
After filtering, kept 94140 out of a possible 210181 Sites
Run Time = 19.00 seconds
```

**Max 2 alleles**

```
$ vcftools --vcf SNP.TRSdp5g5mafMIp9g9dnDNAmaf05.recode.vcf --recode --recode-INFO-all --out SNP.TRSdp5g5mafMIp9g9dnDNAmaf052A --max-alleles 2

output:
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf SNP.TRSdp5g5mafMIp9g9dnDNAmaf05.recode.vcf
	--recode-INFO-all
	--max-alleles 2
	--out SNP.TRSdp5g5mafMIp9g9dnDNAmaf052A
	--recode
After filtering, kept 50 out of 50 Individuals
Outputting VCF file...
After filtering, kept 75402 out of a possible 77419 Sites
Run Time = 15.00 seconds
```

```
Masked haplotig output:
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf SNP.TRSdp5g5mafMIp9g9dnDNAmaf05.recode.vcf
	--recode-INFO-all
	--max-alleles 2
	--out SNP.TRSdp5g5mafMIp9g9dnDNAmaf052Ahap
	--recode
After filtering, kept 50 out of 50 Individuals
Outputting VCF file...
After filtering, kept 91702 out of a possible 94140 Sites
Run Time = 15.00 seconds
```

**Continue on to Outlier Detection using `SNP.TRSdp5g5mafMIp9g9dnDNAmaf052A.recode.vcf` for the original genome**

**Continue on to Outlier Detection using `SNP.TRSdp5g5mafMIp9g9dnDNAmaf052Ahap.recode.vcf` for the masked haplotig genome**

Steps for Outlier Detection can be accessed [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Analysis/Analysis_Part1/EecSeq_Cvirginica_OutlierDetection.md)
