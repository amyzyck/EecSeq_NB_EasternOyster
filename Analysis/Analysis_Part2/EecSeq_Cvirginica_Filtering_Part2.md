# *Crassostrea virginica* EecSeq Seascape Genomics Re - Analysis - SNP Filtering

Author: Amy Zyck

Last Edited: February 22, 2023

EecSeq reads were demultiplexed and underwent quality control, alignment to the Eastern Oyster genome (Haplotig Masked Genome), and variant calling for Single Nucleotide Polymorphisms (SNPs). The bioinformatic steps can be accessed [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Analysis/Analysis_Part2/EecSeq_Cvirginica_dDocent_Part2.md).

# Variant Filtering

First, I created a new directory and linked the VCF file to this new directory for organization purposes.

```
# Within NB_ddhaplo Directory
$ mkdir NB_SNPFiltering_both
$ cd NB_SNPFiltering_both

# Link Raw VCF file to new directory
$ ln -s ../TotalRawSNPs.vcf.gz .
```

Steps taken for filtering are based on steps completed by [J. Puritz](https://github.com/jpuritz) for the [OysterGenomeProject](https://github.com/jpuritz/OysterGenomeProject/blob/master/Bioinformatics/OysterGenome.md) and steps I completed with [RADSeq data for the fiddler crab, _Uca rapax_](https://github.com/amyzyck/RADseq_Uca-rapax_2016/blob/master/Scripts/Filtering/Filtering_UcaRapax.md).

Scripts used for filtering steps are located in the [dDocent repository](https://github.com/jpuritz/dDocent/tree/master/scripts) and the [OysterGenomeProject repository](https://github.com/jpuritz/OysterGenomeProject/tree/master/Bioinformatics/Scripts) on Github.

Before officially starting, I am rerunning the VCF filtering that's a part of the dDocent_ngs script. I just want to confirm that the working VCF file applied that filtering step.

```
$ vcftools --gzvcf TotalRawSNPs.vcf.gz --mac 1 --max-non-ref-af 1 --minQ 30 --max-missing 0.9 --non-ref-af 0.001 --out Final
```

```
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--gzvcf TotalRawSNPs.vcf.gz
	--mac 1
	--max-non-ref-af 1
	--minQ 30
	--max-missing 0.9
	--non-ref-af 0.001
	--out Final

Using zlib version: 1.2.11

After filtering, kept 90 out of 90 Individuals
After filtering, kept 6585214 out of a possible 55253967 Sites
Run Time = 1375.00 seconds
```

**First step for filtering is to remove sites with a QUALITY score below 20 and those sites that have more than 80% missing data.**

```
$ vcftools --gzvcf TotalRawSNPs.vcf.gz --minQ 20 --max-missing 0.20 --recode --recode-INFO-all --out TRSboth
```

```
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--gzvcf TotalRawSNPs.vcf.gz
	--recode-INFO-all
	--minQ 20
	--max-missing 0.2
	--out TRSboth
	--recode

  Outputting VCF file...
  After filtering, kept 41303790 out of a possible 55253967 Sites
  Run Time = 11462.00 seconds
```

**Next step is to mark any genotypes with less than 5 reads as missing.**

```
$ vcftools --vcf TRSboth.recode.vcf --minDP 5 --recode --recode-INFO-all --out TRSBOTHdp5
```

```
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
        --vcf TRSboth.recode.vcf
        --recode-INFO-all
        --minDP 5
        --out TRSBOTHdp5
        --recode

After filtering, kept 90 out of 90 Individuals
Outputting VCF file...
After filtering, kept 41303790 out of a possible 41303790 Sites
Run Time = 12238.00 seconds
```

**Next, remove sites with more than 50% missing data.**

```
$ vcftools --vcf TRSBOTHdp5.recode.vcf --max-missing 0.5 --recode --recode-INFO-all --out TRSBOTHdp5g5
```

```
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
        --vcf TRSBOTHdp5.recode.vcf
        --recode-INFO-all
        --max-missing 0.5
        --out TRSBOTHdp5g5
        --recode

After filtering, kept 90 out of 90 Individuals
Outputting VCF file...
After filtering, kept 2248081 out of a possible 41303790 Sites
Run Time = 1879.00 seconds
```

**MAF filtering for FreeBayes Output.**

```
$ curl -L -O https://raw.githubusercontent.com/jpuritz/dDocent/master/scripts/untested/multi.maf.sh
$ chmod +x multi.maf.sh
$ multi.maf.sh TRSBOTHdp5g5.recode.vcf 0.001 TRSBOTHdp5g5maf
```

```
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf TRSBOTHdp5g5.recode.vcf
	--recode-INFO-all
	--out TRSBOTHdp5g5maf
	--positions maf.loci.to.keep
	--recode

After filtering, kept 90 out of 90 Individuals
Outputting VCF file...
After filtering, kept 2232282 out of a possible 2248081 Sites
Run Time = 619.00 seconds
```

**Use a custom script called `filter_missing_ind.sh` to filter out bad individuals.**

Here I mostly just want to see % of missing data within individuals. We are working with a small number of individuals, so I don't want to remove any.

```
$ curl -L -O https://raw.githubusercontent.com/jpuritz/dDocent/master/scripts/filter_missing_ind.sh
$ chmod +x filter_missing_ind.sh
$ ./filter_missing_ind.sh TRSBOTHdp5g5maf.recode.vcf TRSBOTHdp5g5mafMI
```

```
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf TRSBOTHdp5g5maf.recode.vcf
	--missing-indv
	--out TRSBOTHdp5g5mafMI

After filtering, kept 90 out of 90 Individuals
Outputting Individual Missingness
After filtering, kept 2232282 out of a possible 2232282 Sites
Run Time = 39.00 seconds



                                          Histogram of % missing data per individual
        9 +---------------------------------------------------------------------------------------------------------+
          |          **+             +            +            +            +             +            +            |
          |          **                                       'totalmissing' using (bin($1,binwidth)):(1.0) ******* |
        8 |-+        **   **                                                                                      +-|
          |          **   **                                                                                        |
        7 |-+        **   **                                                                                      +-|
          |          **   **                                                                                        |
          |          **   **                                                                                        |
        6 |-+        **   **                                                                                      +-|
          |          **   **                                                                                        |
        5 |-+        **   **     **                                                                               +-|
          |          **   **     **                                                                                 |
          |          **   **     **                                                                                 |
        4 |-+        ***  **     **                                                                               +-|
          |          ***  **     **                                                                                 |
        3 |-+        ***  **     **        **                  **          **            ***                      +-|
          |          ***  **     **        **                  **          **            * *                        |
          |          ***  **     **        **                  **          **            * *                        |
        2 |-+        *******  ** ***** **  **   ******  *** *******    **  **   ******   * *           ***        +-|
          |          *** ***  ** **  * **  **   *    *  * * * *** *    **  **   *  * *   * *           * *          |
        1 |-+     ****** **********  ************    **** *** *** ***************  * ***** ************* ***      +-|
          |       ** *** *** *** **  * ** ***   *    *  * * * *** * * *** ***   *  * * * * *   *  *  * * * *        |
          |       ** *** *** *** **  * ** ***   * +  *  * * * *** * * *** ***   *  * * * *+*   *  *  * * * *        |
        0 +---------------------------------------------------------------------------------------------------------+
          0           0.1           0.2          0.3          0.4          0.5           0.6          0.7          0.8
                                                       % of missing data

The 85% cutoff would be 0.576223
Would you like to set a different cutoff, yes or no
$ yes #a 85% cutoff would remove 12 individuals
Please enter new cutoff
$ 0.95
All individuals with more than 95.0% missing data will be removed.

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf TRSBOTHdp5g5maf.recode.vcf
	--remove lowDP.indv
	--recode-INFO-all
	--out TRSBOTHdp5g5mafMI
	--recode

  Excluding individuals in 'exclude' list
  After filtering, kept 90 out of 90 Individuals
  Outputting VCF file...
  After filtering, kept 2232282 out of a possible 2232282 Sites
  Run Time = 622.00 seconds
```

Here is a list of the worst 12 samples:

```
$ cat <(head -1 TRSBOTHdp5g5mafMI.imiss ) <(mawk '!/F_MI/' TRSBOTHdp5g5mafMI.imiss | sort -k5 -r -n ) | head -13 | column -t
```

```
INDV    N_DATA   N_GENOTYPES_FILTERED  N_MISS   F_MISS
PVD_4   2232282  0                     1648178  0.738338
PVD_5   2232282  0                     1610297  0.721368
PVD_3   2232282  0                     1593186  0.713703
PVD_1   2232282  0                     1591357  0.712883
PVD_2   2232282  0                     1545093  0.692159
NAR2_1  2232282  0                     1506301  0.674781
MCD_10  2232282  0                     1472556  0.659664
KIC_10  2232282  0                     1386910  0.621297
MCD_1   2232282  0                     1357619  0.608175
KIC_5   2232282  0                     1353790  0.60646
NAR2_4  2232282  0                     1344125  0.60213
KIC_9   2232282  0                     1296270  0.580693
```

Looking back at the multiqc report, these 12 samples had lower numbers of total sequence compared to the other samples. These 12 samples had less than 5 million total sequences.

**Use a second custom `script pop_missing_filter.sh` to filter loci that have high missing data values in a single population.**

This step needs a file that maps individuals to populations `popmap`.

```
$ head popmap
```

```
BAR_10	BAR
BAR_1	BAR
BAR_2	BAR
BAR_3	BAR
BAR_4	BAR
BAR_5	BAR
BAR_6	BAR
BAR_7	BAR
BAR_8	BAR
BAR_9	BAR
```

```
$ curl -L -O https://raw.githubusercontent.com/jpuritz/dDocent/master/scripts/pop_missing_filter.sh
$ chmod +x pop_missing_filter.sh
$ ./pop_missing_filter.sh TRSBOTHdp5g5mafMI.recode.vcf popmap 0.1 1 TRSBOTHdp5g5mafMIp9
```

```
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf TRSBOTHdp5g5mafMI.recode.vcf
	--keep keep.BAR
	--out BAR
	--missing-site

Keeping individuals in 'keep' list
After filtering, kept 10 out of 90 Individuals
Outputting Site Missingness
After filtering, kept 2232282 out of a possible 2232282 Sites
Run Time = 41.00 seconds

# This script will run through each population and print an output like the one above. In all populations, after filtering all sites are retained.

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf TRSBOTHdp5g5mafMI.recode.vcf
	--exclude-positions loci.to.remove
	--recode-INFO-all
	--out TRSBOTHdp5g5mafMIp9
	--recode

After filtering, kept 90 out of 90 Individuals
Outputting VCF file...
After filtering, kept 162972 out of a possible 2232282 Sites
Run Time = 63.00 seconds
```

**Filter out any sites with less than 95% overall call rate and MAF of 0.001.**

```
$ vcftools --vcf TRSBOTHdp5g5mafMIp9.recode.vcf --recode-INFO-all --max-missing 0.95 --maf 0.001 --out TRSBOTHdp5g5mafMIp9g9 --recode
```

```
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf TRSBOTHdp5g5mafMIp9.recode.vcf
	--recode-INFO-all
	--maf 0.001
	--max-missing 0.95
	--out TRSBOTHdp5g5mafMIp9g9
	--recode

After filtering, kept 90 out of 90 Individuals
Outputting VCF file...
After filtering, kept 160137 out of a possible 162972 Sites
Run Time = 48.00 seconds
```

**Next, split the VCF files into nDNA and mtDNA for further filtering**

```
$ cat TRSBOTHdp5g5mafMIp9g9.recode.vcf | mawk '!/NC_007175.2/' > TRSBOTHdp5g5mafMIp9g9.nDNA.vcf

$ cat TRSBOTHdp5g5mafMIp9g9.recode.vcf | head -10000 | mawk '/#/' > header
$ cat TRSBOTHdp5g5mafMIp9g9.recode.vcf | head -10000 | mawk '/#/' > header
$ cat header <(cat TRSBOTHdp5g5mafMIp9g9.recode.vcf |head -20000| mawk '/NC_007175.2/') > TRSBOTHdp5g5mafMIp9g9.mtDNA.vcf
```

The remaining filtering steps will be split between the mtDNA and nDNA

#### mtDNA

**Break complex mutational events (combinations of SNPs and INDELs) into separate SNP and INDEL calls, and then remove INDELs.**

```
$ vcffilter -s -f "AB < 0.001" TRSBOTHdp5g5mafMIp9g9.mtDNA.vcf | vcffilter -s -f "QUAL / DP > 0.25" > TRSBOTHdp5g5mafMIp9g9.mtDNA.F.vcf
$ vcfallelicprimitives -k -g TRSBOTHdp5g5mafMIp9g9.mtDNA.F.vcf | sed 's:\.|\.:\.\/\.:g' > TRSBOTHdp5g5mafMIp9g9mtDNA.F.prim
$ vcftools --vcf TRSBOTHdp5g5mafMIp9g9mtDNA.F.prim --remove-indels --recode --recode-INFO-all --out SNP.TRSBOTHdp5g5mafMIp9g9mtDNA
```

```
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf TRSBOTHdp5g5mafMIp9g9mtDNA.F.prim
	--recode-INFO-all
	--out SNP.TRSBOTHdp5g5mafMIp9g9mtDNA
	--recode
	--remove-indels

After filtering, kept 90 out of 90 Individuals
Outputting VCF file...
After filtering, kept 632 out of a possible 655 Sites
Run Time = 1.00 seconds
```

#### nDNA

**Use a custom filter script `dDocent_ngs_filters`.**

```
$ curl -L -O https://raw.githubusercontent.com/amyzyck/EecSeq_NB_EasternOyster/master/Scripts/dDocent_ngs_filters
$ chmod +x dDocent_ngs_filters
$ ./dDocent_ngs_filters TRSBOTHdp5g5mafMIp9g9.nDNA.vcf TRSBOTHdp5g5mafMIp9g9dnDNA
```

```
This script will automatically filter a FreeBayes generated VCF file using criteria related to site depth,
quality versus depth, allelic balance at heterzygous individuals, and paired read representation.
The script assumes that loci and individuals with low call rates (or depth) have already been removed.

Contact Jon Puritz (jpuritz@gmail.com) for questions and see script comments for more details on particular filters

Is this from a mixture of SE and PE libraries? Enter yes or no.
$ no
Number of sites filtered based on allele balance at heterozygous loci, locus quality, and mapping quality / Depth
 41283 of 159239

Number of additional sites filtered based on properly paired status
 44 of 117956

Number of sites filtered based on high depth and lower than 2*DEPTH quality score
 19169 of 117956




                                               Histogram of mean depth per site
     6000 +---------------------------------------------------------------------------------------------------------+
          |  +     +    +    +    +    +    +    +    +     +    +    +    +    +    +    +    +     +    +    +    |
          |                *** **                         'meandepthpersite' using (bin($1,binwidth)):(1.0) ******* |
          |                * ******                                                                                 |
     5000 |-+            *** * ** ***                                                                             +-|
          |              * * * ** * *                                                                               |
          |             ** * * ** * ****                                                                            |
          |             ** * * ** * ** *                                                                            |
     4000 |-+           ** * * ** * ** ***                                                                        +-|
          |           **** * * ** * ** * **                                                                         |
          |           * ** * * ** * ** * ****                                                                       |
     3000 |-+         * ** * * ** * ** * ** ***                                                                   +-|
          |           * ** * * ** * ** * ** * ***                                                                   |
          |         *** ** * * ** * ** * ** * * **                                                                  |
          |         * * ** * * ** * ** * ** * * ******                                                              |
     2000 |-+       * * ** * * ** * ** * ** * * ** * ****                                                         +-|
          |         * * ** * * ** * ** * ** * * ** * ** *****                                                       |
          |        ** * ** * * ** * ** * ** * * ** * ** * * **                                                      |
          |        ** * ** * * ** * ** * ** * * ** * ** * * **                                                      |
     1000 |-+      ** * ** * * ** * ** * ** * * ** * ** * * **                                                    +-|
          |      **** * ** * * ** * ** * ** * * ** * ** * * **                                                      |
          |      * ** * ** * * ** * ** * ** * * ** * ** * * *******                                                 |
          |  + *** ** * ** * * ** * ** * ** * * ** * ** * * ** * **************************************** ****** ** |
        0 +---------------------------------------------------------------------------------------------------------+
             12    15   18   21   24   27   30   33   36    39   42   45   48   51   54   57   60    63   66   69   72
                                                          Mean Depth

The 95% cutoff would be 58
Would you like to use a different maximum mean depth cutoff than 58, yes or no
$ no
Maximum mean depth cutoff is 58
Number of sites filtered based on maximum mean depth
 5026 of 98751

Total number of sites filtered
 65514 of 159239

Remaining sites
 93725

Filtered VCF file is called Output_prefix.FIL.recode.vcf

Filter stats stored in TRSBOTHdp5g5mafMIp9g9dnDNA.filterstats
```

**Break complex mutational events (combinations of SNPs and INDELs) into separate SNP and INDEL calls, and then remove INDELs.**

```
$ vcfallelicprimitives TRSBOTHdp5g5mafMIp9g9dnDNA.FIL.recode.vcf --keep-info --keep-geno > TRSBOTHdp5g5mafMIp9g9dnDNA.prim.vcf
$ vcftools --vcf TRSBOTHdp5g5mafMIp9g9dnDNA.prim.vcf --remove-indels --recode --recode-INFO-all --out SNP.TRSBOTHdp5g5mafMIp9g9dnDNA
```

```
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf TRSBOTHdp5g5mafMIp9g9dnDNA.prim.vcf
	--recode-INFO-all
	--out SNP.TRSBOTHdp5g5mafMIp9g9dnDNA
	--recode
	--remove-indels

After filtering, kept 90 out of 90 Individuals
Outputting VCF file...
After filtering, kept 121779 out of a possible 139791 Sites
Run Time = 36.00 seconds
```

***
**Filtered by a few different Minor Allele Frequencies**

Minor Allele Frequency greater than 1%

```
$ vcftools --vcf SNP.TRSBOTHdp5g5mafMIp9g9dnDNA.recode.vcf --maf 0.01 --recode --recode-INFO-all --out SNP.TRSBOTHdp5g5mafMIp9g9dnDNAmaf01
```

```
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf SNP.TRSBOTHdp5g5mafMIp9g9dnDNA.recode.vcf
	--recode-INFO-all
	--maf 0.01
	--out SNP.TRSBOTHdp5g5mafMIp9g9dnDNAmaf01
	--recode

After filtering, kept 90 out of 90 Individuals
Outputting VCF file...
After filtering, kept 93934 out of a possible 121779 Sites
Run Time = 27.00 seconds
```

Minor Allele Frequency greater than 2.5%

```
$ vcftools --vcf SNP.TRSBOTHdp5g5mafMIp9g9dnDNA.recode.vcf --maf 0.025 --recode --recode-INFO-all --out SNP.TRSBOTHdp5g5mafMIp9g9dnDNAmaf025
```

```
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf SNP.TRSBOTHdp5g5mafMIp9g9dnDNA.recode.vcf
	--recode-INFO-all
	--maf 0.025
	--out SNP.TRSBOTHdp5g5mafMIp9g9dnDNAmaf025
	--recode

After filtering, kept 90 out of 90 Individuals
Outputting VCF file...
After filtering, kept 62729 out of a possible 121779 Sites
Run Time = 20.00 seconds
```

Minor Allele Frequency greater than 5%

```
$ vcftools --vcf SNP.TRSBOTHdp5g5mafMIp9g9dnDNA.recode.vcf --maf 0.05 --recode --recode-INFO-all --out SNP.TRSBOTHdp5g5mafMIp9g9dnDNAmaf05
```

```
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf SNP.TRSBOTHdp5g5mafMIp9g9dnDNA.recode.vcf
	--recode-INFO-all
	--maf 0.05
	--out SNP.TRSBOTHdp5g5mafMIp9g9dnDNAmaf05
	--recode

After filtering, kept 90 out of 90 Individuals
Outputting VCF file...
After filtering, kept 47550 out of a possible 121779 Sites
Run Time = 15.00 seconds
```

**Max 2 alleles**

```
$ vcftools --vcf SNP.TRSBOTHdp5g5mafMIp9g9dnDNAmaf05.recode.vcf --recode --recode-INFO-all --out SNP.TRSBOTHdp5g5mafMIp9g9dnDNAmaf052A --max-alleles 2
```

```
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf SNP.TRSBOTHdp5g5mafMIp9g9dnDNAmaf05.recode.vcf
	--recode-INFO-all
	--max-alleles 2
	--out SNP.TRSBOTHdp5g5mafMIp9g9dnDNAmaf052A
	--recode

After filtering, kept 90 out of 90 Individuals
Outputting VCF file...
After filtering, kept 46450 out of a possible 47550 Sites
Run Time = 14.00 seconds
```

**Continue on to Outlier Detection using `SNP.TRSBOTHdp5g5mafMIp9g9dnDNAmaf052A.recode.vcf`**

Steps for Outlier Detection can be accessed [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Analysis/Analysis_Part2/EecSeq_Cvirginica_OutlierDetection_Part2.md)
