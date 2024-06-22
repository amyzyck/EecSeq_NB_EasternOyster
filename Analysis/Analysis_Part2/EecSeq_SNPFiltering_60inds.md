# Re-doing SNP Filtering with NAR, NAR2, and GHP individuals removed 
Talked to Jacob and figured it would be good to re-do the SNP filtering steps with NAR, NAR2, and GHP individuals removed at the beginning. Right now, I'm not sure if including all 90 individuals will influence which and how many SNPs are retained during the filtering process. Testing it out now in case there is a significant difference. 

In `NB_SNPFiltering_both` directory, make new directory for SNP Filtering with just 6 populations  

```
$ mkdir NB_SNPFiltering_6pops
```

Copy `TotalRawSNPs.vcf.gz` file from other directory 

```
$ cp ../TotalRawSNPs.vcf.gz . 
```

Make a text file with the individual names for NAR, NAR2, and GHP, it should look like this:

```
$ head NAR_NAR2_GHP

output:
NAR2_1
NAR2_2
NAR2_3
NAR2_4
NAR2_5
NAR2_6
NAR2_7
NAR2_8
NAR2_9
NAR2_10
```

Remove these individuals from the vcf file using vcftools. Usee --gzvcf flag since the VCF file is gzipped. Use --remove followed by the file containing the list of individuals to exclude from the VCF file. 

```
$ vcftools --gzvcf TotalRawSNPs.vcf.gz --remove NAR_NAR2_GHP --recode-INFO-all --out TotalRawSNPs6 --recode

output: 
Excluding individuals in 'exclude' list
After filtering, kept 60 out of 90 Individuals
Outputting VCF file..
After filtering, kept 55253967 out of a possible 55253967 Sites
Run Time = 12645.00 seconds
```

G-zipping file immediately after to reduce computational space 

```
$ gzip TotalRawSNPs6.recode.vcf
```

Now, I'll go through a repeat all the same filtering steps as before, in the same order. I'll document the output of each step here.

**First step for filtering is to remove sites with a QUALITY score below 20 and those sites that have more than 80% missing data.**

```
$ vcftools --gzvcf TotalRawSNPs6.recode.vcf.gz --minQ 20 --max-missing 0.20 --recode --recode-INFO-all --out TRS6
```

```
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--gzvcf TotalRawSNPs6.recode.vcf.gz
	--recode-INFO-all
	--minQ 20
	--max-missing 0.2
	--out TRS6
	--recode

Outputting VCF file...
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 41036259 out of a possible 55253967 Sites                                                                                   
Run Time = 9025.00 seconds
```

**Next step is to mark any genotypes with less than 10 reads as missing.**

```
$ vcftools --vcf TRS6.recode.vcf --minDP 10 --recode --recode-INFO-all --out TRS6dp5

output:
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 41036259 out of a possible 41036259 Sites
Run Time = 8933.00 seconds
```

**Next, remove sites with more than 50% missing data.**

```
$ vcftools --vcf TRS6dp10.recode.vcf --max-missing 0.5 --recode --recode-INFO-all --out TRS6dp10g5

output:
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 506958 out of a possible 41036259 Sites
Run Time = 791.00 seconds
```

**MAF filtering for FreeBayes Output.**

```
$ curl -L -O https://raw.githubusercontent.com/jpuritz/dDocent/master/scripts/untested/multi.maf.sh
$ chmod +x multi.maf.sh
$ multi.maf.sh TRS6dp10g5.recode.vcf 0.001 TRS6dp10g5maf

output:
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 502909 out of a possible 506958 Sites
Run Time = 104.00 seconds
```

**Use a custom script called `filter_missing_ind.sh` to filter out bad individuals.**

Here I mostly just want to see % of missing data within individuals. We are working with a small number of individuals, so I don't want to remove any.

```
$ curl -L -O https://raw.githubusercontent.com/jpuritz/dDocent/master/scripts/filter_missing_ind.sh
$ chmod +x filter_missing_ind.sh
$ ./filter_missing_ind.sh TRSBOTHdp5g5maf.recode.vcf TRSBOTHdp5g5mafMI

output: 
After filtering, kept 60 out of 60 Individuals
Outputting Individual Missingness
After filtering, kept 502909 out of a possible 502909 Sites
Run Time = 11.00 seconds


                                          Histogram of % missing data per individual
      Number of Occurrences
        6 ++--**-------+-------------+------------+------------+------------+-------------+------------+-----------++
          +   **       +             +            +            'totalmissing' using (bin($1,binwidth)):(1.0) ****** +
          |   **                                                                                                    |
          |   **                                                                                                    |
          |   **                                                                                                    |
        5 ++  ** **                                                                                                ++
          |   ** **                                                                                                 |
          |   ** **                                                                                                 |
          |   ** **                                                                                                 |
        4 ++ *** **                                                                                                ++
          |  *** **                                                                                                 |
          |  *** **                                                                                                 |
          |  *** **                                                                                                 |
          |  *** **                                                                                                 |
        3 ++ *********                                ***                       ******                             ++
          |  *** **  *                                * *                       *    *                              |
          |  *** **  *                                * *                       *    *                              |
          |  *** **  *                                * *                       *    *                              |
        2 ++ *** **  * ***       ********         *** * *          ***          *    *                             ++
          |  *** **  * * *       *      *         * * * *          * *          *    *                              |
          |  *** **  * * *       *      *         * * * *          * *          *    *                              |
          |  *** **  * * *       *      *         * * * *          * *          *    *                              |
          +  *** **  * * *       *   +  *         * * * *      +   * *      +   *    *    +            +            +
        1 ++-*********-***-------********---------***-***------+---***------+---******----+------------+-----------++
          0           0.1           0.2          0.3          0.4          0.5           0.6          0.7          0.8
                                                       % of missing data

The 85% cutoff would be 0.584633
Would you like to set a different cutoff, yes or no
yes
Please enter new cutoff
0.95
All individuals with more than 95.0% missing data will be removed.

After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 502909 out of a possible 502909 Sites
Run Time = 106.00 seconds
```

Here is a list of the worst 12 samples

```
$ cat <(head -1 TRS6dp10g5mafMI.imiss ) <(mawk '!/F_MI/' TRS6dp10g5mafMI.imiss | sort -k5 -r -n ) | head -13 | column -t

output: 
INDV    N_DATA  N_GENOTYPES_FILTERED  N_MISS  F_MISS
PVD_4   502909  0                     357192  0.710252
PVD_1   502909  0                     338131  0.67235
PVD_5   502909  0                     328428  0.653057
KIC_10  502909  0                     324692  0.645628
PVD_3   502909  0                     318969  0.634248
MCD_10  502909  0                     308983  0.614391
PVD_2   502909  0                     306765  0.609981
KIC_9   502909  0                     300543  0.597609
BAR_10  502909  0                     294017  0.584633
KIC_5   502909  0                     289109  0.574873
BAR_4   502909  0                     278569  0.553915
MCD_1   502909  0                     278364  0.553508
```

Looking back at the multiqc report, these 12 samples had lower numbers of total sequence compared to the other samples. These 12 samples had less than 5 million total sequences. Again, I'm not going to remove individuals. 

**Use a second custom `script pop_missing_filter.sh` to filter loci that have high missing data values in a single population.**

This step needs a file that maps individuals to populations `popmap`.

```
$ head popmap

output: 
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

This will run through each population before spitting out a final output:

```
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 61556 out of a possible 502909 Sites
Run Time = 17.00 seconds
```

**Filter out any sites with less than 95% overall call rate and MAF of 0.001.**

```
$ vcftools --vcf TRS6dp10g5mafMIp9.recode.vcf --recode-INFO-all --max-missing 0.95 --maf 0.001 --out TRS6dp10g5mafMIp9g95 --recode

output:
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 50863 out of a possible 61556 Sites
Run Time = 11.00 seconds
```

**Next, split the VCF files into nDNA and mtDNA for further filtering**

```
$ cat TRS6dp10g5mafMIp9g95.recode.vcf | mawk '!/NC_007175.2/' > TRS6dp10g5mafMIp9g95.nDNA.vcf

$ cat TRS6dp10g5mafMIp9g95.recode.vcf | head -10000 | mawk '/#/' > header
$ cat TRS6dp10g5mafMIp9g95.recode.vcf | head -10000 | mawk '/#/' > header
$ cat header <(cat TRS6dp10g5mafMIp9g95.recode.vcf |head -20000| mawk '/NC_007175.2/') > TRS6dp10g5mafMIp9g95.mtDNA.vcf
```

The remaining filtering steps will be split between the mtDNA and nDNA

#### mtDNA

**Break complex mutational events (combinations of SNPs and INDELs) into separate SNP and INDEL calls, and then remove INDELs.**

```
$ vcffilter -s -f "AB < 0.001" TRS6dp10g5mafMIp9g95.mtDNA.vcf | vcffilter -s -f "QUAL / DP > 0.25" > TRS6dp10g5mafMIp9g95.mtDNA.F.vcf
$ vcfallelicprimitives -k TRS6dp10g5mafMIp9g95.mtDNA.F.vcf | sed 's:\.|\.:\.\/\.:g' > TRS6dp10g5mafMIp9g95mtDNA.F.prim
$ vcftools --vcf TRS6dp10g5mafMIp9g95mtDNA.F.prim --remove-indels --recode --recode-INFO-all --out SNP.TRS6dp10g5mafMIp9g95mtDNA

output:
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 462 out of a possible 473 Sites
Run Time = 0.00 seconds
```

#### nDNA

**Use a custom filter script `dDocent_ngs_filters`.**

```
$ curl -L -O https://raw.githubusercontent.com/amyzyck/EecSeq_NB_EasternOyster/master/Scripts/dDocent_ngs_filters
$ chmod +x dDocent_ngs_filters
$ ./dDocent_ngs_filters TRS6dp10g5mafMIp9g95.nDNA.vcf TRS6dp10g5mafMIp9g95nDNA

output:
This script will automatically filter a FreeBayes generated VCF file using criteria related to site depth,
quality versus depth, allelic balance at heterzygous individuals, and paired read representation.
The script assumes that loci and individuals with low call rates (or depth) have already been removed.

Contact Jon Puritz (jpuritz@gmail.com) for questions and see script comments for more details on particular filters

Is this from a mixture of SE and PE libraries? Enter yes or no.
no
Number of sites filtered based on allele balance at heterozygous loci, locus quality, and mapping quality / Depth
 12989 of 50165 

Number of additional sites filtered based on properly paired status
 11 of 37176 

Number of sites filtered based on high depth and lower than 2*DEPTH quality score
 6259 of 37176 



                                             Histogram of mean depth per site

  1400 ++----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----++
       +     +     +     +     +     +     +     +     +   'meandepthpersite' using (bin($1,binwidth)):(1.0) ******+|
       |                                                                                                            |
  1200 ++                           ***                                                                            ++
       |                           *******                                                                          |
       |                         ******* **                                                                         |
       |                        ** ***** ****                                                                       |
  1000 ++                      *** ***** *****                                                                     ++
       |                       *** ***** *****                                                                      |
       |                      **** ***** ***********                                                                |
   800 ++                     **** ***** ***** *****                                                               ++
       |                     ***** ***** ***** ********                                                             |
       |                     ***** ***** ***** ***** **                                                             |
   600 ++                    ***** ***** ***** ***** *****                                                         ++
       |                   ******* ***** ***** ***** *******                                                        |
       |                  ** ***** ***** ***** ***** ***** ****                                                     |
   400 ++                 ** ***** ***** ***** ***** ***** ******                                                  ++
       |                  ** ***** ***** ***** ***** ***** ********                                                 |
       |                 *** ***** ***** ***** ***** ***** ****** *                                                 |
       |                 *** ***** ***** ***** ***** ***** ****** **                                                |
   200 ++               **** ***** ***** ***** ***** ***** ****** ****                                             ++
       |               ***** ***** ***** ***** ***** ***** ****** *************                                     |
       +     +     + ******* ***** ***** ***** ***** ***** ****** ***** ***** **********************************  **|
     0 ++----+----***************************************************************************************************
       10    15    20    25    30    35    40    45    50    55    60    65    70    75    80    85    90    95   100
                                                        Mean Depth

The 95% cutoff would be 81
Would you like to use a different maximum mean depth cutoff than 81, yes or no
no
Maximum mean depth cutoff is 81
Number of sites filtered based on maximum mean depth
 1558 of 30909 

Total number of sites filtered
 20814 of 50165 

Remaining sites
 29351 

Filtered VCF file is called Output_prefix.FIL.recode.vcf

Filter stats stored in TRS6dp10g5mafMIp9g95nDNA.filterstats
```

**Break complex mutational events (combinations of SNPs and INDELs) into separate SNP and INDEL calls, and then remove INDELs.**

```
$ vcfallelicprimitives TRS6dp10g5mafMIp9g95nDNA.FIL.recode.vcf --keep-info --keep-geno > TRS6dp10g5mafMIp9g95nDNA.prim.vcf
$ vcftools --vcf TRS6dp10g5mafMIp9g95nDNA.prim.vcf --remove-indels --recode --recode-INFO-all --out SNP.TRS6dp10g5mafMIp9g95nDNA

output:
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 22588 out of a possible 25918 Sites
Run Time = 6.00 seconds
```

***
**Filtered by a few different Minor Allele Frequencies**

Minor Allele Frequency greater than 1%

```
$ vcftools --vcf SNP.TRS6dp10g5mafMIp9g95nDNA.recode.vcf --maf 0.01 --recode --recode-INFO-all --out SNP.TRS6dp10g5mafMIp9g95nDNAmaf01

output:
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 15894 out of a possible 22588 Sites
Run Time = 4.00 seconds
```

Minor Allele Frequency greater than 2.5%

```
$ vcftools --vcf SNP.TRS6dp10g5mafMIp9g95nDNA.recode.vcf --maf 0.025 --recode --recode-INFO-all --out SNP.TRS6dp10g5mafMIp9g95nDNAmaf025

output:
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 12730 out of a possible 22588 Sites
Run Time = 6.00 seconds
```

Minor Allele Frequency greater than 5%

```
$ vcftools --vcf SNP.TRS6dp10g5mafMIp9g95nDNA.recode.vcf --maf 0.05 --recode --recode-INFO-all --out SNP.TRS6dp10g5mafMIp9g95nDNAmaf05

output:
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 8700 out of a possible 22588 Sites
Run Time = 6.00 seconds
```

I'm going to move on to the last filtering step with both the 2.5% maf and 5% maf files. 

**Max 2 alleles**

2.5% maf file

```
$ vcftools --vcf SNP.TRS6dp10g5mafMIp9g95nDNAmaf025.recode.vcf --recode --recode-INFO-all --out SNP.TRS6dp10g5mafMIp9g95nDNAmaf0252A --max-alleles 2

output:
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 12456 out of a possible 12730 Sites
Run Time = 7.00 seconds
```

5% maf file

```
$ vcftools --vcf SNP.TRS6dp10g5mafMIp9g95nDNAmaf05.recode.vcf --recode --recode-INFO-all --out SNP.TRS6dp10g5mafMIp9g95nDNAmaf052A --max-alleles 2

output:
After filtering, kept 60 out of 60 Individuals
Outputting VCF file...
After filtering, kept 8586 out of a possible 8700 Sites
Run Time = 5.00 seconds
```

This is a large difference in final SNP count between the 60 individuals and 90 individuals. It seems like the step where sites with more than 50% missing data are removed is where this major divergence starts. 

