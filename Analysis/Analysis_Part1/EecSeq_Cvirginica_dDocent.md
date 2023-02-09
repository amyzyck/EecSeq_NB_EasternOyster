# *Crassostrea virginica* EecSeq Seascape Genomics Analysis - dDocent pipeline

Author: Amy Zyck

Last Edited: March 2, 2020

Data uploaded and analyzed on KITT. User logged in before following steps are completed.

**Data location: `PATH /home/azyck/NB_capture`**

## In Terminal:

### 1. Setup: Downloading software programs, creating environments anaconda folders

**Downloading Bioconda (skip this step if Bioconda has already been downloaded).**

```javascript
# downloading Miniconda software
$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ chmod +x Miniconda3-latest-Linux-x86_64.sh
$ ./Miniconda3-latest-Linux-x86_64.sh

# restarting window with source command
$ source ~/.bashrc

# adding different channels
$ conda config --add channels defaults
	# should get: "warning: 'defaults' already in 'channels' list, moving to the top"
$ conda config --add channels bioconda
$ conda config --add channels conda-forge

# to see if this worked with cat command
$ cat .condarc
```
> Bioconda is a bioinformatics software package manager. See more at [https://bioconda.github.io](https:///bioconda.github.io).

**Create and activate a dDocent conda environment:**

```
$ conda create -n nb_capture ddocent
$ conda activate nb_capture

# the beginning of the line should then look like: (nb_capture) [username@KITT ~]$
```

**Creating and entering a directory for this project**

```
$ mkdir NB_capture # creates a directory
$ cd NB_capture # navigates into this new directory

# the beginning of the line should now look like: (nb_capture) [username@KITT NB_capture]$
```



**Link raw data to `NB_Capture` directory**

The raw data was uploaded to KITT by [J. Puritz](https://github.com/jpuritz).

```
$ ln -s /RAID_STORAGE2/azyck/NB_capture /home/azyck

# NP_capture directory (containing all files) within RAID_STORAGE2 will be linked to the azyck home directory in KITT
```

### 2. FastQC Analysis

```
$ mkdir multiplexed_fastqc_results
$ cd multiplexed_fastqc_results
$ conda install -c bioconda fastqc
$ fastqc ../*fastq.gz
$ mv *fastqc.* multiplexed_fastqc_results/ # moves fastqc files to fastqc directory
```
> [FastQC](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf) is a program designed to visualize the quality of high throughput sequencing datasets. The report will highlight any areas where the data looks unsual or unreliable.

**Multiqc Analysis**

```
$ conda install -c bioconda multiqc
$ multiqc .

$ mv multiqc_report.html multiplexed_multiqc_report.html # rename file

# copy multiqc report to directory outside of KITT
# run the following commands outside of KITT
$ scp -r -P XXXX x@kitt.uri.edu:/home/azyck/NB_capture/NB_multiplexed/multiplexed_fastqc_results/multiplexed_multiqc_report.html /Users/azyck/Documents/Research/NB_AdultOyster/Analysis_Outputs/fastqc
# XXXX above indicates the password specific for our lab's KITT.
# scp = secure copy, -r = indicates copying a folder/repository

# open multiqc_report to view the quality of raw, multiplexed data files
```
> Multiqc creates a single report from the FastQC results.

**Sequence Quality**

![SeqQual](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Output/multiqc/multiplexed_fastqc_per_base_sequence_quality.png)

This graph shows the mean quality at each position across the read. The quality looks pretty good. Quality is lower at the front and end of the sequence, but still within the green.

**Per Sequence Quality Scores**

![PerSeq](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Output/multiqc/multiplexed_fastqc_per_sequence_quality_scores.png)

This graph looks at the average quality scores per sequence. Again looks pretty good. Any low quality scores will be trimmed.

**Per Sequence GC Content**

![GC](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Output/multiqc/multiplexed_fastqc_per_sequence_gc_content.png)

If all individuals are from the same species, then this graph will show that the sequences are roughly normally distributed. The slight difference in peaks is likely a result of adapter content differences.

Full multiqc report can be viewed [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Output/multiqc/multiplexed_multiqc_report.html).

### 3. Demultiplexing Data Files

Each library contains EecSeq reads for 5 individuals. From raw data files, R1 indicates forward and R2 indicates reverse sequences.

> Demultiplexing refers to the step in processing where you use barcode information to identify which sequences came from what sample after multiple samples are sequenced together. Barcodes refer to the unqiue sequence that were added to each sample before all the samples are mixed together ([Happy Belly Bioinformatics](https://astrobiomike.github.io/amplicon/demultiplexing)).

#### Using the program [STACKS](http://catchenlab.life.illinois.edu/stacks/):
The following commands are already installed on the system. See the above STACKS link for more information.
Each library needs a separate barcodes file in txt format.

```
# example of barcode file. The name and barcode need to be separated by a single tab.
# The programs used in this pipeline (dDocent) need the name format to be Population_SampleID.

PVD_1	ATCACG
PVD_2	CGATGT
PVD_3	TTAGGC
PVD_4	TGGCCA
PVD_5	ACAGTG
```

**Turning barcode files into a list of barcodes**

```
# make sure to be in the directory where the raw data is located.
# the -f2 command selects the second column of the file

$ cut -f2 LibraryA_barcodes.txt > barcodes_A
$ cut -f2 LibraryB_barcodes.txt > barcodes_B
$ cut -f2 LibraryC_barcodes.txt > barcodes_C
$ cut -f2 LibraryD_barcodes.txt > barcodes_D
$ cut -f2 LibraryE_barcodes.txt > barcodes_E
$ cut -f2 LibraryF_barcodes.txt > barcodes_F
$ cut -f2 LibraryG_barcodes.txt > barcodes_G
$ cut -f2 LibraryH_barcodes.txt > barcodes_H
$ cut -f2 LibraryI_barcodes.txt > barcodes_I
$ cut -f2 LibraryJ_barcodes.txt > barcodes_J

# To view the new barcodes file
$ head barcodes_A
```

**Make a new directory for each library. The output files will be saved in the appropriate directory.**

```
$ mkdir LibraryA
$ mkdir LibraryB
$ mkdir LibraryC
$ mkdir LibraryD
$ mkdir LibraryE
$ mkdir LibraryF
$ mkdir LibraryG
$ mkdir LibraryH
$ mkdir LibraryI
$ mkdir LibraryJ
```
**In each library directory download renaming script from J. Puritz**

```
$ curl -L -O https://raw.githubusercontent.com/jpuritz/dDocent/master/Rename_for_dDocent.sh
$ chmod +x Rename_for_dDocent.sh
```

**Move each barcode txt file to the appropriate library directory.**

```
$ mv LibraryA_barcodes.txt LibraryA/
$ mv LibraryB_barcodes.txt LibraryB/
$ mv LibraryC_barcodes.txt LibraryC/
$ mv LibraryD_barcodes.txt LibraryD/
$ mv LibraryE_barcodes.txt LibraryE/
$ mv LibraryF_barcodes.txt LibraryF/
$ mv LibraryG_barcodes.txt LibraryG/
$ mv LibraryH_barcodes.txt LibraryH/
$ mv LibraryI_barcodes.txt LibraryI/
$ mv LibraryJ_barcodes.txt LibraryJ/
```

**Using the function `process_shortreads` to demultiplex each file.**

Program options can be found [here](http://catchenlab.life.illinois.edu/stacks/comp/process_shortreads.php).

#### Library A
```
$ process_shortreads -1 D501_R1_001.fastq.gz -2 D501_R2_001.fastq.gz -b barcodes_A -r -i gzfastq -o /home/azyck/NB_capture/LibraryA/

output:
Using Phred+33 encoding for quality scores.
Reads trimmed shorter than 31 nucleotides will be discarded.
Found 1 paired input file(s).
Searching for single-end, inlined barcodes.
Loaded 5 barcodes (6bp).
Will attempt to recover barcodes with at most 1 mismatches.
Processing file 1 of 1 [D501_R1_001.fastq.gz]
  Reading data from:
  D501_R1_001.fastq.gz and
  D501_R2_001.fastq.gz
  25443144 total reads; -3206062 ambiguous barcodes; +411476 recovered; -0 low quality reads; 22237082 retained reads.
    0 trimmed reads; 0 orphaned paired-ends.
Closing files, flushing buffers...
Outputing details to log: '/home/azyck/NB_capture/LibraryA/process_shortreads.log'

25443144 total sequences;
  3206062 ambiguous barcode drops;
  0 low quality read drops;
  0 trimmed reads;
  0 orphaned paired-end reads;
22237082 retained reads.
```

In `LibraryA` directory

```
$ cd LibraryA

# rename for dDocent
$ bash Rename_for_dDocent.sh LibraryA_barcodes.txt
```
**Create directory for files that failed demultiplexing.**

This has to be done after each step because the same barcode ID numbers were used in each library. If not done after each library, the previous files will be overrided
```
$ mkdir LibA_failed_demultiplexed
$ mv sample* LibA_failed_demultiplexed
```

#### Library B
```
$ process_shortreads -1 D502_R1_001.fastq.gz -2 D502_R2_001.fastq.gz -b barcodes_B -r -i gzfastq -o /home/azyck/NB_capture/LibraryB/

output:
Using Phred+33 encoding for quality scores.
Reads trimmed shorter than 31 nucleotides will be discarded.
Found 1 paired input file(s).
Searching for single-end, inlined barcodes.
Loaded 5 barcodes (6bp).
Will attempt to recover barcodes with at most 1 mismatches.
Processing file 1 of 1 [D502_R1_001.fastq.gz]
  Reading data from:
  D502_R1_001.fastq.gz and
  D502_R2_001.fastq.gz
  125379006 total reads; -9159574 ambiguous barcodes; +1997470 recovered; -0 low quality reads; 116219432 retained reads.
    0 trimmed reads; 0 orphaned paired-ends.
Closing files, flushing buffers...
Outputing details to log: '/home/azyck/NB_capture/LibraryB/process_shortreads.log'

125379006 total sequences;
  9159574 ambiguous barcode drops;
  0 low quality read drops;
  0 trimmed reads;
  0 orphaned paired-end reads;
116219432 retained reads.
```

In `LibraryB` Directory

```
$ cd LibraryB

# rename for dDocent
$ bash Rename_for_dDocent.sh LibraryB_barcodes.txt
```

Create directory for files that failed demultiplexing

```
$ mkdir LibB_failed_demultiplexed
$ mv sample* LibB_failed_demultiplexed
```

#### Library C
```
$ process_shortreads -1 D503_R1_001.fastq.gz -2 D503_R2_001.fastq.gz -b barcodes_C -r -i gzfastq -o /home/azyck/NB_capture/LibraryC/

output:
Using Phred+33 encoding for quality scores.
Reads trimmed shorter than 31 nucleotides will be discarded.
Found 1 paired input file(s).
Searching for single-end, inlined barcodes.
Loaded 5 barcodes (6bp).
Will attempt to recover barcodes with at most 1 mismatches.
Processing file 1 of 1 [D503_R1_001.fastq.gz]
  Reading data from:
  D503_R1_001.fastq.gz and
  D503_R2_001.fastq.gz
  116892130 total reads; -8392236 ambiguous barcodes; +1812170 recovered; -0 low quality reads; 108499894 retained reads.
    0 trimmed reads; 0 orphaned paired-ends.
Closing files, flushing buffers...
Outputing details to log: '/home/azyck/NB_capture/LibraryC/process_shortreads.log'

116892130 total sequences;
  8392236 ambiguous barcode drops;
  0 low quality read drops;
  0 trimmed reads;
  0 orphaned paired-end reads;
108499894 retained reads.
```

In `LibraryC` Directory

```
$ cd LibraryC

# rename for dDocent
$ bash Rename_for_dDocent.sh LibraryC_barcodes.txt
```

Create directory for files that failed demultiplexing

```
$ mkdir LibC_failed_demultiplexed
$ mv sample* LibC_failed_demultiplexed
```

#### Library D
```
$ process_shortreads -1 D504_R1_001.fastq.gz -2 D504_R2_001.fastq.gz -b barcodes_D -r -i gzfastq -o /home/azyck/NB_capture/LibraryD/

output:
Using Phred+33 encoding for quality scores.
Reads trimmed shorter than 31 nucleotides will be discarded.
Found 1 paired input file(s).
Searching for single-end, inlined barcodes.
Loaded 5 barcodes (6bp).
Will attempt to recover barcodes with at most 1 mismatches.
Processing file 1 of 1 [D504_R1_001.fastq.gz]
  Reading data from:
  D504_R1_001.fastq.gz and
  D504_R2_001.fastq.gz
  110062034 total reads; -6772878 ambiguous barcodes; +1733904 recovered; -0 low quality reads; 103289156 retained reads.
    0 trimmed reads; 0 orphaned paired-ends.
Closing files, flushing buffers...
Outputing details to log: '/home/azyck/NB_capture/LibraryD/process_shortreads.log'

110062034 total sequences;
  6772878 ambiguous barcode drops;
  0 low quality read drops;
  0 trimmed reads;
  0 orphaned paired-end reads;
103289156 retained reads.
```

In `LibraryD` Directory

```
$ cd LibraryD

# rename for dDocent
$ bash Rename_for_dDocent.sh LibraryD_barcodes.txt
```

Create directory for files that failed demultiplexing

```
$ mkdir LibD_failed_demultiplexed
$ mv sample* LibD_failed_demultiplexed
```

#### Library E
```
$ process_shortreads -1 D505_R1_001.fastq.gz -2 D505_R2_001.fastq.gz -b barcodes_E -r -i gzfastq -o /home/azyck/NB_capture/LibraryE/

output:
Using Phred+33 encoding for quality scores.
Reads trimmed shorter than 31 nucleotides will be discarded.
Found 1 paired input file(s).
Searching for single-end, inlined barcodes.
Loaded 5 barcodes (6bp).
Will attempt to recover barcodes with at most 1 mismatches.
Processing file 1 of 1 [D505_R1_001.fastq.gz]
  Reading data from:
  D505_R1_001.fastq.gz and
  D505_R2_001.fastq.gz
  96177288 total reads; -4504434 ambiguous barcodes; +1528266 recovered; -0 low quality reads; 91672854 retained reads.
    0 trimmed reads; 0 orphaned paired-ends.
Closing files, flushing buffers...
Outputing details to log: '/home/azyck/NB_capture/LibraryE/process_shortreads.log'

96177288 total sequences;
  4504434 ambiguous barcode drops;
  0 low quality read drops;
  0 trimmed reads;
  0 orphaned paired-end reads;
91672854 retained reads.
```

In `LibraryE` Directory

```
$ cd LibraryE

# rename for dDocent
$ bash Rename_for_dDocent.sh LibraryE_barcodes.txt
```

Create directory for files that failed demultiplexing

```
$ mkdir LibE_failed_demultiplexed
$ mv sample* LibE_failed_demultiplexed
```

#### Library F
```
$ process_shortreads -1 D506_R1_001.fastq.gz -2 D506_R2_001.fastq.gz -b barcodes_F -r -i gzfastq -o /home/azyck/NB_capture/LibraryF/

output:
Using Phred+33 encoding for quality scores.
Reads trimmed shorter than 31 nucleotides will be discarded.
Found 1 paired input file(s).
Searching for single-end, inlined barcodes.
Loaded 5 barcodes (6bp).
Will attempt to recover barcodes with at most 1 mismatches.
Processing file 1 of 1 [D506_R1_001.fastq.gz]
  Reading data from:
  D506_R1_001.fastq.gz and
  D506_R2_001.fastq.gz
  94254412 total reads; -14785760 ambiguous barcodes; +1566730 recovered; -0 low quality reads; 79468652 retained reads.
    0 trimmed reads; 0 orphaned paired-ends.
Closing files, flushing buffers...
Outputing details to log: '/home/azyck/NB_capture/LibraryF/process_shortreads.log'

94254412 total sequences;
  14785760 ambiguous barcode drops;
  0 low quality read drops;
  0 trimmed reads;
  0 orphaned paired-end reads;
79468652 retained reads.
```

In `LibraryF` Directory

```
$ cd LibraryF

# rename for dDocent
$ bash Rename_for_dDocent.sh LibraryF_barcodes.txt
```

Create directory for files that failed demultiplexing

```
$ mkdir LibF_failed_demultiplexed
$ mv sample* LibF_failed_demultiplexed
```

#### Library G
```
$ process_shortreads -1 D507_R1_001.fastq.gz -2 D507_R2_001.fastq.gz -b barcodes_G -r -i gzfastq -o /home/azyck/NB_capture/LibraryG/

output:
Using Phred+33 encoding for quality scores.
Reads trimmed shorter than 31 nucleotides will be discarded.
Found 1 paired input file(s).
Searching for single-end, inlined barcodes.
Loaded 5 barcodes (6bp).
Will attempt to recover barcodes with at most 1 mismatches.
Processing file 1 of 1 [D507_R1_001.fastq.gz]
  Reading data from:
  D507_R1_001.fastq.gz and
  D507_R2_001.fastq.gz
  89097356 total reads; -1522122 ambiguous barcodes; +1375376 recovered; -0 low quality reads; 87575234 retained reads.
    0 trimmed reads; 0 orphaned paired-ends.
Closing files, flushing buffers...
Outputing details to log: '/home/azyck/NB_capture/LibraryG/process_shortreads.log'

89097356 total sequences;
  1522122 ambiguous barcode drops;
  0 low quality read drops;
  0 trimmed reads;
  0 orphaned paired-end reads;
87575234 retained reads.
```

In `LibraryG` Directory

```
$ cd LibraryG

# rename for dDocent
$ bash Rename_for_dDocent.sh LibraryG_barcodes.txt
```

Create directory for files that failed demultiplexing

```
$ mkdir LibG_failed_demultiplexed
$ mv sample* LibG_failed_demultiplexed
```

#### Library H
```
$ process_shortreads -1 D508_R1_001.fastq.gz -2 D508_R2_001.fastq.gz -b barcodes_H -r -i gzfastq -o /home/azyck/NB_capture/LibraryH/

output:
Using Phred+33 encoding for quality scores.
Reads trimmed shorter than 31 nucleotides will be discarded.
Found 1 paired input file(s).
Searching for single-end, inlined barcodes.
Loaded 5 barcodes (6bp).
Will attempt to recover barcodes with at most 1 mismatches.
Processing file 1 of 1 [D508_R1_001.fastq.gz]
  Reading data from:
  D508_R1_001.fastq.gz and
  D508_R2_001.fastq.gz
  102581558 total reads; -6623442 ambiguous barcodes; +1621604 recovered; -0 low quality reads; 95958116 retained reads.
    0 trimmed reads; 0 orphaned paired-ends.
Closing files, flushing buffers...
Outputing details to log: '/home/azyck/NB_capture/LibraryH/process_shortreads.log'

102581558 total sequences;
  6623442 ambiguous barcode drops;
  0 low quality read drops;
  0 trimmed reads;
  0 orphaned paired-end reads;
95958116 retained reads.
```

In `LibraryH` Directory

```
$ cd LibraryH

# rename for dDocent
$ bash Rename_for_dDocent.sh LibraryH_barcodes.txt
```

Create directory for files that failed demultiplexing

```
$ mkdir LibH_failed_demultiplexed
$ mv sample* LibH_failed_demultiplexed
```

#### Library I
```
$ process_shortreads -1 D509_R1_001.fastq.gz -2 D509_R2_001.fastq.gz -b barcodes_I -r -i gzfastq -o /home/azyck/NB_capture/LibraryI/

output:
Using Phred+33 encoding for quality scores.
Reads trimmed shorter than 31 nucleotides will be discarded.
Found 1 paired input file(s).
Searching for single-end, inlined barcodes.
Loaded 5 barcodes (6bp).
Will attempt to recover barcodes with at most 1 mismatches.
Processing file 1 of 1 [D509_R1_001.fastq.gz]
  Reading data from:
  D509_R1_001.fastq.gz and
  D509_R2_001.fastq.gz
  73432602 total reads; -100212 ambiguous barcodes; +1120346 recovered; -0 low quality reads; 73332390 retained reads.
    0 trimmed reads; 0 orphaned paired-ends.
Closing files, flushing buffers...
Outputing details to log: '/home/azyck/NB_capture/LibraryI/process_shortreads.log'

73432602 total sequences;
  100212 ambiguous barcode drops;
  0 low quality read drops;
  0 trimmed reads;
  0 orphaned paired-end reads;
73332390 retained reads.
```

In `LibraryI` Directory

```
$ cd LibraryI

# rename for dDocent
$ bash Rename_for_dDocent.sh LibraryI_barcodes.txt
```

Create directory for files that failed demultiplexing

```
$ mkdir LibI_failed_demultiplexed
$ mv sample* LibI_failed_demultiplexed
```

#### Library J
```
$ process_shortreads -1 D510_R1_001.fastq.gz -2 D510_R2_001.fastq.gz -b barcodes_J -r -i gzfastq -o /home/azyck/NB_capture/LibraryJ/

output:
Using Phred+33 encoding for quality scores.
Reads trimmed shorter than 31 nucleotides will be discarded.
Found 1 paired input file(s).
Searching for single-end, inlined barcodes.
Loaded 5 barcodes (6bp).
Will attempt to recover barcodes with at most 1 mismatches.
Processing file 1 of 1 [D510_R1_001.fastq.gz]
  Reading data from:
  D510_R1_001.fastq.gz and
  D510_R2_001.fastq.gz
  78638348 total reads; -145222 ambiguous barcodes; +1197050 recovered; -0 low quality reads; 78493126 retained reads.
    0 trimmed reads; 0 orphaned paired-ends.
Closing files, flushing buffers...
Outputing details to log: '/home/azyck/NB_capture/LibraryJ/process_shortreads.log'

78638348 total sequences;
  145222 ambiguous barcode drops;
  0 low quality read drops;
  0 trimmed reads;
  0 orphaned paired-end reads;
78493126 retained reads.
```

In `LibraryJ` Directory

```
$ cd LibraryJ

# rename for dDocent
$ bash Rename_for_dDocent.sh LibraryJ_barcodes.txt
```

Create directory for files that failed demultiplexing

```
$ mkdir LibJ_failed_demultiplexed
$ mv sample* LibJ_failed_demultiplexed
```

`-r` tells radtags to fix cut sites/barcodes that have up to 1-2 mutations in them.

`-i gzfastq` = format of the sequence

`-o` specifies path to output the processed files.

**Copy sample files from each `Library` Directory to `NB_capture`.**

```
# In NB_capture Directory
$ cp LibraryA/*.fq.gz . #copying files ending in .fq.gz from LibraryA directory to directory you are currently in

# Repeat for files in directories LibraryB - LibraryJ
```

There should be 100 files total: 10 from each of 5 populations with a forward (F) and reverse (R) sequence.

**Count raw reads for each file**

```
$ for fq in *.fq.gz
> do
> echo $fq
> zcat $fq | echo $((`wc -l`/4))
> done
```

> Each read in a fastq file has 4 lines: 1. identifier 2. sequence 3. description that starts with "+" 4. quality for each base in the second line. Raw reads can be counted by calculating the number of lines in the file and dividing by 4. `echo`	is used to start an argument and writes the argument as a standard output. `zcat` is used for zipped files.

Raw read counts [output](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Output/raw_reads/RawReads_Counts.txt)

**Checking the quality of data post-demultiplexing**

```
$ mkdir demultiplexed_fastqc_results
$ cd demultiplexed_fastqc_results
$ fastqc ../*fq.gz
$ mv *fastqc.* demultiplexed_fastqc_results
```

**Multiqc Analysis**

```
$ multiqc .

$ mv multiqc_report.html demultiplexed_multiqc_report.html

# copy multiqc report to directory outside of KITT
# run the following commands outside of KITT
$ scp -r -P XXXX x@kitt.uri.edu:/home/azyck/NB_capture/demultiplexed_fastqc_results/demultiplexed_multiqc_report.html /Users/azyck/Documents/Research/NB_AdultOyster/Analysis_Outputs/fastqc
# XXXX above indicates the password specific for our lab's KITT.
# scp = secure copy, -r = indicates copying a folder/repository

# open multiqc_report to view the quality of raw, multiplexed data files
```

Multiqc report can be viewed [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Output/multiqc/demultiplexed_multiqc_report.html).

### 4. Read Trimming, Mapping, and SNP Calling

Using [**dDocent**](http://www.ddocent.com/)

> dDocent is a bioinformatics program created by Dr. Jon Purtiz that is specifically designed for different types of RAD sequencing.

**Jon downloaded a version of dDocent on my KITT account `dDocent_ngs` that can be used for Expressed Exome Capture Sequencing (EecSeq). It is located in the `NB_capture` directory**

If this is your first time running dDocent, I recommed going through the [Quick Start Guide](http://www.ddocent.com/quick/). I also recommend:

1. Reading through the [User Guide](http://www.ddocent.com/UserGuide/).
2. Completing the [Assembly Tutorial](http://www.ddocent.com/assembly/), using the simulated dataset.

**Create and activate a dDocent conda environment (if you did not do so previously):**

```
$ conda create -n nb_capture ddocent
$ conda activate nb_capture

# the beginning of the line should then look like: (nb_capture) [azyck@KITT ~]$
```

**Make directory for dDocent and link files and `dDocent_ngs` into directory**

Starting in `NB_capture` directory:

```
$ mkdir NB_ddocent
$ cd NB_ddocent/

$ ln -s ../*.fq.gz .
$ ln -s ../dDocent_ngs .
```

**First run `dDocent_ngs` with just Read Trimming**
```
$ bash dDocent_ngs
```

```
dDocent 2.8.10

Contact jpuritz@uri.edu with any problems


Checking for required software

All required software is installed!

dDocent version 2.8.10 started Fri Feb 21 17:54:22 EST 2020

50 individuals are detected. Is this correct? Enter yes or no and press [ENTER]
$ yes
Proceeding with 50 individuals
dDocent detects 80 processors available on this system.
Please enter the maximum number of processors to use for this analysis.
$ 20

Do you want to quality trim your reads?
Type yes or no and press [ENTER]?
$ yes

Do you want to perform an assembly?
Type yes or no and press [ENTER].
$ no

Reference contigs need to be in a file named reference.fasta

Do you want to map reads?  Type yes or no and press [ENTER]
$ no
Mapping will not be performed
Do you want to use FreeBayes to call SNPs?  Please type yes or no and press [ENTER]
$ no

Please enter your email address.  dDocent will email you when it is finished running.
Don't worry; dDocent has no financial need to sell your email address to spammers.


At this point, all configuration information has been entered and dDocent may take several hours to run.
It is recommended that you move this script to a background operation and disable terminal input and output.
All data and logfiles will still be recorded.
To do this:
Press control and Z simultaneously
Type bg and press enter
Type disown -h and press enter
Now sit back, relax, and wait for your analysis to finish

Trimming reads

dDocent has finished with an analysis in /home/azyck/NB_capture/NB_ddocent

dDocent started Fri Feb 21 17:54:22 EST 2020

dDocent finished Fri Feb 21 18:19:20 EST 2020
```

**Before continuing with dDocent, I am checking the quality of the trimmed reads with FastQC**

Back in `NB_capture` Directory

```
$ mkdir cleanedreads_fastqc_results
$ cd cleanedreads_fastqc_results/

$ ln -s ../NB_ddocent/*.fq.gz .
```

```
# fastqc takes awhile to run on larger data sets, so I recommend running it overnight

$ fastqc *.R1.fq.gz
$ fastqc *.R2.fq.gz
```

```
$ multiqc .

$ mv multiqc_report.html cleanedreads_multiqc_report.html

$ scp -r -P XXXX x@kitt.uri.edu:/home/azyck/NB_capture/cleanedreads_fastqc_results/cleanedreads_multiqc_report.html /Users/azyck/Documents/Research/NB_AdultOyster/Analysis_Outputs/fastqc
```

Multiqc report can be viewed [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Output/multiqc/cleanedreads_multiqc_report.html).


**Download _Crassostrea virginica_ genome**

Jon added the Eastern Oyster genome to the `NB_ddocent` Directory for me. Resources for the Eastern Oyster genome are on [NCBI](https://www.ncbi.nlm.nih.gov/Traces/wgs/MWPT03?display=contigs&page=1).

```
# In NB_ddocent directory
$ wget #paste link here

# Reference contigs need to be in a file named reference.fasta
$ mv #file_name_here reference.fasta
```

**Read Mapping**

Read mapping has to be performed at least once before SNP calling. If trimming and assembly remains the same when using dDocent repeatly, then read mapping is not necessary. A (match score), B (mismatch score), and 0 (gap penalty) parameters are used.

*Match Score*: score for a matching base.

*Mismatch Score*: score for a mismatching base.

*Gap Penalty*: Mutations, caused by either insertions or deletions, are annotated as gaps in the sequence and refered to as indels. These gaps are represented as dashes in the alignment. Optimization of gap penalty allows for the most accurate alignment possible. This will avoid low scores in alignments. Too small of a gap penalty will prevent a high level of alignment, and too large of a gap penalty will prevent accurate alignments.

```
$ bash dDocent_ngs
```

```
dDocent 2.8.10

Contact jpuritz@uri.edu with any problems


Checking for required software

All required software is installed!

dDocent version 2.8.10 started Sun Feb 23 14:03:01 EST 2020

50 individuals are detected. Is this correct? Enter yes or no and press [ENTER]
$ yes
Proceeding with 50 individuals
dDocent detects 80 processors available on this system.
Please enter the maximum number of processors to use for this analysis.
$ 20

Do you want to quality trim your reads?
Type yes or no and press [ENTER]?
$ no # reads only need to be trimmed once

Do you want to perform an assembly?
Type yes or no and press [ENTER].
$ no # We have a reference genome for the Eastern Oyster

Reference contigs need to be in a file named reference.fasta

Do you want to map reads?  Type yes or no and press [ENTER]
$ yes
BWA will be used to map reads.  You may need to adjust -A -B and -O parameters for your taxa.
Would you like to enter a new parameters now? Type yes or no and press [ENTER]
$ yes
Please enter new value for A (match score).  It should be an integer.  Default is 1.
$ 1
Please enter new value for B (mismatch score).  It should be an integer.  Default is 4.
$ 3
Please enter new value for O (gap penalty).  It should be an integer.  Default is 6.
$ 5
Do you want to use FreeBayes to call SNPs?  Please type yes or no and press [ENTER]
$ no

Please enter your email address.  dDocent will email you when it is finished running.
Don't worry; dDocent has no financial need to sell your email address to spammers.


At this point, all configuration information has been entered and dDocent may take several hours to run.
It is recommended that you move this script to a background operation and disable terminal input and output.
All data and logfiles will still be recorded.
To do this:
Press control and Z simultaneously
Type bg and press enter
Type disown -h and press enter
Now sit back, relax, and wait for your analysis to finish
Using BWA to map reads.
```

> Read mapping refers to aligning reads to a reference sequence. [dDocent](http://www.ddocent.com/UserGuide/#read-mapping) uses the program [BWA](http://bio-bwa.sourceforge.net/) to output a SAM (Sequence Alignment/Map) file, which is then converted to a BAM (compact) file using `SAMtools`.
> The default settings are optimized for the human genome and are generally conservative. Finding optimal values for each data set is recommended.

**Mark Duplicates**

**_Note:_** The code for marking duplicates is contained within `dDocent_ngs`, so it is not necessary to re-do these steps. I did not realize this, so I accidentally repeated the steps.

Following the Mark Duplicates section of the [OysterGenomeProject](https://github.com/jpuritz/OysterGenomeProject/blob/master/Bioinformatics/OysterGenome.md) repository on Github, I am marking and removing duplicate sequences.

```
# Downloading Picard script from jpuritz
$ curl -L -O https://raw.githubusercontent.com/jpuritz/OysterGenomeProject/master/Bioinformatics/Scripts/picard.sh
$ chmod +x picard.sh

$ cat picard.sh
```

```
## #!/usr/bin/env bash
##
## cat namelist | parallel -j 8 "java -Xms4g -jar /usr/local/bin/picard.jar MarkDuplicates I={}-RG.bam O={}-RGmd.bam M={}_dup_metrics.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 TAGGING_POLICY=OpticalOnly &> md.{}.log"
##
## echo -e "Picard has finished  in" `pwd` | mailx -s "Analysis has finished" jpuritz@uri.edu #change to your email
```

```
$ bash picard.sh
```

Now that duplicates are marked, duplicate sequences and secondary alignments are removed.

```
$ filter_bam(){
> samtools view -@32 -h -F 0x100 -F 0x400 $1-RGmd.bam | mawk '$6 !~/[8-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/'| samtools view -@ 32 -b
> }

$ export -f filter_bam
$ cat namelist | parallel -j 8 "filter_bam {} > {}.F.bam"
```

Now, megre all the individual bam files into a single bam file for easier parallelized variant calling.

```
$ ls *.F.bam > bam.list
$ samtools merge -@64 -f filter.merged.bam -b bam.list
$ samtools index -@64 filter.merged.bam
```

Use `samtools flagstat` to investigate `bam` files.

https://github.com/bahlolab/bioinfotools/blob/master/SAMtools/flagstat.md

```
$ samtools flagstat GB_6.F.bam

output:
14080016 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
14080016 + 0 mapped (100.00% : N/A)
14080016 + 0 paired in sequencing
7051282 + 0 read1
7028734 + 0 read2
13740006 + 0 properly paired (97.59% : N/A)
14060514 + 0 with itself and mate mapped
19502 + 0 singletons (0.14% : N/A)
209902 + 0 with mate mapped to a different chr
209902 + 0 with mate mapped to a different chr (mapQ>=5)
```

```
$ samtools flagstat NIN_10.F.bam

output:
10251730 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
10251730 + 0 mapped (100.00% : N/A)
10251730 + 0 paired in sequencing
5135084 + 0 read1
5116646 + 0 read2
10030216 + 0 properly paired (97.84% : N/A)
10239207 + 0 with itself and mate mapped
12523 + 0 singletons (0.12% : N/A)
127226 + 0 with mate mapped to a different chr
127226 + 0 with mate mapped to a different chr (mapQ>=5)
```

I looked at a few other `bam` files. All had 100% mapped. A low percentage of the mappings were singletons (only one read from a pair). Percentage of proper pairings was between 97 and 99%.


**Variant Calling**

Run dDocent again to generate the variant calls and then zip the huge VCF file!

```
$ bash dDocent_ngs
```

```
dDocent 2.8.10

Contact jpuritz@uri.edu with any problems


Checking for required software

All required software is installed!

dDocent version 2.8.10 started Tue Feb 25 10:55:46 EST 2020

50 individuals are detected. Is this correct? Enter yes or no and press [ENTER]
$ yes
Proceeding with 50 individuals
dDocent detects 80 processors available on this system.
Please enter the maximum number of processors to use for this analysis.
$ 20

Do you want to quality trim your reads?
Type yes or no and press [ENTER]?
$ no

Do you want to perform an assembly?
Type yes or no and press [ENTER].
$ no

Reference contigs need to be in a file named reference.fasta

Do you want to map reads?  Type yes or no and press [ENTER]
$ no
Mapping will not be performed
Do you want to use FreeBayes to call SNPs?  Please type yes or no and press [ENTER]
$ yes

Please enter your email address.  dDocent will email you when it is finished running.
Don't worry; dDocent has no financial need to sell your email address to spammers.


At this point, all configuration information has been entered and dDocent may take several hours to run.
It is recommended that you move this script to a background operation and disable terminal input and output.
All data and logfiles will still be recorded.
To do this:
Press control and Z simultaneously
Type bg and press enter
Type disown -h and press enter
Now sit back, relax, and wait for your analysis to finish

Using FreeBayes to call SNPs

After filtering, kept 4017690 out of a possible 52332265 Sites
```

> FreeBayes detects variants based on differences in haplotypes, including single nucleotide polymorphisms (SNPs), and indels (insertions and deletions). [FreeBayes](https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionVI/lessons/02_variant-calling.html).

### I re ran `dDocent_ngs` with a new reference.fasta genome where haplotigs are masked.

Resulting output after trimming, read mapping, and SNP calling:
```
Using FreeBayes to call SNPs

After filtering, kept 11314456 out of a possible 48598546 Sites
```


Zip the VCF file. This was completed for the original VCF file and again on the masked haplotig VCF file.

_`dDocent_ngs` may have a step at the end of the script where the VCF file is zipped, so you may not need to do this step_

```
$ gzip TotalRawSNPs.VCF
```

This produced a massive raw variants file!

```
$ vcftools --gzvcf TotalRawSNPs.vcf.gz

output:
VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--gzvcf TotalRawSNPs.vcf.gz

Using zlib version: 1.2.11
After filtering, kept 50 out of 50 Individuals
After filtering, kept 52332265 out of a possible 52332265 Sites
Run Time = 608.00 seconds
```

```
Masked haplotig output:
After filtering, kept 50 out of 50 Individuals
After filtering, kept 48598546 out of a possible 48598546 Sites
Run Time = 653.00 seconds
```

Filtering for each VCF file was then performed and documented [here](https://github.com/amyzyck/EecSeq_NB_EasternOyster/blob/master/Analysis/EecSeq_Cvirginica_Filtering.md)
