# *Crassostrea virginica* EecSeq Demultiplexing, dDocent, and VCF Filtering  

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

# the beginning of the line should then look like: (nb_capture) [azyck@KITT ~]$
```

**Creating and entering a directory for this project**

```
$ mkdir NB_capture # creates a directory
$ cd NB_capture # navigates into this new directory

# the beginning of the line should now look like: (nb_capture) [azyck@KITT NB_capture]$
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
$ samtools flagstat GB_6-RG.bam

output:
15355235 + 0 in total (QC-passed reads + QC-failed reads)
91421 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
15355235 + 0 mapped (100.00% : N/A)
15263814 + 0 paired in sequencing
7638408 + 0 read1
7625406 + 0 read2
14911958 + 0 properly paired (97.69% : N/A)
15245186 + 0 with itself and mate mapped
18628 + 0 singletons (0.12% : N/A)
212270 + 0 with mate mapped to a different chr
212270 + 0 with mate mapped to a different chr (mapQ>=5)
```

```
$ samtools flagstat NIN_10-RG.bam

output:
11255169 + 0 in total (QC-passed reads + QC-failed reads)
63875 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
11255169 + 0 mapped (100.00% : N/A)
11191294 + 0 paired in sequencing
5600706 + 0 read1
5590588 + 0 read2
10963069 + 0 properly paired (97.96% : N/A)
11179262 + 0 with itself and mate mapped
12032 + 0 singletons (0.11% : N/A)
127678 + 0 with mate mapped to a different chr
127678 + 0 with mate mapped to a different chr (mapQ>=5)
```

I looked at a few other `bam` files. All had 100% mapped. A low percentage of the mappings were singletons (only one read from a pair). Percentage of proper pairings was between 97 and 98%.


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

Zip the VCF file

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

### 5. Variant Filtering

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