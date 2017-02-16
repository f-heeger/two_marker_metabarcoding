# Prerequisites

## Python and Python libraries

* Python 3.x
* BioPython
* rpy2
* snakemake (version 3.5.4 or newer)

## External Software

* [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [MultiQC](http://multiqc.info/) (optional)
* [VSEARCH](https://github.com/torognes/vsearch)
* [ITSx](http://microbiology.se/software/itsx/)
* [cutadapt](https://github.com/marcelm/cutadapt)
* [mothur](http://mothur.org/)
* [Flexbar](https://github.com/seqan/flexbar)
* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [Pear](http://sco.h-its.org/exelixis/web/software/pear/)
* [Krona](https://github.com/marbl/Krona/wiki/KronaTools)
* [swarm](https://github.com/torognes/swarm)
* [Lambda](http://seqan.github.io/lambda/) version 0.9 or 1.0 (the 1.9.* version is not supported yet)

## Reference Databases

* UNITE ITS database (will be downloaded automatically)
* 5.8S reference alignment (alignment in fasta format and taxonomy information in tsv format)

# Installation
Make sure that you installed all the prerequisites listed above.

## getting the pipeline files

## preparing your working directory

1. Your working folder should contain two files beside this readme:
    * `metabarcoding.snakemake.py`
    * `config.json`
2. create a folder inside the working directory to hold database files
3. put the 5.8S reference files into the database folder that you created

## setting up your configuration file

The configuration file gives the paths to all necessary resources (software and data) as well as some information about your run.
The config file is in json format and is a list of keys (or names) and values separated by a ":".
For each software that is used during pipeline execution there is an entry that should contain the absolute path to the softwares binary or the command the software can be executed with. There are two exceptions to this. For Lambda the entry should be the path to the folder containing the lambda binary as well as the lambda-indexer binary. For Trimmomatic the path should be given for the Trimmomatic `.jar` file.

In addition to the software paths there are the following entries that give information about your run:

1. **inFolder**: the absolut path to the folder that contains the fastq files of the raw reads (after demultiplexing)
2. **dbFolder**: the name of the folder you created to hold your databases (step 3. in the last section)
3. **samples**: a list of key-value pairs describing your data. The list is enclosed by "{" and each entry is separated by a ",". They key gives the sample ID as it was used used in the sequencing run (what the files are called), while the value gives the names you want the samples identified by in the final result. Example: `{"A1_S1": "lake1_sample1", "A3_S2": "lake1_sample2", "B12_S3": "lake2_sample1"}`
4. **forward_primer**: the sequence of your forward amplification primer
5. **reverse_primer**: the sequence of your reverse amplification primer

# Configuration

The config file also contains optional configurations, which can be used to change the behavior of the pipeline.

* **conflictBehavior** This parameter defines what the pipeline will do if the 5.8S and ITS2 classifications differ from each other. There are three possible values:
    1. `mark`: The conflict is marked in the classification with an entry of the form: <ITS classification>|<5.8S classification>
    2. `5.8S`: Use the 5.8S classification and ignore the ITS classification
    3. `ITS`: Use the ITS2 classification and ignore the 5.8S classification
* **primerError** How many errors should be allowed when checking if forward and reverse primer are present.
* **maxAmplLen** Maximum length of amplicons without the amplification primers. Will be passed to pears `-m` option.
* **minAmplLen** Minimum length of amplicons without the amplification primers. Will be passed to pears `-n` option.

# Running the pipeline

* Open a terminal and navigate to your working directory.
* type the following command: `snakemake -s metabarcoding.snakemake.py -j 6` where "6" is the number of processors to use

# Hidden Features
Some output files are not generated by default. Here is how to get them:

* Concatenated 5.8S and ITS2 sequence for each OTU representative: run `snakemake -s metabarcoding.snakemake.py all.concatMarkers.fasta` to generate `all.concatMarkers.fasta` in the main folder. Some OTUs might be skipped if no 5.8S could be extracted for the representative.

