# Two marker metabarcoding pipeline
Snakemake pipeline for analysis of metabarcoding data of fungi with more than one marker (5.8S and ITS2).

There is a pre-print avialable on BioRxiv descibing this pipeline: [https://doi.org/10.1101/532358 ](https://www.biorxiv.org/content/10.1101/532358v1)

## Prerequisites

### Python, Python libraries and Snakemake

* Python 3.x
* [Biopython](https://biopython.org/)
* [snakemake (version 3.5.4 or newer)](https://snakemake.readthedocs.io/en/stable)

### External Software

* [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [MultiQC](http://multiqc.info/)
* [VSEARCH](https://github.com/torognes/vsearch)
* [ITSx](http://microbiology.se/software/itsx/)
* [cutadapt](https://github.com/marcelm/cutadapt)
* [Flexbar](https://github.com/seqan/flexbar)
* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [Pear](http://sco.h-its.org/exelixis/web/software/pear/)
* [Krona](https://github.com/marbl/Krona/wiki/KronaTools)
* [swarm](https://github.com/torognes/swarm)
* [Lambda](http://seqan.github.io/lambda/) version 0.9 or 1.0 (the 1.9.* version is not supported yet)

### Reference Data
Reference data will be downloaded and processed automatically

* [UNITE](https://unite.ut.ee/) ITS database
* [RFAM](http://rfam.xfam.org/) family RF00002 (5.8S rRNA) 

## Installation
Make sure that you installed all the prerequisites listed above.

### Preparing your working directory

You can directly download the files into your working directory or clone the repository with git. 

### Setting up your configuration file

The config file is in json format and is a list of keys (or names) and values separated by a ":".

The configuration file gives the paths to all necessary resources (software and data) as well as some information about your run and additional configuration. You need to change the information about your run, the software paths and your e-mail address (for querrying the NCBI database).

#### Software paths:

For each software that is used during pipeline execution there is an entry that should contain the absolute path to the softwares binary or the command the software can be executed with. There are two exceptions to this. For Lambda the entry should be the path to the folder containing the lambda binary as well as the lambda-indexer binary. For Trimmomatic the path should be given for the Trimmomatic `.jar` file.

#### Your e-mail adresse:
to build the database the pipeline will query the NCBI taxonomy database. NCBI requires automated API calls like this to also send an e-mail address. Your e-mail address will not be used for anything else.

#### Run the pipeline with test data:
After setting software pathes and you e-mail address, you can run the pipeline on the provided test data to see if it is working correctly. This will also already download the neccessary reference data bases (see Reference Data below on how to set the database version).
To do a test run, open a terminal, navigate to the pipeline working directory and type the following command: `snakemake -s metabarcoding.snakemake.py`. After everything is done you should see the message: `72 of 72  steps (100%) done`.

#### Information about your run:

* **inFolder**: the absolute path to the folder that contains the fastq files of the raw reads (after demultiplexing)
* **dbFolder**: where should the databases be stored
* **samples**: a list of key-value pairs describing your data. The list is enclosed by "{" and each entry is separated by a ",". They key gives the sample ID as it was used used in the sequencing run (what the files are called), while the value gives the names you want the samples identified by in the final result. Example: `{"A1_S1": "lake1_sample1", "A3_S2": "lake1_sample2", "B12_S3": "lake2_sample1"}`
* **forward_primer**: the sequence of your forward amplification primer
* **reverse_primer**: the sequence of your reverse amplification primer

#### Configuration of pipeline behavior

* **conflictBehavior** This parameter defines what the pipeline will do if the 5.8S and ITS2 classifications differ from each other. There are three possible values:
    1. `mark`: The conflict is marked in the classification with an entry of the form: <ITS classification>|<5.8S classification>
    2. `5.8S`: Use the 5.8S classification and ignore the ITS classification
    3. `ITS`: Use the ITS2 classification and ignore the 5.8S classification
* **primerError** How many errors should be allowed when checking if forward and reverse primer are present.
* **maxAmplLen** Maximum length of amplicons without the amplification primers. Will be passed to pears `-m` option.
* **minAmplLen** Minimum length of amplicons without the amplification primers. Will be passed to pears `-n` option.

#### Reference Data
* **unite_version**: Release date of the UNITE database you want to use (default: 01.12.2017 i.e. version 7.2)
* **uniteUrl**: URL were the UNITE database should be downloaded from
* **rfam_version**: Version of the RFAM database you want to use (default: 13.0)



## Running the pipeline

* Open a terminal and navigate to your working directory.
* type the following command: `snakemake -s metabarcoding.snakemake.py -j 6` where "6" is the number of processors to use

## Hidden Features
Some output files are not generated by default. Here is how to get them:

* Concatenated 5.8S and ITS2 sequence for each OTU representative: run `snakemake -s metabarcoding.snakemake.py all.concatMarkers.fasta` to generate `all.concatMarkers.fasta` in the main folder. Some OTUs might be skipped if no 5.8S could be extracted for the representative.

