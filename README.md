# Two marker metabarcoding pipeline
Snakemake pipeline for analysis of metabarcoding data of fungi with more than one marker (5.8S and ITS2).

The pipeline (v1.1) is described in out paper [Combining the 5.8S and ITS2 to improve classification of fungi](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13266). The current version (v1.2) differs from that in that it uses snakemakes conda integration for better portability and easy installation.

There is also a pre-print avialable on BioRxiv descibing version 1.0 of this pipeline: [https://doi.org/10.1101/532358 ](https://www.biorxiv.org/content/10.1101/532358v1)

## Prerequisites


### Software you need to install
You need to install a version of Conda, that will be used to install all other software automatically except Snakemake and Python. You can use Conda to install Snakemake (which will automatically install Python).

* [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
* Python 3.x
* [Snakemake (version 3.5.4 or newer)](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

### Software that will be installed automatically
The following tools will be automatically installed by Snakemake trough conda. Make sure that Conda is in your path, when you run Snakemake.

* [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [MultiQC](http://multiqc.info/)
* [VSEARCH](https://github.com/torognes/vsearch)
* [ITSx](http://microbiology.se/software/itsx/)
* [cutadapt](https://github.com/marcelm/cutadapt)
* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [Krona](https://github.com/marbl/Krona/wiki/KronaTools)
* [swarm](https://github.com/torognes/swarm)
* [Lambda](http://seqan.github.io/lambda/)

### Reference Data
Reference data will be downloaded and processed automatically

* [UNITE](https://unite.ut.ee/) ITS database
* [RFAM](http://rfam.xfam.org/) family RF00002 (5.8S rRNA) 

## Installation

### Preparing your working directory

You can directly download the files into your working directory or clone the repository with git. Your working folder should contain the files beside this readme:

   * `perpDatabases.snakemake.py`
   * `metabarcoding.snakemake.py`
   * `final.snakemake.py`
   * `its.snakemake.py`
   * `r58S.snakemake.py`
   * `readProcessing.snakemake.py`
   * `config.json`

### Setting up your configuration file

The config file is in json format and is a list of keys (or names) and values separated by a ":".

The configuration file gives the paths to all necessary data as well as some information about your run and additional configuration. You need to change the information about your run, the software paths and your e-mail address (for querrying the NCBI database).


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
* **primerError** How many errors (given as error rate i.e. 0.1=10%) should be allowed when checking if forward and reverse primer are present.
* **maxAmplLen** Maximum length of amplicons without the amplification primers. Will be passed to vsearch for merging.
* **minAmplLen** Minimum length of amplicons without the amplification primers. Will be passed to vsearch for merging.

#### Reference Data
* **unite_version**: Release date of the UNITE database you want to use (default: 01.12.2017 i.e. version 7.2)
* **uniteUrl**: URL were the UNITE database should be downloaded from
* **rfam_version**: Version of the RFAM database you want to use (default: 13.0)



## Running the pipeline

* Open a terminal and navigate to your working directory.
* type the following command: `snakemake --use-conda -s metabarcoding.snakemake.py -j 6` where "6" is the number of processors to use

## Hidden Features
Some output files are not generated by default. Here is how to get them:

* Concatenated 5.8S and ITS2 sequence for each OTU representative: run `snakemake -s metabarcoding.snakemake.py all.concatMarkers.fasta` to generate `all.concatMarkers.fasta` in the main folder. Some OTUs might be skipped if no 5.8S could be extracted for the representative.

