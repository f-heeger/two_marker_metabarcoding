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
   * `samples.tsv`

### Setting up your configuration file

The config file is in json format and is a list of keys (or names) and values separated by a ":".

The configuration file gives the paths to all necessary data as well as some information about your run and additional configuration. You need to change the information about your run, the software paths and your e-mail address (for querrying the NCBI database).


#### Your e-mail adresse:
to build the database the pipeline will query the NCBI taxonomy database. NCBI requires automated API calls like this to also send an e-mail address. Your e-mail address will not be used for anything else.

#### Run the pipeline with test data:
After setting your e-mail address, you can run the pipeline on the provided test data to see if it is working correctly. This will also already download the necessary reference data bases (see Reference Data below on how to set the database version) and install needed software via conda.

To do a test run, open a terminal, navigate to the pipeline working directory and type the following command: `snakemake --use-conda -s metabarcoding.snakemake.py`. After everything is done you should see the message: `72 of 72  steps (100%) done`.

#### Information about your run:

* **inFolder**: the absolute path to the folder that contains the fastq files of the raw reads (after demultiplexing)
* **dbFolder**: where should the databases be stored. If you leave this as the default value (`dbs`) it will create a folder called `dbs` in your work directory. This is fine, but because the database preparation step will take a lot of time it might be good to setup the database in a central location and use it for all you runs.
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
* **unite_version**: Release date of the UNITE database you want to use (default: 04.02.2020 i.e. version 8.2). This will only be used for internal naming. The url you give in the next option determines the actual version used.
* **uniteUrl**: URL were the UNITE database should be downloaded from. Note that Unite frequently changes the names and properties of files, which makes it hard to keep the automatic download working for all version. For now changing the version to an older version will unfortunately most likely break the download. See below for a work around to use older versions.
* **rfam_version**: Version of the RFAM database you want to use (default: 14.2)

### Setting up your samples file
The sample file is a table in tab-separated-values format providing information about your input files. Lines with a `#` in front of them will be ignored. One file is described per row with the following columns:

* **sample name**: a human readable sample name that this file is assigned to
* **read number**: does this file contain forward (`1`) or reverse (`2`) reads
* **file name**: the actual name of the input file (it is expected to be located in the folder given in the config file under `inputFolder`, so only give the file name, NOT the full path).

## Running the pipeline

* Open a terminal and navigate to your working directory.
* type the following command: `snakemake --use-conda -s metabarcoding.snakemake.py -j 6` where "6" is the number of processors to use

## Hidden Features
Some output files are not generated by default. Here is how to get them:

* **Concatenated 5.8S and ITS2 sequence for each OTU representative**: run `snakemake -s metabarcoding.snakemake.py all.concatMarkers.fasta` to generate `all.concatMarkers.fasta` in the main folder. Some OTUs might be skipped if no 5.8S could be extracted for the representative.
* **Old UNITE versions**: Unite keeps changing the format of their download, so for now I can only make the automatic download work for the newest version. To use an older version you can try the following workaround:
    * Create the folder you gave for `dbFolder` in the config file (if it does not exists already)
    * Set the `unite_version` in the config file to the correct date string (e.g. 01.08.2015)
    * The`uniteUrl` in the config file will not be used. You can leave it as it is (or change it to something to remind you that it is ignored)
    * Manually download the [general FASTA release](https://unite.ut.ee/repository.php) of UNITE that you want to use and unpack it
    * Create a folder called `sh_general_release_all_<unite_version>` in your `dbFolder` where `<unite_verstion>` is the version string you gave for the `unite_version` in the config file (e.g. 01.08.2015)
    * copy the fasta file, that was in your download from the UNITE website into the new folder and rename it to: `sh_general_release_dynamic_all_<unite_version>.fasta` where `<unite_verstion>` is the version string you gave for the `unite_version` in the config file (e.g. 01.08.2015)

