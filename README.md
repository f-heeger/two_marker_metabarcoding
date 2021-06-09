# Two marker metabarcoding pipeline
Snakemake pipeline for analysis of metabarcoding data of fungi with more than one marker (5.8S and ITS2).

The pipeline (v1.1) is described in out paper [Combining the 5.8S and ITS2 to improve classification of fungi](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13266). 

The current version (v1.4) can be found on the [release page](https://github.com/f-heeger/two_marker_metabarcoding/releases). It differs from version in the paper in that it uses snakemakes conda integration for better portability and easy installation. Additionally input file handling is handled in a way that is less reliant on default Illumina file names and allows for more complicated experimental setups (e.g. the same sample in multiple runs). I also changed some external dependencies. See the release notes for further detail.

There is also a pre-print available on BioRxiv describing version 1.0 of this pipeline: [https://doi.org/10.1101/532358 ](https://www.biorxiv.org/content/10.1101/532358v1)

## Prerequisites


### Software you need to install
You need to install a version of Conda, that will be used to install all other software automatically except Snakemake and Python. You can use Conda to install Snakemake (which will automatically installs Python).

* [Miniconda](https://docs.conda.io/en/latest/miniconda.html) (or another flavor of conda)
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

You can directly download the files into your working directory or clone the repository with git. Your working folder should contain the following files beside this readme:

   * `envs`
   * `scripts`
   * `testData`
   * `LICENSE`
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

The configuration file gives the paths to all necessary data as well as some information about your run and additional configuration. You need to change the information about your run, the software paths and your e-mail address (for querying the NCBI database).


#### Your e-mail adresse:
to build the database the pipeline will query the NCBI taxonomy database. NCBI requires automated API calls like this to also send an e-mail address. Your e-mail address will not be used for anything else and will be send nowhere except NCBI.

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
    1. `mark`: The conflict is marked in the classification with an entry of the form: \<ITS classification>|<5.8S classification>
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

The repository also contains a bash script called `run.sh`. I recommend using the above command unless you are familiar with bash scripting and conda. If you want to use it, you have to

* edit the first line to give the correct path to your conda installation (or have conda in your path anyway)
* create a conda environment called `sankemake5.4.3` with sankemake (v. 5.4.3) installed in it or edit line 4 to activate a environment in which you have installed snakemake

All parameters you give to the `run.sh` script (e.g. `-j 6`) will be passed on to snakemake.

## Output Files
The pipeline creates different folders for the intermediate files after each step and the final result files.
The most important result is `otu_table.tsv` that is created in the working directory. It contains the number of reads that were assigned to the OTUs in each sample as well as the taxonomic classification of the OTUs. Each row represents one OTU. The first column is the OTU ID, the second, third and forth column are the taxonomic classifications by 5.8S, ITS2 and both combined respectively. The following columns contain the number of reads per sample (see header for sample order).

Other interesting outputs are the rarefaction plot (`All.rarefactions.pdf`), the plot of read numbers after each filtering step during the initial read processing (`readNumbers/readNumbers.pdf`) and the [Krona](https://github.com/marbl/Krona/wiki) plots of the read assignments per sample (`krona/All.krona.html`).


## Hidden Features
There are some non-standard features, that you will not need for a typical analysis, but that I will document here:

* **Concatenated 5.8S and ITS2 sequence for each OTU representative**: run `snakemake -s metabarcoding.snakemake.py all.concatMarkers.fasta` to generate `all.concatMarkers.fasta` in the main folder. Some OTUs might be skipped if no 5.8S could be extracted for the representative.
* **Old UNITE versions**: Unite keeps changing the format of their download, so for now I can only make the automatic download work for the newest version. To use an older version you can try the following workaround:
    * Create the folder you gave for `dbFolder` in the config file (if it does not exists already)
    * Set the `unite_version` in the config file to the correct date string (e.g. 01.08.2015)
    * The`uniteUrl` in the config file will not be used. You can leave it as it is (or change it to something to remind you that it is ignored)
    * Manually download the [general FASTA release](https://unite.ut.ee/repository.php) of UNITE that you want to use and unpack it
    * Create a folder called `sh_general_release_all_<unite_version>` in your `dbFolder` where `<unite_verstion>` is the version string you gave for the `unite_version` in the config file (e.g. 01.08.2015)
    * copy the fasta file, that was in your download from the UNITE website into the new folder and rename it to: `sh_general_release_dynamic_all_<unite_version>.fasta` where `<unite_verstion>` is the version string you gave for the `unite_version` in the config file (e.g. 01.08.2015)
* **Frequency Table**: run `snakemake -s metabarcoding.snakemake.py freq_table.tsv` to create `freq_table.tsv`, which is the same as the OTU table except, that is does not contain taxonomic classifications and the fact that the read numbers are normalized to the total number of reads in that sample (i.e. relative read counts considering total number of reads in that sample instead of absolute read counts).