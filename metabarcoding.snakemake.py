#import xml.etree.ElementTree as et

from snakemake.utils import min_version, R
from Bio import SeqIO

min_version("3.5.4")

configfile: "config.json"

shell.prefix("sleep 10; ") #work aorund to desl with "too quck" rule execution and slow NAS


samples = config["samples"].keys()

rule all:
    input: "krona/All.krona.html", "krona/5_8s.krona.html", "krona/ITS2.krona.html", "taxonomy/all.compareClass.tsv", "otu_table.tsv", "All.rarefactions.pdf"

### generate reference data bases
include: "rules/prepDatabases.snakefile.py"

### read processing
include: "rules/readProcessing.snakemake.py"

### analize 5.8S
include: "rules/r58S.snakemake.py"

### analize ITS2
include: "rules/its.snakemake.py"

### combine 5.8S and ITS2 and analize result
include "rules/rules/final.snakemake.py"



