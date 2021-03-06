
from snakemake.utils import min_version

min_version("3.5.4")

configfile: "config.json"

fileInfo = {}
for line in open(config["sampleFile"]):
    if line[0] == "#":
        continue
    sName, rNr, fileName = line.strip("\n").split("\t")
    if sName not in fileInfo:
        fileInfo[sName] = [[], []]
    fileInfo[sName][int(rNr)-1].append(fileName)

samples = list(fileInfo.keys())

rule all:
    input: "krona/All.krona.html", "krona/5_8s.krona.html", "krona/ITS2.krona.html", "taxonomy/all.compareClass.tsv", "otu_table.tsv", "All.rarefactions.pdf", "readNumbers/readNumbers.pdf"

### generate reference data bases
include: "prepDatabases.snakemake.py"

### read processing
include: "readProcessing.snakemake.py"

### analize 5.8S
include: "r58S.snakemake.py"

### analize ITS2
include: "its.snakemake.py"

### combine 5.8S and ITS2 and analize result
include: "final.snakemake.py"



