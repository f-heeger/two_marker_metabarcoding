import gzip

from Bio import SeqIO

with gzip.open(snakemake.output.comb, "wt") as combOut, \
     open(snakemake.output.sample, "w") as sampleOut, \
     open(snakemake.output.name, "w") as nameOut:
    for sample, fileInfo in snakemake.params.files.items():
        for inFile in fileInfo[int(snakemake.wildcards.read[1])]:
            with gzip.open("%s/%s" % (snakemake.config["inFolder"], inFile), "rt") as inStream:
                for rec in SeqIO.parse(inStream, "fastq"):
                    newId = "_".join(rec.id.split(":")[3:])
                    nameOut.write("%s\t%s\n" % (rec.id, newId))
                    rec.id = newId
                    combOut.write(rec.format("fastq"))
                    sampleOut.write("%s\t%s\n" % (newId, sample))
