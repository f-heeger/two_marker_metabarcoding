import gzip
import glob

from Bio import SeqIO

with gzip.open(snakemake.output.comb, "wt") as combOut, \
     open(snakemake.output.sample, "w") as sampleOut, \
     open(snakemake.output.name, "w") as nameOut:
    for sample in snakemake.params.samples:
        path = "%s/%s_L*_%s_*.fastq.gz" % (snakemake.config["inFolder"], sample, snakemake.wildcards.read)
        inFiles=glob.glob(path)
        if not inFiles:
            raise RuntimeError("No file(s) found for sample %s at %s." % (sample, path))
        for inFile in inFiles:
            with gzip.open(inFile, "rt") as inStream:
                for rec in SeqIO.parse(inStream, "fastq"):
                    newId = "_".join(rec.id.split(":")[3:])
                    nameOut.write("%s\t%s\n" % (rec.id, newId))
                    rec.id = newId
                    combOut.write(rec.format("fastq"))
                    sampleOut.write("%s\t%s\n" % (newId, sample))
