import gzip

from Bio import SeqIO

readSample = {}
for line in open(snakemake.input.sample):
    read, sample = line.strip("\n").split("\t")
    readSample[read] = sample
sampleReads = {}
for rec in SeqIO.parse(gzip.open(snakemake.input.reads, "rt"), "fastq"):
    sample = readSample[rec.id]
    try:
        sampleReads[sample] += 1
    except KeyError:
        sampleReads[sample] = 1
with open(snakemake.output[0], "w") as out:
    for sample, value in sampleReads.items():
        out.write("%s\t%s\n" % (sample, value))
