from Bio import SeqIO

r58s = {}
for r in SeqIO.parse(open(snakemake.input.tsu), "fasta"):
    r58s[r.id.split("|", 1)[0]] = r
notFound = 0

with open(snakemake.output[0], "w") as out, open(snakemake.log[0], "w") as logStream:
    for otu in SeqIO.parse(open(snakemake.input.otus), "fasta"):
        newId = otu.id.split("|", 1)[0]
        try:
            tsu = r58s[newId]
        except KeyError:
            logStream.write("Did not find 5.8S for %s\n" % newId)
            notFound += 1
        concat = tsu + otu
        concat.id = "%s|partial5.8S+ITS2" % newId
        concat.description = ""
        out.write(concat.format("fasta"))

print("Could not find 5.8S for %i OTUs. See %s for details." % (notFound, snakemake.log[0]))
