readNr = {}
for sample in snakemake.params.samples:
    readNr[sample] = {}
    for line in open("swarm/%s.ITS2.otus.out" % sample):
        otuStr, nr = line.strip("\n").split("\t")
        otu = otuStr.split("|")[0]
        readNr[sample][otu] = nr

r58sCls = {}
for line in open(snakemake.input.tax58s):
    name, lin, nr = line.strip().split("\t")
    r58sCls[name.split("|")[0]] = lin

itsCls = {}
for line in open(snakemake.input.taxIts):
    name, lin = line.strip().split("\t")
    itsCls[name.split("|")[0]] = lin

with open(output[0], "w") as out:
    out.write("otu_ID\t5.8S classification\tITS2 classification\tfinal classification\t%s\n" % "\t".join(samples))
    for line in open(snakemake.input.taxComb):
        otu, cls = line.strip("\n").split("\t")
        numbers = [readNr[sample].get(otu, "0") for sample in samples]
        out.write("%s\t%s\t%s\t%s\t%s\n" % (otu, r58sCls[otu], itsCls[otu], cls, "\t".join(numbers)))
