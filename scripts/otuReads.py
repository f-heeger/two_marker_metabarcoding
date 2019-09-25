
repSeq = {}
for line in open(snakemake.input.repSeq):
    read, rep = line.strip().split("\t")
    try:
        repSeq[rep].append(read)
    except KeyError:
        repSeq[rep] = [read]

repIts = {}
for line in open(snakemake.input.repIts):
    seq, rep = line.strip().split("\t")
    try:
        repIts[rep.split("|")[0]].append(seq.split("|")[0])
    except KeyError:
        repIts[rep.split("|")[0]] = [seq.split("|")[0]]

otu = {}
for line in open(snakemake.input.otuList):
    memSeqs = line.strip().split(" ")
    otu[memSeqs[0].strip(";").split(";")[0]] = [s.split("|")[0] for s in memSeqs]

with open(snakemake.output[0], "w") as otuOut:
    for otuId in otu.keys():
        for itsSeq in otu[otuId]:
            for repSeqId in repIts[itsSeq]:
                for read in repSeq[repSeqId]:
                    otuOut.write("%s\t%s\n" % (otuId, read))
