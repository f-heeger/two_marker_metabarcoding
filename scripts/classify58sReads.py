repseq = {}
for line in open(snakemake.input.repseq):
    read, rep = line.strip().split("\t")
    try:
        repseq[rep].append(read)
    except KeyError:
        repseq[rep] = [read]

rep58s = {}
for line in open(snakemake.input.rep58s):
    seq, rep = line.strip().split("\t")
    try:
        rep58s[rep].append(seq)
    except KeyError:
        rep58s[rep] = [seq]

readClass = {}
for line in open(snakemake.input.tax):
    rId, classification = line.strip().split("\t")
    cls = []
    for entry in classification.strip(";").split(";"):
        if entry == "unclassified":
            break
        elif entry == "unknown":
            cls = [""]
        else:
            cls.append(entry)
    for r58seq in rep58s["%s|5.8S" % rId]:
        for read in repseq[r58seq.split("|")[0]]:
            readClass[read] = ";".join(cls)

with open(snakemake.output[0], "w") as out:
    for read, cls in readClass.items():
        out.write("%s\t%s\n" % (read, cls))
