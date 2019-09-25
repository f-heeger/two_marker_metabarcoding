sampleReads = {}

for line in open(snakemake.input.sample):
    read, sample = line.strip("\n").split("\t")
    try:
        sampleReads[sample].append(read)
    except KeyError:
        sampleReads[sample] = [read]

readCls = {}
for line in open(snakemake.input.cls):
    read, cls = line.strip("\n").split("\t")
    readCls[read] = cls

for outPath in snakemake.output:
    sample = outPath.split("/")[1][:-9]
    with open(outPath, "w") as out:
        clsCount = {}
        for read in sampleReads[sample]:
            try:
                cls = readCls[read]
            except KeyError:
                #reads that were removed in between
                continue
            try:
                clsCount[cls] += 1
            except KeyError:
                clsCount[cls] = 1
        for cls, count in clsCount.items():
            out.write("%i\t%s\n" % (count, "\t".join(cls.split(";"))))
