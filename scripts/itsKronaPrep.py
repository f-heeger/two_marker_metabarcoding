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
    sample = outPath.split("/")[1][:-15]
    with open(outPath, "w") as out:
        clsCount = {}
        for read in sampleReads[sample]:
            if read not in readCls:
                #this are reads that were filtered out after the sample
                # file was created
                continue
            try:
                clsCount[readCls[read]] += 1
            except KeyError:
                clsCount[readCls[read]] = 1
        for cls, count in clsCount.items():
            out.write("%i\t%s\n" % (count, "\t".join(cls.split(";"))))
