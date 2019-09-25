readOtu = {}
for line in open(snakemake.input.otuReads):
    otuId, readId = line.strip("\n").split("\t")
    readOtu[readId] = otuId

sampleReads = {}
for line in open(snakemake.input.sample):
    read, sample = line.strip("\n").split("\t")
    try:
        sampleReads[sample].append(read)
    except KeyError:
        sampleReads[sample] = [read]

for outPath in snakemake.output:
    sample = outPath.split("/")[1][:-14]
    otuCount = {}
    for read in sampleReads[sample]:
        if read in readOtu:
            #some reads in the sample never made it to an OTU because 
            # they were filtered, but otherwise we do this
            try:
                otuCount[readOtu[read]] += 1
            except KeyError:
                otuCount[readOtu[read]] = 1
    with open(outPath, "w") as out:
        for otuId, count in otuCount.items():
            out.write("%s\t%i\n" % (otuId, count))
