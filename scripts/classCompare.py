readSample = {}
for line in open(snakemake.input.sample):
    read, sample = line.strip("\n").split("\t")
    readSample[read] = sample

otuReads = {}
for line in open(snakemake.input.otuReads):
    otuId, readId = line.strip().split("\t")
    try:
        otuReads[otuId].append(readId)
    except KeyError:
        otuReads[otuId] = [readId]

readCls58S = {}
for line in open(snakemake.input.r58SreadCls):
    read, cls = line.strip("\n").split("\t")
    if cls == "":
        readCls58S[read] = "unknown"
    else:
        readCls58S[read] = cls

otuCls58S = {}
for line in open(snakemake.input.r58sOtuCls):
    name, lin, number = line.strip("\n").split("\t")
    otuCls58S[name.split("|")[0]] = lin

itsClass = {}
for line in open(snakemake.input.itsCls):
    nameStr, cls = line.strip("\n").split("\t")
    itsClass[nameStr.split("|")[0]] = cls

combCls = {}
for line in open(snakemake.input.combCls):
    otuId, cls = line.strip("\n").split("\t")
    combCls[otuId] = cls

with open(snakemake.output.comp, "w") as out, open(snakemake.output.stat, "w") as stat:
    out.write("read\tsample\tOTU\tITS OTU Class\t5.8S read Class\t5.8S OTU Class\tComb OTU Class\n")
    for line in open(input.otuReads):
        otuIdStr, readId = line.strip("\n").split("\t")
        otuId = otuIdStr.split("|")[0]
        out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (readId,
                                        readSample[readId],
                                        otuId,
                                        itsClass.get(otuId, "--"),
                                        readCls58S.get(readId, "--"),
                                        otuCls58S.get(otuId, "--"),
                                        combCls.get(otuId, "--")
                                        )
                 )
        try:
            cls58S = otuCls58S[otuId]
            if cls58S == "unknown":
                depth58S = 0
            else:
                depth58S = len(cls58S.strip(";").split(";"))
        except KeyError:
            depth58S = 0
        stat.write("58S\t%s\t%s\t%i\n" % (readId, otuId, depth58S))
        try:
            clsIts = itsClass[otuId]
            if clsIts == "unknown":
                depthIts = 0
            else:
                depthIts = len(clsIts.strip(";").split(";"))
        except KeyError:
            depthIts = 0
        stat.write("ITS2\t%s\t%s\t%i\n" % (readId, otuId, depthIts))
        try:
            clsComb = combCls[otuId]
            if clsComb.strip(";") == "unknown":
                depthComb = 0
            else:
                cLin = clsComb.strip(";").split(";")
                if "|" in cLin[-1]:
                    #conflicted 
                    depthComb = "NA"
                else:
                    depthComb = len(cLin)
        except KeyError:
            depthComb = 0
        stat.write("comb\t%s\t%s\t%s\n" % (readId, otuId, str(depthComb)))

