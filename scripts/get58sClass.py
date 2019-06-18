from helpers import *

with open(snakemake.log[0], "w") as logFile, open(snakemake.output[0], "w") as out:
    
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

    for otu, reads in otuReads.items():
        otuCls = []
        for readId in reads:
            try:
                tCls = readCls58S[readId]
                if tCls != "unknown":
                    otuCls.append(tCls)
            except KeyError:
                logFile.write("No 5.8S classification for read %s.\n" % readId)
        if otuCls:
            lcaStr =  lca(otuCls, snakemake.params.stringency, unidentified=["unidentified", "unclassified", "unknown"])
            out.write("%s\t%s\t%i\n" % (otu, lcaStr, len(otuCls)))
            otuClsCount = {}
            for cls in otuCls:
                try:
                    otuClsCount[cls] += 1
                except KeyError:
                    otuClsCount[cls] = 1
            logFile.write("%s\t%s\n" % (otu, otuClsCount))
        else:
            out.write("%s\tunknown\t0\n" % (otu))
