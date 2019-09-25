from helpers import *

clu={}

for line in open(snakemake.input.uc):
    lType, cNum, sLen, ident, strand, _1, _2, aln, query, target = line.strip("\n").split("\t")
    if lType == "S":
        clu[query] = [query]
    elif lType == "H":
        clu[target].append(query)
    elif lType == "C":
        assert int(sLen) == len(clu[query])
    else:
        raise ValueError("Unknown record type: %s" % lType)

uTax = {}
for line in open(snakemake.input.uTax):
    tId, tLin = line.strip().split("\t")
    uTax[tId] = tLin

rTax = {}
for line in open(snakemake.input.rTax):
    tId, tLin = line.strip().split("\t")
    rTax[tId] = tLin

with open(snakemake.output[0], "w") as out:
    for rep, memList in clu.items():
        linStrs = []
        for m in memList:
            if m[:2] == "SH":
                linStrs.append(uTax[m.split("|")[0]])
            else:
                lin = rTax[m]
                arr = lin.split(";")
                if arr[1] == "k__Fungi":
                    #ignore everything RFAM has to say about fungi taxonomy
                    linStrs.append("k__Fungi")
                else:
                    linStrs.append(";".join(arr[1:]))
        lcaLin = lca(linStrs, 0.95)
        out.write("%s\t%s\n" % (rep, lcaLin))
