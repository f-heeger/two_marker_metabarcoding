depth = {}
size = {}

for line in open(snakemake.input[0]):
    marker, read, otu, tDepth = line.strip().split("\t")
    try:
        depth[otu][marker] = tDepth
    except KeyError:
        depth[otu] = {"ITS2": None, "58S": None}
        depth[otu][marker] = tDepth
        size[otu] = 0
    if marker == "ITS2":
        size[otu] += 1

with open(snakemake.output[0], "w") as out:
    for oId, mDepth in depth.items():
        out.write("%s\t%i\t%s\t%s\n" % (oId, size[oId], mDepth["58S"], mDepth["ITS2"]))
