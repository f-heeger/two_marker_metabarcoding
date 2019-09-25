seq2cluster = {}
clusterSize = {}

for line in open(snakemake.input.txt):
    arr = line.strip().split("\t")
    if arr[0] == "C":
        cluster = arr[-2].split(";")[0]
        size = arr[2]
        clusterSize[cluster] = int(size)
    elif arr[0] == "S":
        seq = arr[-2].split(";")[0]
        seq2cluster[seq] = seq
    elif arr[0] == "H":
        seq, cluster = arr[-2:]
        seq2cluster[seq.split(";")[0]] = cluster.split(";")[0]
    else:
        raise ValueError("Unknown record type: %s" % arr[0])
    
with open(snakemake.output.tsv, "w") as out:
    for seq, cluster in seq2cluster.items():
        if clusterSize[cluster] >= snakemake.params.minsize:
            out.write("%s\t%s\n" % (seq, cluster))
