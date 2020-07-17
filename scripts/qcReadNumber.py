import json

file2sample = {}
for sample, fileInfo in snakemake.params.files.items():
    for fileName in fileInfo[0]:
        file2sample[fileName.rsplit(".", 2)[0]] = sample

seqNum = {}
with open(snakemake.input[0]) as inStream:
    data = json.load(inStream)
    for fileName, sData in data["report_general_stats_data"][0].items():
        libName, lane, readNum, part = fileName.rsplit("_", 3)
        if readNum != "R1":
            continue # only use R1
        n = int(sData["total_sequences"])
        sample = file2sample[fileName]
        try:
            seqNum[sample] += n
        except KeyError:
            seqNum[sample] = n

with open(snakemake.output[0], "w") as out:
    for sample, n in seqNum.items():
        out.write("%s\t%i\n" % (sample, n))

