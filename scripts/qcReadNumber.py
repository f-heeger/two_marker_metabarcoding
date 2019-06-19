import json

with open(snakemake.output[0], "w") as out, open(snakemake.input[0]) as inStream:
    data = json.load(inStream)
    for sample, sData in data["report_general_stats_data"][0].items():
        sample, lane, readNum, part = sample.rsplit("_", 3)
        if readNum != "R1":
            continue # only use R1
        seqNum = int(sData["total_sequences"])
        out.write("%s\t%i\n" % (sample, seqNum))
