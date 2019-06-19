from Bio import SeqIO

from helpers import *

pos = {}

for seqId, entry in ITSxParser(open(snakemake.input.pos)):
    pos[seqId] = entry

with open(snakemake.output.fasta, "w") as out, open(snakemake.output.gtf, "w") as gtf:
    total = 0
    nothingFound = 0
    no58s = 0
    tooShort = 0
    for rec in SeqIO.parse(open(snakemake.input.seq), "fasta"):
        total += 1
        try:
            anno = pos[rec.id]
        except KeyError:
            nothingFound += 1
            continue
        for key, value in anno.items():
            if not value is None:
                start, end = value
                gtf.write("\t".join([rec.id, "ITSx", key, str(start+1), str(end), ".", "+", ".", "."]) + "\n")
                if key == "ITS2":
                    if start == 0:
                        no58s += 1
                    elif start < snakemake.params.minLen:
                        tooShort +=1
                    else:
                        oldId, size = rec.id.strip(";").split(";")
                        rec.id = "%s|5.8S;%s;" % (oldId, size)
                        rec.description="Extracted 5.8S sequence 0-%i" % (start-1)
                        out.write(rec[0:start].format("fasta"))
propMiss = (nothingFound+no58s+tooShort)/total
if  propMiss > 0.10:
    print("WARNING: in %f%% of all sequences no 5.8S was found. See %s for details." % (propMiss*100, snakmake.log[0]))
with open(snakemake.log[0], "a") as logStream:
    logStream.write("-------- 5.8S Extraction --------\n")
    logStream.write("5.8S was found in %f%% of the %i sequences.\n" % ((1-propMiss)*100, total))
    logStream.write("ITSx found nothing in %i sequences\n" % nothingFound)
    logStream.write("ITSx found no 5.8S (ITS2 starts at 0) in %i sequences\n" % no58s)
    logStream.write("ITSx found ITS too close to the start in %i sequences\n" % tooShort)

