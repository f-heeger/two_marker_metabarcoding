from Bio import SeqIO

from helpers import *

logOut = open(snakemake.log[0], "w")

classifi = {}
seqLength = {}
seqNr = 0
total = 0
evalueFilter = 0
identFilter = 0
covFilter = 0
tax = {}

for line in open(snakmake.input.tax):
    tId, tLin = line.strip().split("\t")
    tax["%s" % tId] = tLin

for rec in SeqIO.parse(open(snakemake.input.otus), "fasta"):
    seqNr += 1
    classifi[rec.id.split("|", 1)[0]] = []
    seqLength[rec.id.split("|", 1)[0]] = len(rec)

for line in open(snakemake.input.lam, encoding="latin-1"):
    total +=1
    qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.strip().split("\t")
    readId = qseqid.split("|",1)[0]
    if float(evalue) > snakemake.params.maxE:
        evalueFilter += 1
        continue
    if float(pident) < snakemake.params.minIdent:
        identFilter +=1
        continue
    if float(length)/seqLength[readId]*100 < snakemake.params.minCov:
        covFilter += 1
        continue
    linStr = tax[sseqid.split(";", 1)[0]]
    classifi[qseqid.split("|", 1)[0]].append((linStr, float(bitscore)))

logOut.write("%i alignmetns for %i sequences\n" % (total, seqNr))
logOut.write("%i excluded, because e-value was higher than %e\n" % (evalueFilter, snakemake.params.maxE))
logOut.write("%i excluded, because identity was lower than %d%%\n" % (identFilter, snakemake.params.minIdent))
logOut.write("%i excluded, because coverage was lower than %d%%\n" % (covFilter, snakemake.params.minCov))

topPerc = snakemake.params.topPerc/100.0
with open(snakemake.output[0], "w") as out:
    for key, hits in classifi.items():
        if not hits:
            out.write("%s\tunknown\n" % (key))
        else:
            sortedHits = sorted(hits, key=lambda x: x[1])[::-1]
            cutoff = 0
            while cutoff < len(sortedHits) and sortedHits[cutoff][1] >= (1.0-snakemake.topPerc)*sortedHits[0][1]:
                cutoff += 1
            lineage = lca([hit[0] for hit in sortedHits[:cutoff]], snakemake.params.stringency)
            out.write("%s\t%s\n" % (key, lineage))

try:
    logOut.close()
except:
    pass
