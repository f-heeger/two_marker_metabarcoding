import gzip

from NcbiMap import NuclId2TaxIdMap, LineageMap

ranks = ["superkingdom","kingdom", "phylum", "class", "order", "family", "genus", "species"]

config = snakemake.config

if not config["email"] or not "@" in config["email"]:
    raise ValueError("'%s' is not a valid email. Set email in config file for NCBI querries." % config["email"])

nuc2tax = NuclId2TaxIdMap(config["email"], cachePath="%(dbFolder)s/nuc2tax.csv" % config, retry=3)
tax2lin = LineageMap(config["email"], cachePath="%(dbFolder)s/tax2lin.csv" % config, retry=3)

noId = 0
noLin = 0

with open(snakemake.output.fasta, "w") as out, \
     open(snakmake.output.tax, "w") as taxOut, \
     open(snakmake.log[0], "w") as logFile:
    for rec in SeqIO.parse(gzip.open(input.fasta, "rt"), "fasta"):
        
        
        sId, reg = rec.id.split("/")
        rec.id="%s|%s" % (sId, reg)
        try:
            tId = nuc2tax[sId]
        except:
            logFile.write("no ncbi tax ID\t%s\n" % sId)
            noId += 1
            continue
        try:
            lin = tax2lin[tId]
        except:
            logFile.write("no ncbi lineage\t%s\n" % Id)
            noLin += 1
            continue
        linStr = ";".join(["%s__%s" %(l[0][0], l[2]) for l in lin if l[0] in ranks])
        taxOut.write("%s\t%s\n" % (rec.id, linStr))
        out.write(rec.format("fasta"))

sys.stderr.write("No taxon ID: %i, no lineage: %i\n" % (noId, noLin))
nuc2tax.save()
tax2lin.save()
