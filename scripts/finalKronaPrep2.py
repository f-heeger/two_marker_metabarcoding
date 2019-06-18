import json
from TaxTreeNode import *

data={}
for inFile in snakemake.input:
    sample = inFile.split("/")[1][:-15]
    data[sample] = json.load(open(inFile))

samples = list(data.keys())
samples.sort()
#create tree
tree = TaxTreeNode("root")
for s in samples:
    tree.loadSubTree(data[s], s)
lines = []
lines.append("<krona>")
lines.append('    <attributes magnitude="reads">')
lines.append('        <attribute display="# Reads">reads</attribute>')
lines.append('        <attribute display="# OTUs">otus</attribute>')
lines.append('    </attributes>')
lines.append('    <datasets>')
for sample in samples:
    lines.append('        <dataset>%s</dataset>' % sample)
lines.append('    </datasets>')
lines.append('    <node name="root">')
lines.extend(tree.kronaXml(samples, "        "))
lines.append('    </node>')
lines.append('</krona>')
with open(snakemake.output[0], "w", encoding="utf-8") as out:
    out.write("\n".join(lines))

