import random

otutab = []
for line in open(snakemake.input[0]):
    otuId, count = line.strip("\n").split("\t")
    otutab.extend([otuId]*int(count))

points = {}
for r in range(snakemake.params.rep):
    for x in range(1, len(otutab), snakemake.params.step):
        y = len(set(random.sample(otutab, x)))
        try:
            points[x].append(y)
        except KeyError:
            points[x] = [y]

with open(snakemake.output[0], "w") as out:
    for x, yList in points.items():
        for y in yList:
            out.write("%i\t%i\n" % (x,y))
