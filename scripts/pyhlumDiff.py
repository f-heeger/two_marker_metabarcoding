with open(snakemake.input[0]) as inStream, open(snakemake.output[0], "w") as out:
    next(inStream) #header
    out.write("otu\tsize\tphylum58s\tphylumIts2\n")
    for line in inStream:
        oId, cls58s, clsits2, clscomb, counts = line.strip("\n").split("\t", 4)
        size = sum([int(c) for c in counts.split("\t")])
        if size < 2:
            #skip singletons
            continue
        if cls58s.split(";")[0] != "k__Fungi" and clsits2.split(";")[0] != "k__Fungi":
            continue
            #skip non-fungi
        try:
            p58s = cls58s.split(";")[1]
        except IndexError:
            p58s = None
        try:
            pits2 = clsits2.split(";")[1]
        except IndexError:
            pits2 = None
        out.write("%s\t%i\t%s\t%s\n" % (oId, size, p58s, pits2))
