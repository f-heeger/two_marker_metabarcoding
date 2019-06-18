otuReads = {}
for line in open(snakmake.input.otuReads):
    otuId, readId = line.strip().split("\t")
    try:
        otuReads[otuId].append(readId)
    except KeyError:
        otuReads[otuId] = [readId]

readClass = {}
for line in open(input.tax):
    otuName, classification = line.strip().split("\t")
    otuId, sizeStr = otuName.strip(";").split(";")
    count = int(sizeStr.split("=")[1])
    assert len(otuReads[otuId]) == count
    cls = []
    for entry in classification.strip(";").split(";"):
        if entry == "unclassified":
            break
        elif entry == "unknown":
            cls = [""]
        else:
            cls.append(entry)
    for read in otuReads[otuId]:
        readClass[read] = ";".join(cls)

with open(snakmake.output.cls, "w") as out:
    for read, cls in readClass.items():
        out.write("%s\t%s\n" % (read, cls))
