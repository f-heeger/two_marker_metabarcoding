otuCount = {}
for line in open(snakemake.input.otuReads):
    otuIdStr, count = line.strip("\n").split("\t")
    otuId = otuIdStr.split("|")[0]
    otuCount[otuId] = int(count)

data = {"$": []}
for line in open(snakemake.input.tax):
    name, linStr = line.strip("\n").split("\t")
    lin = linStr.strip(";").split(";")
    if name not in otuCount:
        #not all OTUs are in all samples
        continue
    if lin[0] == "unknown":
        data["$"].append(otuCount[name])
        continue
    here = data
    for entry in lin:
        try:
            here["$"].append(otuCount[name])
        except KeyError:
            here["$"] = [otuCount[name]]
        try:
            here = here[entry]
        except KeyError:
            here[entry] = {}
            here = here[entry]
    #add count for the last one
    try:
        here["$"].append(otuCount[name])
    except KeyError:
        here["$"] = [otuCount[name]]

with open(snakemake.output[0], "w") as out:
    out.write(json.dumps(data))
