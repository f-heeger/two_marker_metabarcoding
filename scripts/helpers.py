
def lca(lineageStrings, stringency=1.0, 
        unidentified=["unidentified", "unclassified", "unknown"],
        ignoreIncertaeSedis=True):
    lineage = []
    mLineages = []
    #remove bootstrap values ("(100)", "(75)", etc.) if any
    for mLin in [l.strip(";").split(";") for l in lineageStrings]:
        mLineages.append([])
        for entry in mLin:
             mLineages[-1].append(entry.split("(")[0])
    i=0
    maxLinLen = max([len(m) for m in mLineages])
    active = [True]*len(mLineages)
    for i in range(maxLinLen):
        total = 0.0
        counts = {}
        for m, memberLin in enumerate(mLineages):
            if not active[m]:
                continue #ignore lineages that were deactivated further up in the tree
            if len(memberLin) <= i:
                active[m] = False
                continue #ignore lineages that are not this long
            name = memberLin[i].split("__")[-1]
            if name in unidentified:
                continue # ignoring unidentified entrys
            if ignoreIncertaeSedis and name.startswith("Incertae"):
                continue # ignoring Incertae sedis entries.
                         # NOTE: this will mean lineages end at the first Incerta sedis
            total += 1
            try:
                counts[memberLin[i]] += 1
            except KeyError:
                counts[memberLin[i]] = 1
        if not counts:
            #no valid lineage entrys found in this level
            break
        most=sorted(counts.items(), key=lambda x: x[1], reverse=True)[0]
        different = total - most[1]
        #accept the lineage entry if its proportion of all (valid) classifications higher than stringency setting
        if different/total <= (1.0-stringency):
            lineage.append(most[0]) #add the most apearing entry to the new lineage
            #deactivate all lineages that were different at this level
            for m, memberLin in enumerate(mLineages):
                if active[m] and memberLin[i] != most[0]:
                    active[m] = False
        else:
            break
    if len(lineage) == 0:
        lineage = ["unknown"]
    return ";".join(lineage)

def ITSxParser(inStream):
    for line in inStream:
        entry = {}
        seqId, sedLen, ssuStr, its1Str, tsuStr, its2Str, lsuStr, status = line.strip("\n").split("\t")
        names = ["SSU", "ITS1", "5.8S", "ITS2", "LSU"]
        data = [ssuStr, its1Str, tsuStr, its2Str, lsuStr]
        for name, dataStr in zip(names, data):
            entry[name] = readKeyValuePair(dataStr)
        yield seqId, entry
        
def readKeyValuePair(inString):
    key, value = inString.split(": ")
    if value == "Not found" or value == "No start" or value == "No end":
        return None
    else:
        try:
            start, end = value.split("-")
        except ValueError:
            print(key, value)
            raise
        return int(start)-1, int(end)
