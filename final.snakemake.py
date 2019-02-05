import random
import json

rule final_concatMarkers:
    input: otus="swarm/all.ITS2.otus.fasta", tsu="primerremoved/all.5_8S_primerRemoved.fasta"
    output: "all.concatMarkers.fasta"
    log: "logs/all.concatMarker.log"
    run:
        r58s = {}
        for r in SeqIO.parse(open(input.tsu), "fasta"):
            r58s[r.id.split("|", 1)[0]] = r
        notFound = 0
        with open(output[0], "w") as out, open(log[0], "w") as logStream:
            for otu in SeqIO.parse(open(input.otus), "fasta"):
                newId = otu.id.split("|", 1)[0]
                try:
                    tsu = r58s[newId]
                except KeyError:
                    logStream.write("Did not find 5.8S for %s\n" % newId)
                    notFound += 1
                concat = tsu + otu
                concat.id = "%s|partial5.8S+ITS2" % newId
                concat.description = ""
                out.write(concat.format("fasta"))
        print("Could not find 5.8S for %i OTUs. See %s for details." % (notFound, log[0]))

rule final_combineClassification:
    input: itsCls="taxonomy/all.ITS2.otus.class.tsv", r58sCls="taxonomy/all_ITS2.otus_5.8sClass.tsv"
    output: otuComb="taxonomy/all.ITS2.otus.combClass.tsv", conflict="taxonomy/all.ITS2.otus.conflictingClass.tsv"
    run:
        if config["conflictBehavior"] == "mark":
            ifConf = 0
        elif config["conflictBehavior"] == "5.8S":
            ifConf = 1
        elif config["conflictBehavior"] == "ITS":
            ifConf = 2
        else:
            raise RuntimeError('unknown conflict behavior setting: "%s". Please change the conflict behavior in the config file to "mark", "ITS" or "5.8S".' % config["conflictBehavior"])

        tsu = {}
        for line in open(input.r58sCls):
            name, lin, number = line.strip("\n").split("\t")
            tsu[name.split("|")[0]] = lin
        conflict = {}
        with open(output.otuComb, "w") as out:
            for line in open(input.itsCls):
                nameStr, linStr = line.strip("\n").split("\t")
                name = nameStr.split("|")[0]
                lin = linStr.split(";")
                tsuLin = tsu[name].split(";")[:2] # trim down to phylum
                if len(lin) <= 1:
                    #if the ITS classification is only one level deep 
                    # (this could also be "unknown") compare it to the first level
                    # of the 5.8S classification
                    if lin[0] == tsuLin[0] or (lin[0] == "unknown" and tsuLin[0] != "unknown"):
                        #accept 5.8S classification if it is the same as the ITS or 
                        # if the ITS is unknown and the 5.8S is not
                        out.write("%s\t%s;\n" % (name, ";".join(tsuLin)))
                    else:
                        #if this is not the cas this is a conflict at the first level:
                        # act accordingly (nothing else has to be done)
                        if ifConf==0:
                            out.write("%s\t%s|%s;\n" % (name, lin[0], tsuLin[0]))
                        elif ifConf==1:
                            out.write("%s\t%s;\n" % (name, ";".join(tsuLin)))
                        else:
                            out.write("%s\t%s;\n" % (name, ";".join(lin)))
                        try:
                            conflict[(lin[0], tsuLin[0])] += 1
                        except:
                            conflict[(lin[0], tsuLin[0])] = 1
                elif len(tsuLin) == 1:
                    #if the 5.8S classifaction is shorter than 2
                    # (and the ITS classification is longer than 1;
                    # otherwise the first case would have acted)
                    if lin[0] == tsuLin[0]:
                        #if they are the same accept the ITS classification
                        # (including all levels after the first)
                        out.write("%s\t%s;\n" % (name, ";".join(lin)))
                    else:
                        #otherwise write the the conflict and stop there
                        if ifConf==0:
                            out.write("%s\t%s|%s;\n" % (name, lin[0], tsuLin[0]))
                        elif ifConf==1:
                            out.write("%s\t%s;\n" % (name, ";".join(tsuLin)))
                        else:
                            out.write("%s\t%s;\n" % (name, ";".join(lin)))
                        try:
                            conflict[(lin[0], tsuLin[0])] += 1
                        except:
                            conflict[(lin[0], tsuLin[0])] = 1
                else:
                    #if the ITS classification is more than one level deep
                    # check if the first three levels of classifications are identical
                    if lin[:2] == tsuLin[:2]:
                        #if they are the same accept the ITS classification
                        # (including all levels after the third)
                        out.write("%s\t%s;\n" % (name, linStr))
                    else:
                        #otherwise write the classifcation
                        # if a conflict occurs and the 'mark' conflict behavior is used
                        # mark the conflict and truncate the classification on that level
                        # This is writen for a general case, but the conflict can 
                        # only occut in the first two levels, because we ignore the
                        # 5.8S after that.
                        comLin = []
                        for a,b in zip(lin, tsuLin):
                            if a == b:
                                comLin.append(a)
                            else:
                                if ifConf==0:
                                    comLin.append("%s|%s" % (a,b))
                                    break
                                elif ifConf==1:
                                    comLin.append(b)
                                else:
                                    comLin.append(a)
                        out.write("%s\t%s;\n" % (name, ";".join(comLin)))
                        try:
                            conflict[(";".join(lin[:2]), ";".join(tsuLin[:2]))] += 1
                        except:
                            conflict[(";".join(lin[:2]), ";".join(tsuLin[:2]))] = 1
        with open(output.conflict, "w") as out:
            for conf, number in conflict.items():
                out.write("%s | %s\t%i\n" % (conf[0], conf[1], number))

rule final_classCompare:
    input: itsCls="taxonomy/all.ITS2.otus.class.tsv", otuReads="swarm/all.otuReads.tsv", r58SreadCls="taxonomy/all_5_8S_classification.tsv", r58sOtuCls="taxonomy/all_ITS2.otus_5.8sClass.tsv", combCls="taxonomy/all.ITS2.otus.combClass.tsv", sample="readInfo/sample_R1.tsv"
    output: comp="taxonomy/all.compareClass.tsv", stat="taxonomy/all.clsStat.tsv"
    run:
        readSample = {}
        for line in open(input.sample):
            read, sample = line.strip("\n").split("\t")
            readSample[read] = sample
        otuReads = {}
        for line in open(input.otuReads):
            otuId, readId = line.strip().split("\t")
            try:
                otuReads[otuId].append(readId)
            except KeyError:
                otuReads[otuId] = [readId]
        readCls58S = {}
        for line in open(input.r58SreadCls):
            read, cls = line.strip("\n").split("\t")
            if cls == "":
                readCls58S[read] = "unknown"
            else:
                readCls58S[read] = cls
        otuCls58S = {}
        for line in open(input.r58sOtuCls):
            name, lin, number = line.strip("\n").split("\t")
            otuCls58S[name.split("|")[0]] = lin
        itsClass = {}
        for line in open(input.itsCls):
            nameStr, cls = line.strip("\n").split("\t")
            itsClass[nameStr.split("|")[0]] = cls
        combCls = {}
        for line in open(input.combCls):
            otuId, cls = line.strip("\n").split("\t")
            combCls[otuId] = cls
        with open(output.comp, "w") as out, open(output.stat, "w") as stat:
            out.write("read\tsample\tOTU\tITS OTU Class\t5.8S read Class\t5.8S OTU Class\tComb OTU Class\n")
            for line in open(input.otuReads):
                otuIdStr, readId = line.strip("\n").split("\t")
                otuId = otuIdStr.split("|")[0]
                out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (readId,
                                                readSample[readId],
                                                otuId,
                                                itsClass.get(otuId, "--"),
                                                readCls58S.get(readId, "--"),
                                                otuCls58S.get(otuId, "--"),
                                                combCls.get(otuId, "--")
                                                )
                         )
                try:
                    cls58S = otuCls58S[otuId]
                    if cls58S == "unknown":
                        depth58S = 0
                    else:
                        depth58S = len(cls58S.strip(";").split(";"))
                except KeyError:
                    depth58S = 0
                stat.write("58S\t%s\t%s\t%i\n" % (readId, otuId, depth58S))
                try:
                    clsIts = itsClass[otuId]
                    if clsIts == "unknown":
                        depthIts = 0
                    else:
                        depthIts = len(clsIts.strip(";").split(";"))
                except KeyError:
                    depthIts = 0
                stat.write("ITS2\t%s\t%s\t%i\n" % (readId, otuId, depthIts))
                try:
                    clsComb = combCls[otuId]
                    if clsComb.strip(";") == "unknown":
                        depthComb = 0
                    else:
                        cLin = clsComb.strip(";").split(";")
                        if "|" in cLin[-1]:
                            #conflicted 
                            depthComb = "NA"
                        else:
                            depthComb = len(cLin)
                except KeyError:
                    depthComb = 0
                stat.write("comb\t%s\t%s\t%s\n" % (readId, otuId, str(depthComb)))

rule final_diffClsDepth:
    input: "taxonomy/all.clsStat.tsv"
    output: "taxonomy/all.otuClsDepth.tsv"
    run:
        depth = {}
        size = {}
        for line in open(input[0]):
            marker, read, otu, tDepth = line.strip().split("\t")
            try:
                depth[otu][marker] = tDepth
            except KeyError:
                depth[otu] = {"ITS2": None, "58S": None}
                depth[otu][marker] = tDepth
                size[otu] = 0
            if marker == "ITS2":
                size[otu] += 1
        
        with open(output[0], "w") as out:
            for oId, mDepth in depth.items():
                out.write("%s\t%i\t%s\t%s\n" % (oId, size[oId], mDepth["58S"], mDepth["ITS2"]))

rule final_phylumDiff:
    input: "otu_table.tsv"
    output: "fungiPhylumDiff.tsv"
    run:
        with open(input[0]) as inStream, open(output[0], "w") as out:
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

rule final_plotPhylumDiff:
    input: "fungiPhylumDiff.tsv"
    output: "gainFrom58s.pdf"
    conda:
        "envs/ggplot.yaml"
    script:
        "scripts/plotPhylumDiff.R"

rule final_computeRarefaction:
    input: "swarm/{sample}.ITS2.otus.out"
    output: "rarefaction/{sample}.rarefaction.tsv"
    params: rep=100, step=1000
    run:
        otutab = []
        for line in open(input[0]):
            otuId, count = line.strip("\n").split("\t")
            otutab.extend([otuId]*int(count))
        points = {}
        for r in range(params.rep):
            for x in range(1, len(otutab), params.step):
                y = len(set(random.sample(otutab, x)))
                try:
                    points[x].append(y)
                except KeyError:
                    points[x] = [y]
        with open(output[0], "w") as out:
            for x, yList in points.items():
                for y in yList:
                    out.write("%i\t%i\n" % (x,y))

rule final_combineRarefactions:
    input: expand("rarefaction/{sample}.rarefaction.tsv", sample=samples)
    output: "rarefaction/All.rarefactions.tsv"
    run:
        with open(output[0], "w") as out:
            for sample in samples:
                for line in open("rarefaction/%s.rarefaction.tsv" % sample):
                    out.write("%s\t%s" % (sample, line)) #no \n at the end as it was not striped from the line

rule final_plotRarefactions:
    input: "rarefaction/All.rarefactions.tsv"
    output: "All.rarefactions.pdf"
    conda:
        "envs/ggplot.yaml"
    script:
        "scripts/plotRarefaction.R"

rule final_creatOtuTable:
    input: tax58s="taxonomy/all_ITS2.otus_5.8sClass.tsv", taxIts="taxonomy/all.ITS2.otus.class.tsv", taxComb="taxonomy/all.ITS2.otus.combClass.tsv", otus=expand("swarm/{sample}.ITS2.otus.out", sample=samples)
    output: "otu_table.tsv"
    run:
        readNr = {}
        for sample in samples:
            readNr[sample] = {}
            for line in open("swarm/%s.ITS2.otus.out" % sample):
                otuStr, nr = line.strip("\n").split("\t")
                otu = otuStr.split("|")[0]
                readNr[sample][otu] = nr
        r58sCls = {}
        for line in open(input.tax58s):
            name, lin, nr = line.strip().split("\t")
            r58sCls[name.split("|")[0]] = lin
        itsCls = {}
        for line in open(input.taxIts):
            name, lin = line.strip().split("\t")
            itsCls[name.split("|")[0]] = lin
        with open(output[0], "w") as out:
            out.write("otu_ID\t5.8S classification\tITS2 classification\tfinal classification\t%s\n" % "\t".join(samples))
            for line in open(input.taxComb):
                otu, cls = line.strip("\n").split("\t")
                numbers = [readNr[sample].get(otu, "0") for sample in samples]
                out.write("%s\t%s\t%s\t%s\t%s\n" % (otu, r58sCls[otu], itsCls[otu], cls, "\t".join(numbers)))

rule final_kronaPrepStep1:
    input: tax="taxonomy/all.ITS2.otus.combClass.tsv", otuReads = "swarm/{sample}.ITS2.otus.out"
    output: "krona/{sample}_all.krona.json"
    run:
        otuCount = {}
        for line in open(input.otuReads):
            otuIdStr, count = line.strip("\n").split("\t")
            otuId = otuIdStr.split("|")[0]
            otuCount[otuId] = int(count)
        data = {"$": []}
        for line in open(input.tax):
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
        
        with open(output[0], "w") as out:
            out.write(json.dumps(data))


rule final_kronaPrepStep2:
    input:  expand("krona/{sample}_all.krona.json", sample=samples)
    output: "krona/All.krona.xml"
    run:
        data={}
        for inFile in input:
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
        with open(output[0], "w", encoding="utf-8") as out:
            out.write("\n".join(lines))

rule final_kronaAll:
    input: "krona/All.krona.xml"
    output: "krona/All.krona.html"
    conda:
        "envs/krona.yaml"
    shell:
        "ktImportXML -o {output} {input}"
