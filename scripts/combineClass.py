
if snakemake.config["conflictBehavior"] == "mark":
    ifConf = 0
elif snakemake.config["conflictBehavior"] == "5.8S":
    ifConf = 1
elif snakemake.config["conflictBehavior"] == "ITS":
    ifConf = 2
else:
    raise RuntimeError('unknown conflict behavior setting: "%s". Please change the conflict behavior in the config file to "mark", "ITS" or "5.8S".' % snakemake.config["conflictBehavior"])

tsu = {}
for line in open(snakemake.input.r58sCls):
    name, lin, number = line.strip("\n").split("\t")
    tsu[name.split("|")[0]] = lin

conflict = {}
with open(snakemake.output.otuComb, "w") as out:
    for line in open(snakemake.input.itsCls):
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

with open(snakemake.output.conflict, "w") as out:
    for conf, number in conflict.items():
        out.write("%s | %s\t%i\n" % (conf[0], conf[1], number))
