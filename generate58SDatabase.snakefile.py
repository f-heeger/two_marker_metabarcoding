configfile: "config.json"

rule db_all:
    input: aln="dbs/UNITE.5_8S.aln", tax="dbs/UNITE.5_8S.tax"

rule db_getRfamFile:
    output: "dbs/RF00002.seed.fasta"
    shell:
        "wget -O dbs/RF00002.seed.fasta \"http://rfam.xfam.org/family/RF00002/alignment?acc=RF00002&format=fasta&download=1\""
        
rule db_makeRfamDna:
    input: "dbs/RF00002.seed.fasta"
    output: "dbs/RF00002.seed.dna.fasta"
    run:
        with open(output[0], "w") as out:
            for line in open(input[0]):
                if line[0] == ">":
                    out.write(line)
                else:
                    out.write(line.replace("U", "T"))

rule db_getUniteFile:
    output: "dbs/sh_general_release_s_31.01.2016.fasta"
    shell:
        "cd dbs;"\
        "wget https://unite.ut.ee/sh_files/sh_general_release_s_31.01.2016.zip;"\
        "unzip sh_general_release_s_31.01.2016.zip;"\
        "rm sh_general_release_s_31.01.2016.zip"

rule db_extract58S:
    input: "dbs/sh_general_release_s_31.01.2016.fasta"
    output: "dbs/sh_general_release_s_31.01.2016.5_8S.fasta"
    threads: 6
    shell:
        "%(itsx)s -t . -i {input} -o dbs/sh_general_release_s_31.01.2016 --save_regions 5.8S --cpu {threads} --graphical F" % config
        
rule db_mothur_uniq:
    input: "dbs/sh_general_release_s_31.01.2016.5_8S.fasta"
    output: "dbs/sh_general_release_s_31.01.2016.5_8S.unique.fasta", "dbs/sh_general_release_s_31.01.2016.5_8S.names"
    shell:
        "%(mothur)s \"#unique.seqs(fasta={input});\"" % config
        
rule db_mothur_align:
    input: seqs="dbs/sh_general_release_s_31.01.2016.5_8S.unique.fasta", tmpl="dbs/RF00002.seed.dna.fasta"
    output: "dbs/sh_general_release_s_31.01.2016.5_8S.unique.align", "dbs/sh_general_release_s_31.01.2016.5_8S.unique.align.report", "dbs/sh_general_release_s_31.01.2016.5_8S.unique.flip.accnos"
    shell:
        "%(mothur)s \"#align.seqs(candidate={input.seqs}, template={input.tmpl}, flip=T);\"" % config
        
rule db_mothur_screen:
    input: aln="dbs/sh_general_release_s_31.01.2016.5_8S.unique.align", names="dbs/sh_general_release_s_31.01.2016.5_8S.names"
    output: "dbs/sh_general_release_s_31.01.2016.5_8S.unique.good.align", "dbs/sh_general_release_s_31.01.2016.5_8S.unique.bad.accnos", "dbs/sh_general_release_s_31.01.2016.5_8S.good.names"
    shell:
         "%(mothur)s \"#screen.seqs(fasta={input.aln}, name={input.names}, start=1, end=207, maxn=0, maxambig=0);\"" % config

rule db_create_dbfiles:
    input: aln="dbs/sh_general_release_s_31.01.2016.5_8S.unique.good.align", names="dbs/sh_general_release_s_31.01.2016.5_8S.good.names"
    output: aln="dbs/UNITE.5_8S.aln", tax="dbs/UNITE.5_8S.tax"
    run:
        #get a consensus taxonomy for the entries in my alignment which are 
        # collapsed from identic 5.8S sequences from different species
        taxonomy = {}
        with open(input.names, encoding="latin-1") as nameFile:
            for line in nameFile:
                rId, membersStr = line.split("\t")
                lineage = []
                members = membersStr.split(",")
                taxonomy[rId] = mothurLCA([m.split("|")[4] for m in members])
                    
        with open(output.tax, "w") as tax, open(output.aln, "w") as aln:
            for rec in SeqIO.parse(open(input.aln, encoding="latin-1"), "fasta"):
                newId = "_".join(rec.id.split("|")[:4])
                tax.write("%s\t%s;\n" % (newId, taxonomy[rec.id]))
                rec.description=""
                rec.id = newId
                aln.write(rec.format("fasta"))
                
def mothurLCA(lineageStrings, stringency=1.0):
    lineage = []
    mLineages = []
    #remove bootstrap values ("(100)", "(75)", etc.) if any
    for mLin in [l.strip(";").split(";") for l in lineageStrings]:
        mLineages.append([])
        for entry in mLin:
             mLineages[-1].append(entry.split("(")[0])
    i=0
    minLinLen = min([len(m) for m in mLineages])
    while i<minLinLen:
        tLin = None # <- this will save the first non "unknown" entry in this lineage level
        different = 0.0
        total = 0.0
        for memberLin in mLineages:
            if memberLin[i].split("__")[-1] in ["unidentified", "unclassified", "unknown"]:
                # ignoring unidentified entrys
                continue
            total += 1
            if tLin is None:
                tLin = memberLin[i] # if we have not seen a non "unknown" entry at this level, this becomes our reference
            else:
                if memberLin[i] != tLin:
                    different += 1
        #accept the lineage entry if its proportion of all (valid) classifications higher than stringency setting
        if not tLin is None and different/total <= (1.0-stringency):
            lineage.append(tLin)
        else:
            break
        i += 1
    if len(lineage) == 0:
        lineage = ["unknown"]
    return ";".join(lineage)
