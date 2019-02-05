
rule its_copyIts:
    input: "itsx/all.ITS2.fasta"
    output: "itsx/all.ITS2_extracted.fasta"
    run:
        with open(output[0], "w") as out:
            for rec in SeqIO.parse(open(input[0]), "fasta"):
                rId, sizeStr = rec.id.split("|")[0].strip(";").split(";")
                rec.id = "%s|ITS2;%s;" % (rId, sizeStr)
                out.write(rec.format("fasta"))

rule its_goodReads:
    input: "itsx/all.ITS2_extracted.fasta"
    output: "itsx/all.ITS2_extracted.good.fasta"
    log: "logs/all_good_ITS.log"
    run:
        with open(output[0], "w") as out:
            for rec in SeqIO.parse(open(input[0]), "fasta"):
                if rec.seq.count("N") == 0:
                    out.write(rec.format("fasta"))

rule its_dereplicate:
    input: "itsx/all.ITS2_extracted.good.fasta"
    output: fasta="its_derep/all.ITS2_extracted.derep.fasta", tsv="readInfo/all.repITS2.tsv", txt="its_derep/all.uc.txt"
    log: "logs/all_repIts.log"
    params: minsize=2
    conda:
        "envs/vsearch.yaml"
    run:
        shell("vsearch --derep_fulllength {input} --output {output.fasta} --uc {output.txt} --sizein --sizeout --minuniquesize {params.minsize} --log {log}"
        
        seq2cluster = {}
        clusterSize = {}
        for line in open(output.txt):
            arr = line.strip().split("\t")
            if arr[0] == "C":
                cluster = arr[-2].split(";")[0]
                size = arr[2]
                clusterSize[cluster] = int(size)
            elif arr[0] == "S":
                seq = arr[-2].split(";")[0]
                seq2cluster[seq] = seq
            elif arr[0] == "H":
                seq, cluster = arr[-2:]
                seq2cluster[seq.split(";")[0]] = cluster.split(";")[0]
            else:
                raise ValueError("Unknown record type: %s" % arr[0])
            
        with open(output.tsv, "w") as out:
            for seq, cluster in seq2cluster.items():
                if clusterSize[cluster] >= params.minsize:
                    out.write("%s\t%s\n" % (seq, cluster))

rule its_clustering:
    input: "its_derep/all.ITS2_extracted.derep.fasta"
    output: seeds="swarm/all.ITS2.otus.fasta", otuList="swarm/all.ITS2.otus.out"
    log: "logs/all_swarm.log"
    threads: 3
    conda:
        "envs/swarm2.yaml"
    shell:
        "swarm -f -z -t {threads} -w {output.seeds} {input} -o {output.otuList} &> {log}"


rule its_alignToUnite:
    input: otus="swarm/all.ITS2.otus.fasta", db="%(dbFolder)s/unite_%(unite_version)s.fasta" % config, dbFlag="%(dbFolder)s/unite_%(unite_version)s.fasta.lambdaIndexCreated" % config
    output: "lambda/all.ITS2.otus_vs_UNITE.m8"
    log: "logs/all_lambda.log"
    threads: 3
    conda:
        "envs/lambda.yaml"
    shell:
        "lambda -q {input.otus} -d {input.db} -o {output} -p blastn -t {threads} &> {log}"

rule its_classify:
    input: lam="lambda/all.ITS2.otus_vs_UNITE.m8", otus="swarm/all.ITS2.otus.fasta",  tax="%(dbFolder)s/unite_%(unite_version)s.tsv" % config
    output: "taxonomy/all.ITS2.otus.class.tsv"
    params: maxE=1e-6, topPerc=5.0, minIdent=80.0, minCov=85.0, stringency=.90
    log: "logs/its_class.log"
    run:
        logOut = open(log[0], "w")
        tax = {}
        for line in open(input.tax):
            tId, tLin = line.strip().split("\t")
            tax["%s" % tId] = tLin
        classifi = {}
        itsLength = {}
        seqNr = 0
        total = 0
        evalueFilter = 0
        identFilter = 0
        covFilter = 0
        for rec in SeqIO.parse(open(input.otus), "fasta"):
            seqNr += 1
            classifi[rec.id] = []
            itsLength[rec.id.split("|", 1)[0]] = len(rec)
        for line in open(input.lam, encoding="latin-1"):
            total +=1
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.strip().split("\t")
            readId = qseqid.split("|",1)[0]
            if float(evalue) > params.maxE:
                evalueFilter += 1
                continue
            if float(pident) < params.minIdent:
                identFilter +=1
                continue
            if float(length)/itsLength[readId]*100 < params.minCov:
                covFilter += 1
                continue
            linStr = tax[sseqid]
            classifi[qseqid].append((linStr, float(bitscore)))
        logOut.write("%i alignmetns for %i sequences\n" % (total, seqNr))
        logOut.write("%i excluded, because e-value was higher than %e\n" % (evalueFilter, params.maxE))
        logOut.write("%i excluded, because identity was lower than %d%%\n" % (identFilter, params.minIdent))
        logOut.write("%i excluded, because coverage was lower than %d%%\n" % (covFilter, params.minCov))
        topPerc = params.topPerc/100.0
        with open(output[0], "w") as out:
            for key, hits in classifi.items():
                if not hits:
                    out.write("%s\tunknown\n" % (key))
                else:
                    sortedHits = sorted(hits, key=lambda x: x[1])[::-1]
                    cutoff = 0
                    while cutoff < len(sortedHits) and sortedHits[cutoff][1] >= (1.0-topPerc)*sortedHits[0][1]:
                        cutoff += 1
                    lineage = lca([hit[0] for hit in sortedHits[:cutoff]], params.stringency)
                    out.write("%s\t%s\n" % (key, lineage))
        try:
            logOut.close()
        except:
            pass

rule createOtuReads:
    input: otuList="swarm/all.ITS2.otus.out", repIts="readInfo/all.repITS2.tsv", repSeq="readInfo/all.repseq.tsv"
    output: "swarm/all.otuReads.tsv"
    run:
        repSeq = {}
        for line in open(input.repSeq):
            read, rep = line.strip().split("\t")
            try:
                repSeq[rep].append(read)
            except KeyError:
                repSeq[rep] = [read]
        repIts = {}
        for line in open(input.repIts):
            seq, rep = line.strip().split("\t")
            try:
                repIts[rep.split("|")[0]].append(seq.split("|")[0])
            except KeyError:
                repIts[rep.split("|")[0]] = [seq.split("|")[0]]
        otu = {}
        for line in open(input.otuList):
            memSeqs = line.strip().split(" ")
            otu[memSeqs[0].strip(";").split(";")[0]] = [s.split("|")[0] for s in memSeqs]
        with open(output[0], "w") as otuOut:
            for otuId in otu.keys():
                for itsSeq in otu[otuId]:
                    for repSeqId in repIts[itsSeq]:
                        for read in repSeq[repSeqId]:
                            otuOut.write("%s\t%s\n" % (otuId, read))

rule its_readClassification:
    input: tax="taxonomy/all.ITS2.otus.class.tsv", otuList="swarm/all.ITS2.otus.out", otuReads="swarm/all.otuReads.tsv"
    output: cls="taxonomy/all.ITS2.classification.tsv"
    run:
        otuReads = {}
        for line in open(input.otuReads):
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
        with open(output.cls, "w") as out:
            for read, cls in readClass.items():
                out.write("%s\t%s\n" % (read, cls))

rule its_get58sClassifications:
    input: otuReads="swarm/all.otuReads.tsv", r58SreadCls="taxonomy/all_5_8S_classification.tsv"
    output: "taxonomy/all_ITS2.otus_5.8sClass.tsv"
    log: "logs/all_get58sClass.log"
    params: stringency=0.90
    run:
        with open(log[0], "w") as logFile, open(output[0], "w") as out:
            
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
            
            for otu, reads in otuReads.items():
                otuCls = []
                for readId in reads:
                    try:
                        tCls = readCls58S[readId]
                        if tCls != "unknown":
                            otuCls.append(tCls)
                    except KeyError:
                        logFile.write("No 5.8S classification for read %s.\n" % readId)
                if otuCls:
                    lcaStr =  lca(otuCls, params.stringency, unidentified=["unidentified", "unclassified", "unknown"])
                    out.write("%s\t%s\t%i\n" % (otu, lcaStr, len(otuCls)))
                    otuClsCount = {}
                    for cls in otuCls:
                        try:
                            otuClsCount[cls] += 1
                        except KeyError:
                            otuClsCount[cls] = 1
                    logFile.write("%s\t%s\n" % (otu, otuClsCount))
                else:
                    out.write("%s\tunknown\t0\n" % (otu))

rule its_kronaPrep:
    input: cls="taxonomy/all.ITS2.classification.tsv", sample="readInfo/sample_R1.tsv"
    output: tab=expand("krona/{sample}_ITS2.krona.tsv", sample=samples)
    run:
        sampleReads = {}
        for line in open(input.sample):
            read, sample = line.strip("\n").split("\t")
            try:
                sampleReads[sample].append(read)
            except KeyError:
                sampleReads[sample] = [read]
        readCls = {}
        for line in open(input.cls):
            read, cls = line.strip("\n").split("\t")
            readCls[read] = cls
        for outPath in output:
            sample = outPath.split("/")[1][:-15]
            with open(outPath, "w") as out:
                clsCount = {}
                for read in sampleReads[sample]:
                    if read not in readCls:
                        #this are reads that were filtered out after the sample
                        # file was created
                        continue
                    try:
                        clsCount[readCls[read]] += 1
                    except KeyError:
                        clsCount[readCls[read]] = 1
                for cls, count in clsCount.items():
                    out.write("%i\t%s\n" % (count, "\t".join(cls.split(";"))))

rule its_krona:
    input: expand("krona/{sample}_ITS2.krona.tsv", sample=samples)
    output: "krona/ITS2.krona.html"
    conda:
        "envs/krona.yaml"
    shell:
        "ktImportText -o {output} {input}"


rule its_perSampleOtuReads:
    input: otuReads="swarm/all.otuReads.tsv", sample="readInfo/sample_R1.tsv"
    output: expand("swarm/{sample}.ITS2.otus.out", sample=samples)
    run:
        readOtu = {}
        for line in open(input.otuReads):
            otuId, readId = line.strip("\n").split("\t")
            readOtu[readId] = otuId
        sampleReads = {}
        for line in open(input.sample):
            read, sample = line.strip("\n").split("\t")
            try:
                sampleReads[sample].append(read)
            except KeyError:
                sampleReads[sample] = [read]
        for outPath in output:
            sample = outPath.split("/")[1][:-14]
            otuCount = {}
            for read in sampleReads[sample]:
                if read in readOtu:
                    #some reads in the sample never made it to an OTU because 
                    # they were filtered, but otherwise we do this
                    try:
                        otuCount[readOtu[read]] += 1
                    except KeyError:
                        otuCount[readOtu[read]] = 1
            with open(outPath, "w") as out:
                for otuId, count in otuCount.items():
                    out.write("%s\t%i\n" % (otuId, count))


