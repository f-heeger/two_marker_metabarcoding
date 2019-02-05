
rule r58S_extract_58S:
    input: pos="itsx/all.positions.txt", seq="chimera/all.nochimera.fasta"
    output: fasta="itsx/all.5_8S_extracted.fasta", gtf="itsx/all.gtf"
    log: "logs/all_5_8Sextraction.log"
    params: minLen=75
    run:
        pos = {}
        for seqId, entry in ITSxParser(open(input.pos)):
            pos[seqId] = entry
        with open(output.fasta, "w") as out, open(output.gtf, "w") as gtf:
            total = 0
            nothingFound = 0
            no58s = 0
            tooShort = 0
            for rec in SeqIO.parse(open(input.seq), "fasta"):
                total += 1
                try:
                    anno = pos[rec.id]
                except KeyError:
                    nothingFound += 1
                    continue
                for key, value in anno.items():
                    if not value is None:
                        start, end = value
                        gtf.write("\t".join([rec.id, "ITSx", key, str(start+1), str(end), ".", "+", ".", "."]) + "\n")
                        if key == "ITS2":
                            if start == 0:
                                no58s += 1
                            elif start < params.minLen:
                                tooShort +=1
                            else:
                                oldId, size = rec.id.strip(";").split(";")
                                rec.id = "%s|5.8S;%s;" % (oldId, size)
                                rec.description="Extracted 5.8S sequence 0-%i" % (start-1)
                                out.write(rec[0:start].format("fasta"))
        propMiss = (nothingFound+no58s+tooShort)/total
        if  propMiss > 0.10:
            print("WARNING: in %f%% of all sequences no 5.8S was found. See %s for details." % (propMiss*100, log[0]))
        with open(log[0], "a") as logStream:
            logStream.write("-------- 5.8S Extraction --------\n")
            logStream.write("5.8S was found in %f%% of the %i sequences.\n" % ((1-propMiss)*100, total))
            logStream.write("ITSx found nothing in %i sequences\n" % nothingFound)
            logStream.write("ITSx found no 5.8S (ITS2 starts at 0) in %i sequences\n" % no58s)
            logStream.write("ITSx found ITS too close to the start in %i sequences\n" % tooShort)


rule r58S_removePrimerAndNs:
    input: "itsx/all.5_8S_extracted.fasta"
    output: "primerremoved/all.5_8S_primerRemoved.fasta"
    log: "logs/all_58S_cutadapt.log"
    params: minOverlap=10, maxN=1
    conda:
        "envs/cutadapt.yaml"
    shell:
        "cutadapt -g ^%(forward_primer)s --trimmed-only -O {params.minOverlap} --max-n={params.maxN} -o {output} {input} &> {log}" % config

rule r58S_dereplicate1:
    """dereplicate 5.8S sequences and only retain "clusters" with more than one sequence"""
    input: "primerremoved/all.5_8S_primerRemoved.fasta"
    output: fasta="r58S_derep/all.5_8S_derep.fasta", txt="r58S_derep/all.uc.txt"
    log: "logs/all_rep58S.log"
    params: minsize=2
    conda:
        "envs/vsearch.yaml"
    shell:
        "vsearch --derep_fulllength {input} --output {output.fasta} --uc {output.txt} --sizein --sizeout --minuniquesize {params.minsize} --log {log}"

rule r58S_dereplicate2:
    input: txt="r58S_derep/all.uc.txt"
    output: tsv="readInfo/all.rep58S.tsv"
    run:
        seq2cluster = {}
        clusterSize = {}
        for line in open(input.txt):
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

rule r58S_align:
    input: otus="r58S_derep/all.5_8S_derep.fasta", db="%(dbFolder)s/58S_derep.fasta" % config, dbFlag="%(dbFolder)s/58S_derep.fasta.lambdaIndexCreated" % config
    output: "lambda/all.58S.derep_vs_58SRef.m8"
    log: "logs/all_58s_lambda.log"
    threads: 3
    conda:
        "envs/lambda.yaml"
    shell:
        "lambda -q {input.otus} -d {input.db} -o {output} -p blastn -t {threads} &> {log}"

rule r58S_classify:
    input: lam="lambda/all.58S.derep_vs_58SRef.m8", otus="r58S_derep/all.5_8S_derep.fasta", tax="%(dbFolder)s/58S_tax.tsv" % config
    output: "taxonomy/all.58S.derep.class.tsv"
    params: maxE=1e-6, topPerc=5.0, minIdent=80.0, minCov=85.0, stringency=.90
    log: "logs/58s_class.log"
    run:
        logOut = open(log[0], "w")
        classifi = {}
        seqLength = {}
        seqNr = 0
        total = 0
        evalueFilter = 0
        identFilter = 0
        covFilter = 0
        tax = {}
        for line in open(input.tax):
            tId, tLin = line.strip().split("\t")
            tax["%s" % tId] = tLin
        for rec in SeqIO.parse(open(input.otus), "fasta"):
            seqNr += 1
            classifi[rec.id.split("|", 1)[0]] = []
            seqLength[rec.id.split("|", 1)[0]] = len(rec)
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
            if float(length)/seqLength[readId]*100 < params.minCov:
                covFilter += 1
                continue
            linStr = tax[sseqid.split(";", 1)[0]]
            classifi[qseqid.split("|", 1)[0]].append((linStr, float(bitscore)))
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

rule r58S_readClassification:
    input: tax="taxonomy/all.58S.derep.class.tsv", repseq="readInfo/all.repseq.tsv", rep58s="readInfo/all.rep58S.tsv", sample="readInfo/sample_R1.tsv"
    output: "taxonomy/all_5_8S_classification.tsv"
    run:
        repseq = {}
        for line in open(input.repseq):
            read, rep = line.strip().split("\t")
            try:
                repseq[rep].append(read)
            except KeyError:
                repseq[rep] = [read]
        rep58s = {}
        for line in open(input.rep58s):
            seq, rep = line.strip().split("\t")
            try:
                rep58s[rep].append(seq)
            except KeyError:
                rep58s[rep] = [seq]
            
        readClass = {}
        for line in open(input.tax):
            rId, classification = line.strip().split("\t")
            cls = []
            for entry in classification.strip(";").split(";"):
                if entry == "unclassified":
                    break
                elif entry == "unknown":
                    cls = [""]
                else:
                    cls.append(entry)
            for r58seq in rep58s["%s|5.8S" % rId]:
                for read in repseq[r58seq.split("|")[0]]:
                    readClass[read] = ";".join(cls)
        with open(output[0], "w") as out:
            for read, cls in readClass.items():
                out.write("%s\t%s\n" % (read, cls))

rule r58S_prepKronaInput:
    input: cls="taxonomy/all_5_8S_classification.tsv", sample="readInfo/sample_R1.tsv"
    output: expand("krona/{sample}_5_8S.tsv", sample=samples)
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
            sample = outPath.split("/")[1][:-9]
            with open(outPath, "w") as out:
                clsCount = {}
                for read in sampleReads[sample]:
                    try:
                        cls = readCls[read]
                    except KeyError:
                        #reads that were removed in between
                        continue
                    try:
                        clsCount[cls] += 1
                    except KeyError:
                        clsCount[cls] = 1
                for cls, count in clsCount.items():
                    out.write("%i\t%s\n" % (count, "\t".join(cls.split(";"))))

rule r58S_krona:
    input: expand("krona/{sample}_5_8S.tsv", sample=samples)
    output: "krona/5_8s.krona.html"
    conda:
        "envs/krona.yaml"
    shell:
        "ktImportText -o {output} {input}"

