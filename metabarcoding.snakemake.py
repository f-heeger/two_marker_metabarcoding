import itertools
import json
import xml.etree.ElementTree as et
import xml.dom.minidom as minidom

from snakemake.utils import min_version
from Bio import SeqIO

min_version("3.5.4")

configfile: "config.json"

shell.prefix("sleep 60; ") #work aorund to desl with "too quck" rule execution and slow NAS

samples = ["%s%s_S%i" % (a,b,c) for ((b, a), c) in zip(itertools.product(range(1,13), "ABCDEFGH"), range(1,97))]
#samples =["A1_S1", "B1_S2", "H5_S40", "A5_S33"]


rule all:
    input: "krona/All.krona.html", "krona/5_8s.krona.html", "krona/ITS2.krona.html"
#    input: "taxonomy/All.krona.html", "taxonomy/5_8s.krona.html", expand(["taxonomy/{sample}_ITS2.otus_5.8sClass.tsv", "taxonomy/{sample}.ITS2.otus.class.tsv"], sample=samples)

################ generate reference sequence for 5.8S ##########################

include: "generate58SDatabase.snakefile.py"

################ generate lamda db from UNTIE ##################################

rule db_creatUniteIndex:
    input: "dbs/sh_general_release_s_31.01.2016.fasta"
    output: touch("dbs/sh_general_release_s_31.01.2016.fasta.lambdaIndexCreated")
    threads: 6
    shell:
        "%(lambdaFolder)s/lambda_indexer -d {input} -p blastn -t {threads}" % config

################ quality control ###############################################

rule qc_concat:
    params: inFolder="/home/heeger/spree/raw_data/2016/illumina/160202_2015-37-CW-Japan-Repeat/Data/Intensities/BaseCalls"
    output: "QC/all_R{read_number}.fastq.gz"
    shell:
        "cat {params.inFolder}/*_L001_R{wildcards.read_number}_001.fastq.gz > {output}"

rule qc_fastqc:
    input: "QC/all_R{read_number}.fastq.gz"
    output: "QC/all_R{read_number}.html"
    threads: 6
    shell:
        "%(fastqc)s --nogroup -o QC --threads {threads} {input}" % config

################ initial sequence processing ###################################

rule init_separatePrimer:
    input: read1="/home/heeger/spree/raw_data/2016/illumina/160202_2015-37-CW-Japan-Repeat/Data/Intensities/BaseCalls/{sample}_L001_R1_001.fastq.gz", read2="/home/heeger/spree/raw_data/2016/illumina/160202_2015-37-CW-Japan-Repeat/Data/Intensities/BaseCalls/{sample}_L001_R2_001.fastq.gz"
    output: "demultiplexed/{sample}_barcode_ITS_1.fastq.gz", "demultiplexed/{sample}_barcode_ITS_2.fastq.gz", "demultiplexed/{sample}_barcode_other_1.fastq.gz", "demultiplexed/{sample}_barcode_other_2.fastq.gz"
    log: "logs/{sample}_flexbar.log"
    shell:
        "echo \">ITS\nCGATGAAGAACG\n>other\nACAAACT\" > barcodes/barcode_{wildcards.sample}.fasta; %(flexbar)s -n 4 -r {input.read1} -p {input.read2} -t demultiplexed/{wildcards.sample} -b barcodes/barcode_{wildcards.sample}.fasta -bk -be LEFT_TAIL -bn 20 -bt 1 -f i1.8 -u 100 -z GZ &> {log}" % config

rule init_trimming:
    input: r1="demultiplexed/{sample}_barcode_ITS_1.fastq.gz", r2="demultiplexed/{sample}_barcode_ITS_2.fastq.gz"
    output: r1="trimmed/{sample}_trimmed_R1.fastq.gz", r2="trimmed/{sample}_trimmed_R2.fastq.gz"
    threads: 3
    params: windowLen=8, minQual=20, minLen=200, avgQual=30
    shell:
        "java -jar %(trimmomatic)s PE -threads {threads} -phred33 {input.r1} {input.r2} {output.r1} {output.r1}.unpaired {output.r2} {output.r2}.unpaired SLIDINGWINDOW:{params.windowLen}:{params.minQual} TRAILING:{params.minQual} MINLEN:{params.minLen} AVGQUAL:{params.avgQual}" % config

rule init_merge:
    input: r1="trimmed/{sample}_trimmed_R1.fastq.gz", r2="trimmed/{sample}_trimmed_R2.fastq.gz"
    output: "merged/{sample}.assembled.fastq"
    params: minLen=300, maxLen=550, minOverlap=10
    threads: 3
    log: "logs/{sample}_pear.log"
    shell:
        "%(pear)s -j {threads} -f {input.r1} -r {input.r2} -o merged/{wildcards.sample} -n {params.minLen} -m {params.maxLen} -v {params.minOverlap} &> {log}" % config

rule init_convertMerged:
    input: "merged/{sample}.assembled.fastq"
    output: "merged/{sample}.fasta"
    run:
        with open(output[0], "w") as out:
            for record in SeqIO.parse(open(input[0]), "fastq"):
                out.write(record.format("fasta"))

rule init_dereplicate:
    input: "merged/{sample}.fasta"
    output: "merged/{sample}.unique.fasta", "merged/{sample}.names"
    log: "logs/{sample}_mothur.log"
    shell:
        "%(mothur)s -q \"#set.logfile(name={log}, append=T); unique.seqs(fasta={input});\" > /dev/null" % config

rule init_cpMergerd:
    input:  seq="merged/{sample}.unique.fasta", name="merged/{sample}.names"
    output:  seq="chimera/{sample}.unique.fasta", name="chimera/{sample}.names"
    shell:
        "cp {input.seq} {output.seq}; cp {input.name} {output.name}"

rule init_removeChimeras:
    input: seq="chimera/{sample}.unique.fasta", name="chimera/{sample}.names"
    output: chimeras="chimera/{sample}.unique.uchime.accnos", fasta="chimera/{sample}.unique.pick.fasta", name="chimera/{sample}.pick.names"
    log: "logs/{sample}_mothur.log"
    shell: 
        "%(mothur)s -q \"#set.logfile(name={log}, append=T); chimera.uchime(fasta={input.seq}, reference=self, name={input.name}); remove.seqs(fasta={input.seq}, accnos={output.chimeras}, name={input.name});\" 2> /dev/null" % config

rule init_undereplicate:
    input: fasta="chimera/{sample}.unique.pick.fasta", name="chimera/{sample}.pick.names"
    output: fasta="chimera/{sample}.unique.pick.redundant.fasta"
    log: "logs/{sample}_mothur.log"
    shell: 
        "%(mothur)s -q \"#set.logfile(name={log}, append=T); deunique.seqs(fasta={input.fasta}, name={input.name});\" > /dev/null" % config

rule init_itsx:
    input: "chimera/{sample}.unique.pick.redundant.fasta"
    output: "itsx/{sample}.5_8S.fasta", "itsx/{sample}.ITS2.fasta", "itsx/{sample}.LSU.fasta", "itsx/{sample}.summary.txt", "itsx/{sample}.positions.txt"
    threads: 3
    log: "logs/{sample}_itsx.log"
    shell:
        "%(itsx)s -t . -i {input} -o itsx/{wildcards.sample} --save_regions 5.8S,ITS2,LSU --complement F --cpu {threads} --graphical F --detailed_results T --partial 50 2> {log}" % config

################ 5.8S processing

rule r58S_extract_58S:
    input: pos="itsx/{sample}.positions.txt", seq="chimera/{sample}.unique.pick.redundant.fasta"
    output: fasta="mothur/{sample}.5_8S_extracted.fasta", gtf="itsx/{sample}.gtf"
    log: "logs/{sample}_5_8Sextraction.log"
    params: minLen=75
    run:
        pos = {}
        for seqId, entry in ITSxParser(open(input.pos)):
            pos[seqId] = entry
        with open(output.fasta, "w") as out, open(output.gtf, "w") as gtf:
            nothingFound = 0
            no58s = 0
            tooShort = 0
            for rec in SeqIO.parse(open(input.seq), "fasta"):
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
                                out.write(rec[0:start].format("fasta"))
        with open(log[0], "a") as logStream:
            logStream.write("-------- 5.8S Extraction --------\n")
            logStream.write("ITSx found nothing in %i sequences\n" % nothingFound)
            logStream.write("ITSx found no 5.8S (ITS2 starts at 0) in %i sequences\n" % no58s)
            logStream.write("ITSx found ITS too close to the start in %i sequences\n" % tooShort)
    

rule r58S_goodReads:
    input: "mothur/{sample}.5_8S_extracted.fasta"
    output: "mothur/{sample}.5_8S_extracted.good.fasta"
    log: "logs/{sample}_mothur.log"
    shell:
        "%(mothur)s -q \"#set.logfile(name={log}); screen.seqs(fasta={input}, maxambig=0, maxn=0);\" > /dev/null " % config
    
rule r58S_uniqueReads:
    input: "mothur/{sample}.5_8S_extracted.good.fasta"
    output: "mothur/{sample}.5_8S_extracted.good.unique.fasta", "mothur/{sample}.5_8S_extracted.good.names"
    log: "logs/{sample}_mothur.log"
    shell:
        "%(mothur)s -q \"#set.logfile(name={log}, append=T); unique.seqs(fasta={input});\" > /dev/null" % config 
    
rule r58S_nonSigletonReads:
    input: fasta="mothur/{sample}.5_8S_extracted.good.unique.fasta", names="mothur/{sample}.5_8S_extracted.good.names"
    output:"mothur/{sample}.5_8S_extracted.good.unique.abund.fasta", "mothur/{sample}.5_8S_extracted.good.abund.names", "mothur/{sample}.5_8S_extracted.good.unique.rare.fasta", "mothur/{sample}.5_8S_extracted.good.rare.names"
    log: "logs/{sample}_mothur.log"
    shell:
        "%(mothur)s -q \"#set.logfile(name={log}, append=T); split.abund(fasta={input.fasta}, name={input.names}, cutoff=1);\" > /dev/null" % config

rule r58S_align:
    input: reads="mothur/{sample}.5_8S_extracted.good.unique.abund.fasta", db="dbs/UNITE.5_8S.aln"
    output: align="mothur/{sample}.5_8S_extracted.good.unique.abund.align", report="mothur/{sample}.5_8S_extracted.good.unique.abund.align.report"
    threads: 3
    log: "logs/{sample}_mothur.log"
    shell:
        "%(mothur)s -q \"#set.logfile(name={log}, append=T); align.seqs(candidate={input.reads}, template={input.db}, processors={threads});\" > /dev/null" % config
        
rule r58S_classify:
    input: aln="mothur/{sample}.5_8S_extracted.good.unique.abund.align", names="mothur/{sample}.5_8S_extracted.good.abund.names", ref="dbs/UNITE.5_8S.aln", tax="dbs/UNITE.5_8S.tax"
    output: "mothur/{sample}.5_8S_extracted.good.unique.abund.5_8S.wang.taxonomy", "mothur/{sample}.5_8S_extracted.good.unique.abund.5_8S.wang.tax.summary"
    log: "logs/{sample}_mothur.log"
    params: cutoff=60
    shell:
        "%(mothur)s -q \"#set.logfile(name={log}, append=T); classify.seqs(fasta={input.aln}, name={input.names}, template={input.ref}, taxonomy={input.tax}, processors={threads}, cutoff={params.cutoff});\" > /dev/null" % config


rule r58S_prepKronaFile:
    input: tax="mothur/{sample}.5_8S_extracted.good.unique.abund.5_8S.wang.taxonomy", name="mothur/{sample}.5_8S_extracted.good.abund.names"
    output: tab="krona/{sample}_5_8S.tsv"
    run:
        count = {}
        for line in open(input.name):
            rId, members = line.strip().split("\t")
            count[rId] = len(members.split(","))
        data = {}
        for line in open(input.tax):
            rId, classification = line.strip().split("\t")
            key = []
            for entry in classification.strip(";").split(";"):
                if entry == "unclassified":
                    break
                elif entry == "unknown":
                    key = [""]
                else:
                    try:
                        name, _ = entry.split("(")
                    except ValueError:
                        print(entry)
                        print(line)
                        raise
                    key.append(name)
            try:
                data["\t".join(key)] += count[rId]
            except KeyError:
                data["\t".join(key)] = count[rId]
        with open(output.tab, "w") as out:
            for key, value in data.items():
                out.write("%i\t%s\n" % (value, key))

def r58S_kronaInput(wildcards):
    return ["krona/%s_5_8S.tsv" %s for s in samples]

rule r58S_krona:
    input: r58S_kronaInput
    output: "krona/5_8s.krona.html"
    shell:
        "%(ktImportText)s -o {output} {input}" % config

################ ITS processing

rule its_copyIts:
    input: "itsx/{sample}.ITS2.fasta"
    output: "mothur/{sample}.ITS2.fasta"
    shell:
        "cp {input} {output}"

rule its_goodReads:
    input: "mothur/{sample}.ITS2.fasta"
    output: "mothur/{sample}.ITS2.good.fasta"
    log: "logs/{sample}_mothur_ITS.log"
    shell:
        "%(mothur)s -q \"#set.logfile(name={log}); screen.seqs(fasta={input}, maxambig=0, maxn=0);\" > /dev/null" % config
        
rule its_uniqueReads:
    input: "mothur/{sample}.ITS2.good.fasta"
    output: "mothur/{sample}.ITS2.good.unique.fasta", "mothur/{sample}.ITS2.good.names"
    log: "logs/{sample}_mothur_ITS.log"
    shell:
        "%(mothur)s -q \"#set.logfile(name={log}, append=T); unique.seqs(fasta={input});\" > /dev/null" % config 
    
rule its_nonSigletonReads:
    input: fasta="mothur/{sample}.ITS2.good.unique.fasta", names="mothur/{sample}.ITS2.good.names"
    output:"mothur/{sample}.ITS2.good.unique.abund.fasta", "mothur/{sample}.ITS2.good.abund.names", "mothur/{sample}.ITS2.good.unique.rare.fasta", "mothur/{sample}.ITS2.good.rare.names"
    log: "logs/{sample}_mothur_ITS.log"
    shell:
        "%(mothur)s -q \"#set.logfile(name={log}, append=T); split.abund(fasta={input.fasta}, name={input.names}, cutoff=1);\" > /dev/null" % config

rule its_prepForClustering:
    input: seq="mothur/{sample}.ITS2.good.unique.abund.fasta", abund="mothur/{sample}.ITS2.good.abund.names"
    output: "swarm/{sample}.ITS2.input.fasta"
    run:
        abund={}
        for line in open(input.abund):
            rep, members = line.strip().split("\t")
            abund[rep] = len(members.split(","))
        with open(output[0], "w") as out:
            for rec in SeqIO.parse(open(input.seq), "fasta"):
                rec.id = "%s_%i" % (rec.id, abund[rec.id])
                out.write(rec.format("fasta"))

rule its_clustering:
    input: "swarm/{sample}.ITS2.input.fasta"
    output: seeds="swarm/{sample}.ITS2.otus.fasta", otuList="swarm/{sample}.ITS2.otus.out"
    log: "logs/{sample}_swarm.log"
    threads: 3
    shell:
        "%(swarm)s -f -t {threads} -w {output.seeds} {input} -o {output.otuList} &> {log}" % config

rule its_get58sClassifications:
    input: its="swarm/{sample}.ITS2.otus.out", tax="mothur/{sample}.5_8S_extracted.good.unique.abund.5_8S.wang.taxonomy", name="mothur/{sample}.5_8S_extracted.good.abund.names"
    output: "taxonomy/{sample}_ITS2.otus_5.8sClass.tsv"
    params: stringency=1.0
    run:
        with open(output[0], "w") as out:
            tsu2tax = {}
            read2tsu = {}
            for line in open(input.name):
                rep, mem = line.strip().split("\t")
                for m in mem.split(","):
                    read2tsu[m] = rep
            for line in open(input.tax):
                seq, lineage = line.strip().split("\t")
                tsu2tax[seq] = lineage  
            for line in open(input.its):
                members = [m.split("|",1)[0] for m in line.strip().split(" ")]
                name = members[0]
                rep58Sclass = []
                for read in members:
                    try:
                         tsu = read2tsu[read]
                    except KeyError:
                        #some reads might have ITS but no 5.8S (eg. because of quality filtering)
                        pass
                    else:
                        rep58Sclass.append(tsu2tax[tsu])
                if rep58Sclass:
                    lcaStr = mothurLCA(rep58Sclass, params.stringency)
                    lcPhylum = ";".join(lcaStr.split(";")[:2]) #trim the lineage doen wot phylum (2. level)
                    out.write("%s\t%s\t%i\n" % (name, lcPhylum, len(rep58Sclass)))
                else:
                    out.write("%s\tunknown\t0\n" % (name))
        

rule its_alignToUnite:
    input: otus="swarm/{sample}.ITS2.otus.fasta", db="dbs/sh_general_release_s_31.01.2016.fasta", dbFlag="dbs/sh_general_release_s_31.01.2016.fasta.lambdaIndexCreated"
    output: "lambda/{sample}.ITS2.otus_vs_UNITE.m8"
    log: "logs/{sample}_lambda.log"
    threads: 3
    shell:
        "%(lambdaFolder)s/lambda -q {input.otus} -d {input.db} -o {output} -p blastn -t {threads} &> {log}" % config

rule its_classify:
    input: lam="lambda/{sample}.ITS2.otus_vs_UNITE.m8", otus="swarm/{sample}.ITS2.otus.fasta"
    output: "taxonomy/{sample}.ITS2.otus.class.tsv"
    params: maxE=1e-6, topPerc=5.0, minIdent=90.0
    run:
        #FIXME: ATTENTION! This is just quick and dirty classification by top 10% rule
        classifi = {}
        for rec in SeqIO.parse(open(input.otus), "fasta"):
            classifi[rec.id.split("|",1)[0]] = []
        for line in open(input.lam, encoding="latin-1"):
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.strip().split("\t")
            readId = qseqid.split("|",1)[0]
            if float(evalue) > params.maxE or float(pident) < params.minIdent:
                continue
            linStr = sseqid.rsplit("|", 1)[-1]
            classifi[readId].append((linStr, float(bitscore)))
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
                    lineage = mothurLCA([hit[0] for hit in sortedHits[:cutoff]])
                    out.write("%s\t%s\n" % (key, lineage))

rule its_kronaPrep:
    input: tax="taxonomy/{sample}.ITS2.otus.class.tsv", otus="swarm/{sample}.ITS2.otus.fasta"
    output: "krona/{sample}.ITS2.otus.krona.tsv"
    run:
        count={}
        for rec in SeqIO.parse(open(input.otus), "fasta"):
            size = int(rec.id.rsplit("_",1)[-1])
            name = rec.id.split("|", 1)[0]
            count[name] = size
        data = {}
        for line in open(input.tax):
            name, lin = line.strip().split("\t")
            try:
                data[lin] += count[name]
            except KeyError:
                data[lin] = count[name]
        with open(output[0], "w") as out:
            for lin, number in data.items():
                if lin == "unknown":
                    out.write("%i\n" % (number))
                else:
                    out.write("%i\t%s\n" % (number, "\t".join(lin.split(";"))))

def its_kronaInput(wildcards):
    return ["krona/%s.ITS2.otus.krona.tsv" % s for s in samples]

rule its_krona:
    input: its_kronaInput
    output: "krona/ITS2.krona.html"
    shell:
        "%(ktImportText)s -o {output} {input}" % config

rule final_combineClassification:
    input: its="taxonomy/{sample}.ITS2.otus.class.tsv", tsu="taxonomy/{sample}_ITS2.otus_5.8sClass.tsv"
    output: "taxonomy/{sample}.ITS2.otus.combClass.tsv", "taxonomy/{sample}.ITS2.otus.conflictingClass.tsv"
    run:
        tsu = {}
        conflict = {}
        for line in open(input.tsu):
            name, lin, number = line.strip().split("\t")
            tsu[name] = lin
        with open(output[0], "w") as out:
            for line in open(input.its):
                name, linStr = line.strip().split("\t")
                lin = linStr.split(";")
                tsuLin = tsu[name].split(";")
                if len(lin) < 2:
                    if lin[0] == tsuLin[0] or (lin[0] == "unknown" and tsuLin[0] != "unknown"):
                        out.write("%s\t%s\n" % (name, tsu[name]))
                    else:
                        out.write("%s\t%s|%s\n" % (name,lin[0],tsuLin[0]))
                        try:
                            conflict[(lin[0], tsuLin[0])] += 1
                        except:
                            conflict[(lin[0], tsuLin[0])] = 1
                else: 
                    if lin[:2] == tsuLin[:2]:
                        out.write("%s\t%s\n" % (name, linStr))
                    else:
                        comLin = []
                        for a,b in zip(lin, tsuLin):
                            if a == b:
                                comLin.append(a)
                            else:
                                comLin.append("%s|%s" % (a,b))
                        out.write("%s\t%s\n" % (name, ";".join(comLin)))
                        try:
                            conflict[(";".join(lin[:2]), ";".join(tsuLin[:2]))] += 1
                        except:
                            conflict[(";".join(lin[:2]), ";".join(tsuLin[:2]))] = 1
        with open(output[1], "w") as out:
            for conf, number in conflict.items():
                out.write("%s | %s\t%i\n" % (conf[0], conf[1], number))

rule final_kronaPrepStep1:
    input: tax="taxonomy/{sample}.ITS2.otus.combClass.tsv", otus="swarm/{sample}.ITS2.otus.fasta"
    output: "krona/{sample}_all.krona.json"
    run:
        count={}
        for rec in SeqIO.parse(open(input.otus), "fasta"):
            size = int(rec.id.rsplit("_",1)[-1])
            name = rec.id.split("|", 1)[0]
            count[name] = size
        data = {"$": []}
        for line in open(input.tax):
            name, linStr = line.strip().split("\t")
            lin = linStr.split(";")
            if lin[0] == "unknown":
                data["$"].append(count[name])
                continue
            here = data
            for entry in lin:
                try:
                    here["$"].append(count[name])
                except KeyError:
                    here["$"] = [count[name]]
                try:
                    here = here[entry]
                except KeyError:
                    here[entry] = {}
                    here = here[entry]
            #add count for the last one
            try:
                here["$"].append(count[name])
            except KeyError:
                here["$"] = [count[name]]
        
        with open(output[0], "w") as out:
            out.write(json.dumps(data))

def final_kronaInput(wildcards):
    return ["krona/%s_all.krona.json" % s for s in samples]

rule final_kronaPrepStep2:
    input:  final_kronaInput
    output: "krona/All.krona.xml"
    run:
        data={}
        for inFile in input:
            sample = "_".join(inFile.split("_")[:2]).split("/")[-1]
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
    shell:
        "%(ktImportXML)s -o {output} {input}" % config

################ helper functions ###################################

class TaxTreeNode(object):
    def __init__(self, name):
        self.name = name
        self.children = {}
        self.reads = {}
        self.otus = {}
        
    def loadSubTree(self, inDict, sample):
        for key, value in inDict.items():
            if key == "$":
                try:
                    self.reads[sample] += sum(value)
                except KeyError:
                    self.reads[sample] = sum(value)
                try:
                    self.otus[sample] += len(value)
                except KeyError:
                    self.otus[sample] = len(value)
            else:
                if key not in self.children:
                    newChild = TaxTreeNode(key)
                    self.children[key] = newChild
                self.children[key].loadSubTree(value, sample)
    
            
    def kronaXml(self, order, indent=""):
        lines=[]
        lines.append('%s<reads>' % indent)
        for sample in order:
            lines.append('%s    <val sample="%s">%i</val>' % (indent, sample, 
                                                              self.reads.get(sample, 0)))
        lines.append('%s</reads>' % indent)
        lines.append('%s<otus>' % indent)
        for sample in order:
            lines.append('%s    <val sample="%s">%i</val>' % (indent, sample, 
                                                              self.otus.get(sample, 0)))
        lines.append('%s</otus>' % indent)
        for name, child in self.children.items():
            lines.append('%s<node name="%s">' % (indent, name.split(";")[-1]))
            lines.extend(child.kronaXml(order, indent+"    "))
            lines.append('%s</node>' % indent)
        return lines


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
