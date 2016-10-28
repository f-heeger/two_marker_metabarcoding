import itertools
import json
import xml.etree.ElementTree as et
import gzip
import glob
import random

from snakemake.utils import min_version, R
from Bio import SeqIO

min_version("3.5.4")

configfile: "config.json"

shell.prefix("sleep 10; ") #work aorund to desl with "too quck" rule execution and slow NAS

#samples =  {"A1_S1": "A1_S1", "B1_S2": "B1_S2", "H5_S40": "H5_S40", "A5_S33": "A5_S33"}

samples = config["samples"].keys()

rule all:
    input: "krona/All.krona.html", "krona/5_8s.krona.html", "krona/ITS2.krona.html", "taxonomy/all.compareClass.tsv", "otu_table.tsv", "All.rarefactions.pdf"
#    input: "taxonomy/All.krona.html", "taxonomy/5_8s.krona.html", expand(["taxonomy/{sample}_ITS2.otus_5.8sClass.tsv", "taxonomy/{sample}.ITS2.otus.class.tsv"], sample=samples)

################ generate reference sequence for 5.8S ##########################

#include: "generate58SDatabase.snakefile.py"

################ generate lamda db from UNTIE ##################################

rule db_getUniteFile:
    output: "%(dbFolder)s/sh_general_release_dynamic_22.08.2016.fasta" % config
    shell:
        "cd %(dbFolder)s;" \
        "wget https://unite.ut.ee/sh_files/sh_general_release_22.08.2016.zip;" \
        "unzip sh_general_release_22.08.2016.zip;" \
        "rm sh_general_release_22.08.2016.zip" % config

rule db_creatUniteIndex:
    input: "%(dbFolder)s/sh_general_release_dynamic_22.08.2016.fasta" % config
    output: touch("%(dbFolder)s/sh_general_release_dynamic_22.08.2016.fasta.lambdaIndexCreated" % config)
    threads: 6
    shell:
        "%(lambdaFolder)s/lambda_indexer -d {input} -p blastn -t {threads}" % config

################ quality control ###############################################

rule init_concat:
    output: comb="raw/all_{read}.fastq.gz", sample="readInfo/sample_{read}.tsv", name="readInfo/name_{read}.tsv"
    run:
        with gzip.open(output.comb, "wt") as combOut, \
             open(output.sample, "w") as sampleOut, \
             open(output.name, "w") as nameOut:
            for sample in samples:
                path = "%s/%s_L*_%s_*.fastq.gz" % (config["inFolder"], sample, wildcards.read)
                inFiles=glob.glob(path)
                if not inFiles:
                    raise RuntimeError("No file(s) found for sample %s at %s." % (sample, path))
                for inFile in inFiles:
                    with gzip.open(inFile, "rt") as inStream:
                        for rec in SeqIO.parse(inStream, "fastq"):
                            newId = "_".join(rec.id.split(":")[3:])
                            nameOut.write("%s\t%s\n" % (rec.id, newId))
                            rec.id = newId
                            combOut.write(rec.format("fastq"))
                            sampleOut.write("%s\t%s\n" % (newId, sample))

rule qc_fastqc:
    input: "%(inFolder)s/{sample}_L001_R{read_number}_001.fastq.gz" % config
    output: "QC/{sample}_L001_R{read_number}_001_fastqc.zip"
    threads: 6
    shell:
        "%(fastqc)s --nogroup -o QC --threads {threads} {input}" % config

def qc_multiqc_input(wildcards):
    return ["QC/%s_L001_R%s_001_fastqc.zip" % (s,r) for s,r in itertools.product(samples, ["1","2"])]

rule qc_multiqc:
    input: qc_multiqc_input
    output: "QC/multiqc_report.html", "QC/multiqc_data/multiqc_fastqc.txt"
    shell:
        "%(multiqc)s -f --interactive -o QC QC/*_fastqc.zip" % config

rule qc_readCounts:
    input: "QC/multiqc_data/multiqc_fastqc.txt"
    output: "readNumbers/rawReadNumbers.tsv"
    run:
        with open(output[0], "w") as out:
            for l, line in enumerate(open(input[0])):
                if l==0:
                    continue #skip header
                arr = line.strip("\n").split("\t")
                #Sample	avg_sequence_length	percent_dedup	percent_duplicates	seq_len_read_count	seq_len_bp	sequence_length	total_sequences	percent_gc
                sample, lane, readNum, part = arr[0].rsplit("_", 3)
                if readNum != "R1":
                    continue # only use R1
                seqNum = int(float(arr[-2]))
                out.write("%s\t%i\n" % (sample, seqNum))

rule init_filterByBarcodeQuality:
    input: read1="raw/all_R1.fastq.gz", read2="raw/all_R2.fastq.gz", index1="raw/all_I1.fastq.gz", index2="raw/all_I2.fastq.gz"
    output: "raw/all_goodIndex_R1.fastq.gz", "raw/all_goodIndex_R2.fastq.gz"
    log: "logs/all_indexQualFilter.log"
    run:
        total=0
        written=0
        with gzip.open(output[0], "wt") as out1, gzip.open(output[1], "wt") as out2:
            with gzip.open(input.read1, "rt") as read1File, gzip.open(input.read2, "rt") as read2File, gzip.open(input.index1, "rt") as index1File, gzip.open(input.index2, "rt") as index2File:
                read1 = SeqIO.parse(read1File, "fastq")
                read2 = SeqIO.parse(read2File, "fastq")
                index1 = SeqIO.parse(index1File, "fastq")
                index2 = SeqIO.parse(index2File, "fastq")
                while True:
                    try:
                        r1 = next(read1)
                    except StopIteration:
                        break
                    r2 = next(read2)
                    i1 = next(index1)
                    i2 = next(index2)
                    total += 1
                    minQi1 = min(i1._per_letter_annotations["phred_quality"])
                    minQi2 = min(i2._per_letter_annotations["phred_quality"])
                    if min(minQi1, minQi2) >= 20:
                        out1.write(r1.format("fastq"))
                        out2.write(r2.format("fastq"))
                        written += 1
        open(log[0], "w").write("%i of %i remain after index filtering (%.2f%%)" % (written, total, float(written)/total*100))


rule init_indexQualityReadNumbers:
    input: reads="raw/all_goodIndex_R1.fastq.gz", sample="readInfo/sample_R1.tsv"
    output: "readNumbers/indexQualityReadNumbers.tsv"
    run:
        readSample = {}
        for line in open(input.sample):
            read, sample = line.strip("\n").split("\t")
            readSample[read] = sample
        sampleReads = {}
        for rec in SeqIO.parse(gzip.open(input.reads, "rt"), "fastq"):
            sample = readSample[rec.id]
            try:
                sampleReads[sample] += 1
            except KeyError:
                sampleReads[sample] = 1
        with open(output[0], "w") as out:
            for sample, value in sampleReads.items():
                out.write("%s\t%s\n" % (sample, value))

rule init_filterPrimer:
    input: read1="raw/all_goodIndex_R1.fastq.gz", read2="raw/all_goodIndex_R2.fastq.gz"
    output: "primers/all_barcode_ITS-ITS_1.fastq.gz", "primers/all_barcode_ITS-ITS_2.fastq.gz"
    log: "logs/all_flexbar.log"
    threads: 6
    shell:
        "echo \">ITS\n%(forward_primer)s\" > primers/fwd_primer.fasta; echo \">ITS\n%(reverse_primer)s\" > primers/rev_primer.fasta; %(flexbar)s -n {threads} -r {input.read1} -p {input.read2} -t primers/all -b primers/fwd_primer.fasta -b2 primers/rev_primer.fasta -bk -be LEFT_TAIL -bn 25 -bt %(primerErr)f -f i1.8 -u 100 -z GZ &> {log}" % config

rule init_primerReadNumbers:
    input: reads="primers/all_barcode_ITS-ITS_1.fastq.gz", sample="readInfo/sample_R1.tsv"
    output: "readNumbers/primerReadNumbers.tsv"
    run:
        readSample = {}
        for line in open(input.sample):
            read, sample = line.strip("\n").split("\t")
            readSample[read] = sample
        sampleReads = {}
        for rec in SeqIO.parse(gzip.open(input.reads, "rt"), "fastq"):
            sample = readSample[rec.id]
            try:
                sampleReads[sample] += 1
            except KeyError:
                sampleReads[sample] = 1
        with open(output[0], "w") as out:
            for sample, value in sampleReads.items():
                out.write("%s\t%s\n" % (sample, value))

rule init_trimming:
    input: r1="primers/all_barcode_ITS-ITS_1.fastq.gz", r2="primers/all_barcode_ITS-ITS_2.fastq.gz"
    output: r1="trimmed/all_trimmed_R1.fastq.gz", r2="trimmed/all_trimmed_R2.fastq.gz"
    log: "logs/all_timming.log"
    threads: 3
    params: windowLen=8, minQual=20, minLen=200, avgQual=30
    shell:
        "java -jar %(trimmomatic)s PE -threads {threads} -phred33 {input.r1} {input.r2} {output.r1} {output.r1}.unpaired {output.r2} {output.r2}.unpaired SLIDINGWINDOW:{params.windowLen}:{params.minQual} TRAILING:{params.minQual} MINLEN:{params.minLen} AVGQUAL:{params.avgQual} &> {log}" % config

#rule init_trimmStats:
#    input: "logs/{sample}_trimmomatic.log"
#    output: "logs/{sample}_trimmStats.pdf"
#    run:
#        R("""
#        library(ggplot2)
#        d = read.table("{input}", header=F)
#        colnames(d) = c("seqName", "read", "len", "firstBase", "lastBase", "trimmed")
#        d$readNr = matrix(unlist(strsplit(as.character(d$read), ":")), ncol=4, byrow=T)[,1]
#        
#        ggplot(d) + geom_histogram(aes(len, fill=readNr), position="dodge", binwidth=5)
#        ggsave("{output}")
#        """)

rule init_trimmedReadNumbers:
    input: reads="trimmed/all_trimmed_R1.fastq.gz", sample="readInfo/sample_R1.tsv"
    output: "readNumbers/trimmedReadNumbers.tsv"
    run:
        readSample = {}
        for line in open(input.sample):
            read, sample = line.strip("\n").split("\t")
            readSample[read] = sample
        sampleReads = {}
        for rec in SeqIO.parse(gzip.open(input.reads, "rt"), "fastq"):
            sample = readSample[rec.id]
            try:
                sampleReads[sample] += 1
            except KeyError:
                sampleReads[sample] = 1
        with open(output[0], "w") as out:
            for sample, value in sampleReads.items():
                out.write("%s\t%s\n" % (sample, value))

rule init_merge:
    input: r1="trimmed/all_trimmed_R1.fastq.gz", r2="trimmed/all_trimmed_R2.fastq.gz"
    output: "merged/all.assembled.fastq"
    params: minLen=300, maxLen=550, minOverlap=10
    threads: 3
    log: "logs/all_pear.log"
    shell:
        "%(pear)s -j {threads} -f {input.r1} -r {input.r2} -o merged/all -n {params.minLen} -m {params.maxLen} -v {params.minOverlap} &> {log}" % config

rule init_convertMerged:
    input: "merged/all.assembled.fastq"
    output: "merged/all.fasta"
    run:
        with open(output[0], "w") as out:
            for record in SeqIO.parse(open(input[0]), "fastq"):
                out.write(record.format("fasta"))

rule init_mergedReadNumbers:
    input: reads="merged/all.assembled.fastq", sample="readInfo/sample_R1.tsv"
    output: "readNumbers/mergedReadNumbers.tsv"
    run:
        readSample = {}
        for line in open(input.sample):
            read, sample = line.strip("\n").split("\t")
            readSample[read] = sample
        sampleReads = {}
        for rec in SeqIO.parse(open(input.reads), "fastq"):
            sample = readSample[rec.id]
            try:
                sampleReads[sample] += 1
            except KeyError:
                sampleReads[sample] = 1
        with open(output[0], "w") as out:
            for sample, value in sampleReads.items():
                out.write("%s\t%s\n" % (sample, value))

rule init_readNumberOverview:
    input: raw="readNumbers/rawReadNumbers.tsv", indexQual="readNumbers/indexQualityReadNumbers.tsv", primer="readNumbers/primerReadNumbers.tsv", trimmed="readNumbers/trimmedReadNumbers.tsv", merged="readNumbers/mergedReadNumbers.tsv"
    output: "readNumbers/readNumbers.pdf"
    run:
        R("""
        library(ggplot2)
        raw = read.table("{input.raw}", header=F)
        raw = raw[order(raw$V1),]
        raw$stage = "raw"
        d = raw
        iQual = read.table("{input.indexQual}", header=F)
        iQual = iQual[order(iQual$V1),]
        iQual$stage = "indexQualityFiltered"
        d = rbind(d, iQual)
        primer = read.table("{input.primer}", header=F)
        primer = primer[order(primer$V1),]
        primer$stage = "primerFound"
        d = rbind(d, primer)
        trimmed = read.table("{input.trimmed}", header=F)
        trimmed = trimmed[order(trimmed$V1),]
        trimmed$stage = "trimmed"
        d=rbind(d, trimmed)
        merged = read.table("{input.merged}", header=F)
        merged = merged[order(merged$V1),]
        merged$stage = "merged"
        d=rbind(d,merged)
        colnames(d) = c("sample", "readNum", "stage")
        d$stage = factor(d$stage, levels=c("merged", "trimmed", "primerFound", "indexQualityFiltered", "raw"))
        ggplot(d) + geom_bar(aes(sample, readNum, fill=stage), stat="identity", position="dodge") + coord_flip() + geom_hline(yintercept = 10000, linetype="dashed")
        ggsave("{output}", width=7, height=28, units="in")
        """)

rule init_dereplicate:
    input: "merged/all.fasta"
    output: fasta="init_derep/all.derep.fasta", tsv="readInfo/all.repseq.tsv", txt="init_derep/all.uc.txt"
    log: "logs/all_repSeq.log"
    run:
        shell("%(vsearch)s --derep_fulllength {input} --output {output.fasta} --uc {output.txt} --sizeout --log {log}" % config)
        with open(output.tsv, "w") as out:
            for line in open(output.txt):
                arr = line.strip().split("\t")
                if arr[0] == "C":
                    pass
                elif arr[0] == "S":
                    seq = arr[-2]
                    out.write("%s\t%s\n" % (seq, seq))
                elif arr[0] == "H":
                    seq, cluster = arr[-2:]
                    out.write("%s\t%s\n" % (seq, cluster))
                else:
                    raise ValueError("Unknown record type: %s" % arr[0])
                

rule init_removeChimera:
    input: "init_derep/all.derep.fasta"
    output: fasta="chimera/all.nochimera.fasta", tsv="chimera/all.chimeraReport.tsv"
    log: "logs/all_chimera.log"
    shell:
        "%(vsearch)s --uchime_denovo {input} --nonchimeras {output.fasta} --uchimeout {output.tsv} --log {log}" % config
    

rule init_itsx:
    input: "chimera/all.nochimera.fasta"
    output: "itsx/all.5_8S.fasta", "itsx/all.ITS2.fasta", "itsx/all.LSU.fasta", "itsx/all.summary.txt", "itsx/all.positions.txt"
    threads: 6
    log: "logs/all_itsx.log"
    shell:
        "%(itsx)s -t . -i {input} -o itsx/all --save_regions 5.8S,ITS2,LSU --complement F --cpu {threads} --graphical F --detailed_results T --partial 50 2> {log}" % config

################ 5.8S processing

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
                                oldId, size = rec.id.strip(";").split(";")
                                rec.id = "%s|5.8S;%s;" % (oldId, size)
                                rec.description="Extracted 5.8S sequence 0-%i" % (start-1)
                                out.write(rec[0:start].format("fasta"))
        with open(log[0], "a") as logStream:
            logStream.write("-------- 5.8S Extraction --------\n")
            logStream.write("ITSx found nothing in %i sequences\n" % nothingFound)
            logStream.write("ITSx found no 5.8S (ITS2 starts at 0) in %i sequences\n" % no58s)
            logStream.write("ITSx found ITS too close to the start in %i sequences\n" % tooShort)


rule r58S_removePrimerAndNs:
    input: "itsx/all.5_8S_extracted.fasta"
    output: "primerremoved/all.5_8S_primerRemoved.fasta"
    log: "logs/all_58S_cutadapt.log"
    params: minOverlap=10, maxN=1
    shell:
        "%(cutadapt)s -g ^%(forward_primer)s --trimmed-only -O {params.minOverlap} --max-n={params.maxN} -o {output} {input} &> {log}" % config

rule r58S_dereplicate:
    """dereplicate 5.8S sequences and only retain "clusters" with more than one sequence"""
    input: "primerremoved/all.5_8S_primerRemoved.fasta"
    output: fasta="mothur/all.5_8S_derep.fasta", tsv="readInfo/all.rep58S.tsv", txt="r58S_derep/all.uc.txt"
    log: "logs/all_rep58S.log"
    params: minsize=2
    run:
        shell("%(vsearch)s --derep_fulllength {input} --output {output.fasta} --uc {output.txt} --sizein --sizeout --minuniquesize {params.minsize} --log {log}" % config)
        
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

rule r58S_align:
    input: reads="mothur/all.5_8S_derep.fasta", db="%(dbFolder)s/UNITE.5_8S.aln" % config
    output: align="mothur/all.5_8S_derep.align", report="mothur/all.5_8S_derep.align.report"
    threads: 3
    log: "logs/58S_mothur.log"
    shell:
        "%(mothur)s -q \"#set.logfile(name={log}, append=T); align.seqs(candidate={input.reads}, template={input.db}, processors={threads});\" > /dev/null" % config
        
rule r58S_classify:
    input: aln="mothur/all.5_8S_derep.align", ref="%(dbFolder)s/UNITE.5_8S.aln" % config, tax="%(dbFolder)s/UNITE.5_8S.tax" % config
    output: "mothur/all.5_8S_derep.5_8S.wang.taxonomy", "mothur/all.5_8S_derep.5_8S.wang.tax.summary"
    log: "logs/58s_mothur.log"
    params: cutoff=60
    shell:
        "%(mothur)s -q \"#set.logfile(name={log}, append=T); classify.seqs(fasta={input.aln}, template={input.ref}, taxonomy={input.tax}, processors={threads}, cutoff={params.cutoff});\" > /dev/null" % config


rule r58S_readClassification:
    input: tax="mothur/all.5_8S_derep.5_8S.wang.taxonomy", repseq="readInfo/all.repseq.tsv", rep58s="readInfo/all.rep58S.tsv", sample="readInfo/sample_R1.tsv"
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
            rName, classification = line.strip().split("\t")
            rId, sizeStr = rName.strip(";").split(";")
            count = int(sizeStr.split("=")[1])
            cls = []
            for entry in classification.strip(";").split(";"):
                if entry == "unclassified":
                    break
                elif entry == "unknown":
                    cls = [""]
                else:
                    try:
                        name, _ = entry.split("(")
                    except ValueError:
                        print(entry)
                        print(line)
                        raise
                    cls.append(name)
            i=0
            for r58seq in rep58s[rId]:
                for read in repseq[r58seq.split("|")[0]]:
                    readClass[read] = ";".join(cls)
                    i+=1
            assert i==count
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
    shell:
        "%(ktImportText)s -o {output} {input}" % config

################ ITS processing

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
    run:
        shell("%(vsearch)s --derep_fulllength {input} --output {output.fasta} --uc {output.txt} --sizein --sizeout --minuniquesize {params.minsize} --log {log}" % config)
        
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
    shell:
        "%(swarm)s -f -z -t {threads} -w {output.seeds} {input} -o {output.otuList} &> {log}" % config


rule its_alignToUnite:
    input: otus="swarm/all.ITS2.otus.fasta", db="%(dbFolder)s/UNITE_public_31.01.2016.ascii.good.fasta" % config, dbFlag="%(dbFolder)s/UNITE_public_31.01.2016.ascii.fasta.lambdaIndexCreated" % config
    output: "lambda/all.ITS2.otus_vs_UNITE.m8"
    log: "logs/all_lambda.log"
    threads: 3
    shell:
        "%(lambdaFolder)s/lambda -q {input.otus} -d {input.db} -o {output} -p blastn -t {threads} &> {log}" % config

rule its_classify:
    input: lam="lambda/all.ITS2.otus_vs_UNITE.m8", otus="swarm/all.ITS2.otus.fasta"
    output: "taxonomy/all.ITS2.otus.class.tsv"
    params: maxE=1e-6, topPerc=5.0, minIdent=80.0, minCov=85.0, stringency=.90
    log: "logs/its_class.log"
    run:
        logOut = open(log[0], "w")
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
            itsLength[rec.id.split("|",1)[0]] = len(rec)
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
            linStr = sseqid.rsplit("|", 1)[-1]
            if linStr.endswith("Incertae"):
                linStr += "_sedis" #FIXME: workaroud for taking the tayonomy from the fasta header which ends at a space. Wither fix fasta headers or take taxonomy from UNITE taxonomy file
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
                    lcaStr =  lca(otuCls, params.stringency, unidentified=["unclassified", "unknown"])
                    lcPhylum = ";".join(lcaStr.split(";")[:2]) #trim the lineage down to phylum (2. level)
                    out.write("%s\t%s\t%i\n" % (otu, lcPhylum, len(otuCls)))
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
    shell:
        "%(ktImportText)s -o {output} {input}" % config

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
                            out.write("%s\t%s|%s\n" % (name, lin[0], tsuLin[0]))
                        elif ifConf==1:
                            out.write("%s\t%s\n" % (name, ";".join(tsuLin)))
                        else:
                            out.write("%s\t%s\n" % (name, ";".join(lin)))
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
                            out.write("%s\t%s|%s\n" % (name, lin[0], tsuLin[0]))
                        elif ifConf==1:
                            out.write("%s\t%s\n" % (name, ";".join(tsuLin)))
                        else:
                            out.write("%s\t%s\n" % (name, ";".join(lin)))
                        try:
                            conflict[(lin[0], tsuLin[0])] += 1
                        except:
                            conflict[(lin[0], tsuLin[0])] = 1
                else:
                    #if the ITS classification is more than one level deep
                    # check if the first two levels of classifications are identical
                    if lin[:2] == tsuLin[:2]:
                        #if they are the same accept the ITS classification
                        # (including all levels after the second)
                        out.write("%s\t%s\n" % (name, linStr))
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
    run:
        R("""library(ggplot2)
        
            raw=read.table("{input}")
            colnames(raw) = c("sample", "x", "y")

            allPlotData = data.frame()
            anno=data.frame()

            for (tSample in unique(raw$sample)){{

                print(tSample)
                d = subset(raw, sample==tSample)
                plotData = data.frame(unique(d$x))
                colnames(plotData) = c("x")
                N=length(plotData$x)

                plotData$mean=NA
                plotData$cmin=NA
                plotData$cmax=NA
                plotData$sample=tSample

                conf = 0.95

                for (i in 1:length(plotData$x)) {{
                    x_i = plotData$x[i]
                    s = subset(d, x==x_i)
                    plotData$mean[i] = mean(s$y)
                    plotData$cmin[i] = sort(s$y)[floor(length(s$y)*(1-conf)/2)]
                    plotData$cmax[i] = sort(s$y)[ceiling(length(s$y)*(1-(1-conf)/2))]
                }}
                anno=rbind(anno, plotData[plotData$x==max(plotData$x),])

            allPlotData = rbind(allPlotData, plotData)

            }}

            ggplot(allPlotData, aes(x, mean)) + geom_point(aes(color=sample), size=0.1) + geom_ribbon(aes(ymin=cmin, ymax=cmax, fill=sample), alpha=0.2) + annotate("text", x=anno$x, y=anno$mean, label=anno$sample)
            ggsave("{output}", width=16, height=9)""")

rule final_creatOtuTable:
    input: tax="taxonomy/all.ITS2.otus.combClass.tsv", otus=expand("swarm/{sample}.ITS2.otus.out", sample=samples)
    output: "otu_table.tsv"
    run:
        readNr = {}
        for sample in samples:
            readNr[sample] = {}
            for line in open("swarm/%s.ITS2.otus.out" % sample):
                otuStr, nr = line.strip("\n").split("\t")
                otu = otuStr.split("|")[0]
                readNr[sample][otu] = nr
        with open(output[0], "w") as out:
            out.write("otu_ID\tclassification\t%s\n" % "\t".join(samples))
            for line in open(input.tax):
                otu, cls = line.strip("\n").split("\t")
                numbers = [readNr[sample].get(otu, "0") for sample in samples]
                out.write("%s\t%s\t%s\n" % (otu, cls, "\t".join(numbers)))

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
            lin = linStr.split(";")
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
