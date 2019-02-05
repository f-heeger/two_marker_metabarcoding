import glob
import gzip
import itertools

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
    log: "logs/fastqc_{sample}.txt"
    threads: 6
    shell:
        "%(fastqc)s --nogroup -o QC --threads {threads} {input} &> {log}" % config

def qc_multiqc_input(wildcards):
    return ["QC/%s_L001_R%s_001_fastqc.zip" % (s,r) for s,r in itertools.product(samples, ["1","2"])]

rule qc_multiqc:
    input: qc_multiqc_input
    output: "QC/multiqc_report.html", "QC/multiqc_data/multiqc_fastqc.txt"
    log: "logs/multiqc.txt"
    shell:
        "%(multiqc)s -f --interactive -o QC QC/*_fastqc.zip &> {log}" % config

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
        total = 0
        written = 0
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
    input: read1="raw/all_R1.fastq.gz", read2="raw/all_R2.fastq.gz"
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

rule init_trimmStats:
    input: "logs/{sample}_trimmomatic.log"
    output: "logs/{sample}_trimmStats.pdf"
    run:
        R("""
        library(ggplot2)
        d = read.table("{input}", header=F)
        colnames(d) = c("seqName", "read", "len", "firstBase", "lastBase", "trimmed")
        d$readNr = matrix(unlist(strsplit(as.character(d$read), ":")), ncol=4, byrow=T)[,1]
        
        ggplot(d) + geom_histogram(aes(len, fill=readNr), position="dodge", binwidth=5)
        ggsave("{output}")
        """)

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
    params: minOverlap=10
    threads: 3
    log: "logs/all_pear.log"
    shell:
        "%(pear)s -j {threads} -f {input.r1} -r {input.r2} -o merged/all -n %(minAmplLen)s -m %(maxAmplLen)s -v {params.minOverlap} &> {log}" % config

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


