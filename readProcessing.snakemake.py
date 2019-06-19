import itertools

rule init_concat:
    output: comb="raw/all_{read}.fastq.gz", sample="readInfo/sample_{read}.tsv", name="readInfo/name_{read}.tsv"
    params: samples=samples
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/concat.py"

rule qc_fastqc:
    input: "%(inFolder)s/{sample}_L001_R{read_number}_001.fastq.gz" % config
    output: "QC/{sample}_L001_R{read_number}_001_fastqc.zip"
    log: "logs/fastqc_{sample}_R{read_number}.txt"
    threads: 6
    conda:
        "envs/fastqc.yaml"
    shell:
        "fastqc --nogroup -o QC --threads {threads} {input} &> {log}" % config

def qc_multiqc_input(wildcards):
    return ["QC/%s_L001_R%s_001_fastqc.zip" % (s,r) for s,r in itertools.product(samples, ["1","2"])]

rule qc_multiqc:
    input: qc_multiqc_input
    output: "QC/multiqc_report.html", "QC/multiqc_data/multiqc_data.json"
    log: "logs/multiqc.txt"
    conda:
        "envs/multiqc.yaml"
    shell:
        "multiqc -f --interactive -o QC QC/*_fastqc.zip &> {log}" % config

rule qc_readCounts:
    input: "QC/multiqc_data/multiqc_data.json"
    output: "readNumbers/rawReadNumbers.tsv"
    script:
        "scripts/qcReadNumber.py"


rule init_filterPrimer:
    input: r1="raw/all_R1.fastq.gz", r2="raw/all_R2.fastq.gz"
    output: r1="primers/all_primers_R1.fastq.gz", r2="primers/all_primers_R2.fastq.gz"
    log: "logs/all_cutadapt.log"
    threads: 6
    conda:
        "envs/cutadapt.yaml"
    shell:
        "cutadapt -g X%(forward_primer)s -G X%(reverse_primer)s --action none --discard-untrimmed --cores {threads} --error-rate %(primerErr)f -o {output.r1} -p {output.r2} {input.r1} {input.r2} > {log}" % config

rule init_primerReadNumbers:
    input: reads="primers/all_primers_R1.fastq.gz", sample="readInfo/sample_R1.tsv"
    output: "readNumbers/primerReadNumbers.tsv"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/primerReadNumber.py"

rule init_trimming:
    input: r1="primers/all_primers_R1.fastq.gz", r2="primers/all_primers_R2.fastq.gz"
    output: r1="trimmed/all_trimmed_R1.fastq.gz", r2="trimmed/all_trimmed_R2.fastq.gz"
    log: "logs/all_timming.log"
    threads: 3
    params: windowLen=8, minQual=20, minLen=200, avgQual=30
    conda:
        "envs/trimmomatic.yaml"
    shell:
        "trimmomatic PE -threads {threads} -phred33 {input.r1} {input.r2} {output.r1} {output.r1}.unpaired {output.r2} {output.r2}.unpaired SLIDINGWINDOW:{params.windowLen}:{params.minQual} TRAILING:{params.minQual} MINLEN:{params.minLen} AVGQUAL:{params.avgQual} &> {log}" % config

rule init_trimmStats:
    input: "logs/{sample}_trimmomatic.log"
    output: "logs/{sample}_trimmStats.pdf"
    conda:
        "envs/ggplot.yaml"
    script:
        "scripts/trimStats.R"

rule init_trimmedReadNumbers:
    input: reads="trimmed/all_trimmed_R1.fastq.gz", sample="readInfo/sample_R1.tsv"
    output: "readNumbers/trimmedReadNumbers.tsv"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/trimmedReadNumber.py"


rule init_merge:
    input: r1="trimmed/all_trimmed_R1.fastq.gz", r2="trimmed/all_trimmed_R2.fastq.gz"
    output: "merged/all.assembled.fastq"
    params: minOverlap=10
    threads: 3
    log: "logs/all_pear.log"
    conda:
        "envs/pear.yaml"
    shell:
        "pear -j {threads} -f {input.r1} -r {input.r2} -o merged/all -n %(minAmplLen)s -m %(maxAmplLen)s -v {params.minOverlap} &> {log}" % config

rule init_convertMerged:
    input: "merged/all.assembled.fastq"
    output: "merged/all.fasta"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/convertMerged.py"

rule init_mergedReadNumbers:
    input: reads="merged/all.assembled.fastq", sample="readInfo/sample_R1.tsv"
    output: "readNumbers/mergedReadNumbers.tsv"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/mergeReadNum.py"

rule init_readNumberOverview:
    input: raw="readNumbers/rawReadNumbers.tsv", primer="readNumbers/primerReadNumbers.tsv", trimmed="readNumbers/trimmedReadNumbers.tsv", merged="readNumbers/mergedReadNumbers.tsv"
    output: "readNumbers/readNumbers.pdf"
    conda:
        "envs/ggplot.yaml"
    script:
        "scripts/plotReadNumber.R"

rule init_dereplicate1:
    input: "merged/all.fasta"
    output: fasta="init_derep/all.derep.fasta", txt="init_derep/all.uc.txt"
    log: "logs/all_repSeq.log"
    conda:
        "envs/vsearch.yaml"
    shell:
        "vsearch --derep_fulllength {input} --output {output.fasta} --uc {output.txt} --sizeout --log {log}"

rule init_dereplicate2:
    input: txt="init_derep/all.uc.txt"
    output: tsv="readInfo/all.repseq.tsv"
    run:
        with open(output.tsv, "w") as out:
            for line in open(input.txt):
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
    conda:
        "envs/vsearch.yaml"
    shell:
        "vsearch --uchime_denovo {input} --nonchimeras {output.fasta} --uchimeout {output.tsv} --log {log}" % config
    

rule init_itsx:
    input: "chimera/all.nochimera.fasta"
    output: "itsx/all.5_8S.fasta", "itsx/all.ITS2.fasta", "itsx/all.LSU.fasta", "itsx/all.summary.txt", "itsx/all.positions.txt"
    threads: 6
    log: "logs/all_itsx.log"
    conda:
        "envs/itsx.yaml"
    shell:
        "ITSx -t . -i {input} -o itsx/all --save_regions 5.8S,ITS2,LSU --complement F --cpu {threads} --graphical F --detailed_results T --partial 50 2> {log}" % config


