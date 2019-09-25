
rule r58S_extract_58S:
    input: pos="itsx/all.positions.txt", seq="chimera/all.nochimera.fasta"
    output: fasta="itsx/all.5_8S_extracted.fasta", gtf="itsx/all.gtf"
    log: "logs/all_5_8Sextraction.log"
    params: minLen=75
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/extract_58s.py"


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
    params: minsize=2
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
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/classify58s.py"


rule r58S_readClassification:
    input: tax="taxonomy/all.58S.derep.class.tsv", repseq="readInfo/all.repseq.tsv", rep58s="readInfo/all.rep58S.tsv", sample="readInfo/sample_R1.tsv"
    output: "taxonomy/all_5_8S_classification.tsv"
    script:
        "scripts/classify58sReads.py"


rule r58S_prepKronaInput:
    input: cls="taxonomy/all_5_8S_classification.tsv", sample="readInfo/sample_R1.tsv"
    output: expand("krona/{sample}_5_8S.tsv", sample=samples)
    script:
        "scripts/prep58sKronaInput.py"

rule r58S_krona:
    input: expand("krona/{sample}_5_8S.tsv", sample=samples)
    output: "krona/5_8s.krona.html"
    log:
        "logs/58S_krona.log"
    conda:
        "envs/krona.yaml"
    shell:
        "ktImportText -o {output} {input} &> {log}"

