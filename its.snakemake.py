
rule its_copyIts:
    input: "itsx/all.ITS2.fasta"
    output: "itsx/all.ITS2_extracted.fasta"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/copyIts.py"

rule its_goodReads:
    input: "itsx/all.ITS2_extracted.fasta"
    output: "itsx/all.ITS2_extracted.good.fasta"
    log: "logs/all_good_ITS.log"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/goodItsReads.py"

rule its_dereplicate1:
    input: "itsx/all.ITS2_extracted.good.fasta"
    output: fasta="its_derep/all.ITS2_extracted.derep.fasta", txt="its_derep/all.uc.txt"
    log: "logs/all_repIts.log"
    params: minsize=2
    conda:
        "envs/vsearch.yaml"
    shell:
        "vsearch --derep_fulllength {input} --output {output.fasta} --uc {output.txt} --sizein --sizeout --minuniquesize {params.minsize} --log {log}"

rule its_dereplicate2:
    input: txt="its_derep/all.uc.txt"
    output: tsv="readInfo/all.repITS2.tsv"
    params: minsize=2
    script:
        "scripts/itsDereplicate.py"


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
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/itsClassify.py"

rule its_createOtuReads:
    input: otuList="swarm/all.ITS2.otus.out", repIts="readInfo/all.repITS2.tsv", repSeq="readInfo/all.repseq.tsv"
    output: "swarm/all.otuReads.tsv"
    script:
        "scripts/otuReads.py"


rule its_readClassification:
    input: tax="taxonomy/all.ITS2.otus.class.tsv", otuList="swarm/all.ITS2.otus.out", otuReads="swarm/all.otuReads.tsv"
    output: cls="taxonomy/all.ITS2.classification.tsv"
    script:
        "scripts/itsReadsClassification.py"


rule its_get58sClassifications:
    input: otuReads="swarm/all.otuReads.tsv", r58SreadCls="taxonomy/all_5_8S_classification.tsv"
    output: "taxonomy/all_ITS2.otus_5.8sClass.tsv"
    log: "logs/all_get58sClass.log"
    params: stringency=0.90
    script:
        "scripts/get58sClass.py"

rule its_kronaPrep:
    input: cls="taxonomy/all.ITS2.classification.tsv", sample="readInfo/sample_R1.tsv"
    output: tab=expand("krona/{sample}_ITS2.krona.tsv", sample=samples)
    script:
        "scripts/itsKronaPrep.py"


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
    script:
        "scripts/perSampleOtuReads.py"


