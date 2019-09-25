import random
import json

rule final_concatMarkers:
    input: otus="swarm/all.ITS2.otus.fasta", tsu="primerremoved/all.5_8S_primerRemoved.fasta"
    output: "all.concatMarkers.fasta"
    log: "logs/all.concatMarker.log"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/concatMarkers.py"

rule final_combineClassification:
    input: itsCls="taxonomy/all.ITS2.otus.class.tsv", r58sCls="taxonomy/all_ITS2.otus_5.8sClass.tsv"
    output: otuComb="taxonomy/all.ITS2.otus.combClass.tsv", conflict="taxonomy/all.ITS2.otus.conflictingClass.tsv"
    script:
        "scripts/combineClass.py"

rule final_classCompare:
    input: itsCls="taxonomy/all.ITS2.otus.class.tsv", otuReads="swarm/all.otuReads.tsv", r58SreadCls="taxonomy/all_5_8S_classification.tsv", r58sOtuCls="taxonomy/all_ITS2.otus_5.8sClass.tsv", combCls="taxonomy/all.ITS2.otus.combClass.tsv", sample="readInfo/sample_R1.tsv"
    output: comp="taxonomy/all.compareClass.tsv", stat="taxonomy/all.clsStat.tsv"
    script:
        "scripts/classCompare.py"

rule final_diffClsDepth:
    input: "taxonomy/all.clsStat.tsv"
    output: "taxonomy/all.otuClsDepth.tsv"
    script:
        "scripts/diffClsDepth.py"

rule final_phylumDiff:
    input: "otu_table.tsv"
    output: "fungiPhylumDiff.tsv"
    script:
        "scripts/phylumDiff.py"

rule final_plotPhylumDiff:
    input: "fungiPhylumDiff.tsv"
    output: "gainFrom58s.pdf"
    conda:
        "envs/ggplot.yaml"
    script:
        "scripts/plotPhylumDiff.R"

rule final_computeRarefaction:
    input: "swarm/{sample}.ITS2.otus.out"
    output: "rarefaction/{sample}.rarefaction.tsv"
    params: rep=100, step=1000
    script:
        "scripts/computeRarefaction.py"

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
    conda:
        "envs/ggplot.yaml"
    script:
        "scripts/plotRarefaction.R"

rule final_creatOtuTable:
    input: tax58s="taxonomy/all_ITS2.otus_5.8sClass.tsv", taxIts="taxonomy/all.ITS2.otus.class.tsv", taxComb="taxonomy/all.ITS2.otus.combClass.tsv", otus=expand("swarm/{sample}.ITS2.otus.out", sample=samples)
    output: "otu_table.tsv"
    params: samples=samples
    script:
        "scripts/createOtuTab.py"

rule final_kronaPrepStep1:
    input: tax="taxonomy/all.ITS2.otus.combClass.tsv", otuReads = "swarm/{sample}.ITS2.otus.out"
    output: "krona/{sample}_all.krona.json"
    script:
        "scripts/finalKronaPrep1.py"


rule final_kronaPrepStep2:
    input:  expand("krona/{sample}_all.krona.json", sample=samples)
    output: "krona/All.krona.xml"
    script:
        "scripts/finalKronaPrep2.py"

rule final_kronaAll:
    input: "krona/All.krona.xml"
    output: "krona/All.krona.html"
    log: "logs/kronaAll.log"
    conda:
        "envs/krona.yaml"
    shell:
        "ktImportXML -o {output} {input} &> {log}"
