import gzip
import time
import csv
from Bio import Entrez
from Bio.Entrez.Parser import ValidationError
from Bio import SeqIO

configfile: "config.json"


rule db_all:
    input: "%(dbFolder)s/58S_derep.fasta" % config, "%(dbFolder)s/58S_tax.tsv" % config

rule db_getRfamFile:
    output: "%(dbFolder)s/RF00002_%(rfam_version)s.fasta.gz" % config
    log: "%(dbFolder)s/logs/rfam_dl.log" % config
    shell:
        "wget -o {log} -O %(dbFolder)s/RF00002_%(rfam_version)s.fasta.gz \"ftp://ftp.ebi.ac.uk/pub/databases/Rfam/%(rfam_version)s/fasta_files/RF00002.fa.gz\"" % config
        
rule db_makeRfamFiles:
    input: fasta="%(dbFolder)s/RF00002_%(rfam_version)s.fasta.gz" % config
    output: fasta="%(dbFolder)s/RF00002_%(rfam_version)s_dna.fasta" % config, tax="%(dbFolder)s/RF00002_%(rfam_version)s_tax.tsv" % config
    log: "%(dbFolder)s/logs/rfam_tax_log.txt" % config
    run:
        ranks = ["superkingdom","kingdom", "phylum", "class", "order", "family", "genus", "species"]
        if not config["email"] or not "@" in config["email"]:
            raise ValueError("'%s' is not a valid email. Set email in config file for NCBI querries." % config["email"])
        nuc2tax = NuclId2TaxIdMap(config["email"], cachePath="%(dbFolder)s/nuc2tax.csv" % config, retry=3)
        tax2lin = LineageMap(config["email"], cachePath="%(dbFolder)s/tax2lin.csv" % config, retry=3)
        noId = 0
        noLin = 0
        with open(output.fasta, "w") as out, open(output.tax, "w") as taxOut, open(log[0], "w") as logFile:
            for rec in SeqIO.parse(gzip.open(input.fasta, "rt"), "fasta"):
                
                
                sId, reg = rec.id.split("/")
                rec.id="%s|%s" % (sId, reg)
                try:
                    tId = nuc2tax[sId]
                except:
                    logFile.write("no ncbi tax ID\t%s\n" % sId)
                    noId += 1
                    continue
                try:
                    lin = tax2lin[tId]
                except:
                    logFile.write("no ncbi lineage\t%s\n" % Id)
                    noLin += 1
                    continue
                linStr = ";".join(["%s__%s" %(l[0][0], l[2]) for l in lin if l[0] in ranks])
                taxOut.write("%s\t%s\n" % (rec.id, linStr))
                out.write(rec.format("fasta"))
        sys.stderr.write("No taxon ID: %i, no lineage: %i\n" % (noId, noLin))
        nuc2tax.save()
        tax2lin.save()


rule db_getUniteFile:
    output: "%(dbFolder)s/sh_general_release_dynamic_%(unite_version)s.fasta" % config
    log: "logs/unite_dl.log"
    shell:
        "cd %(dbFolder)s;"\
        "wget -o {log} -O sh_general_release_dynamic_%(unite_version)s.zip %(uniteUrl)s;"\
        "unzip sh_general_release_dynamic_%(unite_version)s.zip;"\
        "rm sh_general_release_dynamic_%(unite_version)s.zip" % config

rule db_makeUniteFiles:
    input: "%(dbFolder)s/sh_general_release_dynamic_%(unite_version)s.fasta" % config
    output: fasta="%(dbFolder)s/unite_%(unite_version)s.fasta" % config, tax="%(dbFolder)s/unite_%(unite_version)s.tsv" % config, sh2gId="%(dbFolder)s/unite_%(unite_version)s_gIds.tsv" % config
    run:
        with open(output.fasta, "w") as fasta, open(output.tax, "w") as tax, open(output.sh2gId, "w") as sh2gId:
            for line in open(input[0]):
                if line[0] == ">":
                    name, gId, sh, typ, lin = line[1:].strip("\n").split("|")
                    fasta.write(">%s %s\n" % (sh, name))
                    tax.write("%s\t%s\n" % (sh, lin))
                    sh2gId.write("%s\t%s\n" % (sh, gId))
                else:
                    fasta.write(line)

rule db_extract58S:
    input: "%(dbFolder)s/unite_%(unite_version)s.fasta" % config
    output: "%(dbFolder)s/ITSx/unite_%(unite_version)s.5_8S.fasta" % config
    log: "%(dbFolder)s/logs/db_itsx.log" % config
    threads: 6
    conda:
        "envs/itsx.yaml"
    shell:
        "ITSx -t . -i {input} -o %(dbFolder)s/ITSx/unite_%(unite_version)s --save_regions 5.8S --cpu {threads} --graphical F &> {log}"

rule cat58S:
    input: "%(dbFolder)s/ITSx/unite_%(unite_version)s.5_8S.fasta" % config, "%(dbFolder)s/RF00002_%(rfam_version)s_dna.fasta" % config
    output: "%(dbFolder)s/58S.fasta" % config
    shell:
        "cat {input} > {output}"

rule derep:
    input: "%(dbFolder)s/58S.fasta" % config
    output: fasta="%(dbFolder)s/58S_derep.fasta" % config, uc="%(dbFolder)s/58S_derep.uc.txt" % config
    log: "%(dbFolder)s/logs/derep58S.log" % config
    conda:
        "envs/vsearch.yaml"
    shell:
        "vsearch --derep_fulllength {input} --output {output.fasta} --uc {output.uc} --sizeout --log {log} &> /dev/null"

rule createTax:
    input: uc="%(dbFolder)s/58S_derep.uc.txt" % config, uTax="%(dbFolder)s/unite_%(unite_version)s.tsv" % config, rTax="%(dbFolder)s/RF00002_%(rfam_version)s_tax.tsv" % config
    output: "%(dbFolder)s/58S_tax.tsv" % config
    run:
        clu={}
        for line in open(input.uc):
            lType, cNum, sLen, ident, strand, _1, _2, aln, query, target = line.strip("\n").split("\t")
            if lType == "S":
                clu[query] = [query]
            elif lType == "H":
                clu[target].append(query)
            elif lType == "C":
                assert int(sLen) == len(clu[query])
            else:
                raise ValueError("Unknown record type: %s" % lType)
        uTax = {}
        for line in open(input.uTax):
            tId, tLin = line.strip().split("\t")
            uTax[tId] = tLin
        rTax = {}
        for line in open(input.rTax):
            tId, tLin = line.strip().split("\t")
            rTax[tId] = tLin
        with open(output[0], "w") as out:
            for rep, memList in clu.items():
                linStrs = []
                for m in memList:
                    if m[:2] == "SH":
                        linStrs.append(uTax[m.split("|")[0]])
                    else:
                        lin = rTax[m]
                        arr = lin.split(";")
                        if arr[1] == "k__Fungi":
                            #ignore everything RFAM has to say about fungi taxonomy
                            linStrs.append("k__Fungi")
                        else:
                            linStrs.append(";".join(arr[1:]))
                lcaLin = lca(linStrs, 0.95)
                out.write("%s\t%s\n" % (rep, lcaLin))

rule db_creat58SIndex:
    input: "%(dbFolder)s/58S_derep.fasta" % config
    output: touch("%(dbFolder)s/58S_derep.fasta.lambdaIndexCreated" % config)
    threads: 6
    conda:
        "envs/lambda.yaml"
    shell:
        "lambda_indexer -d {input} -p blastn -t {threads}"

rule db_creatUniteIndex:
    input: "%(dbFolder)s/unite_%(unite_version)s.fasta" % config
    output: touch("%(dbFolder)s/unite_%(unite_version)s.fasta.lambdaIndexCreated" % config)
    threads: 6
    conda:
        "envs/lambda.yaml"
    shell:
        "lambda_indexer -d {input} -p blastn -t {threads}"

