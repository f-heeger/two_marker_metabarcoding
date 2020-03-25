
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
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/getRfamFile.py"


rule db_getUniteFile:
    output: "%(dbFolder)s/sh_general_release_dynamic_%(unite_version)s.zip" % config
    log: "logs/unite_dl.log"
    shell:
        "wget -o {log} -O {output} %(uniteUrl)s;" % config

rule db_unpackUniteFile:
    input: "%(dbFolder)s/sh_general_release_dynamic_%(unite_version)s.zip" % config
    output: "%(dbFolder)s/sh_general_release_dynamic_all_%(unite_version)s.fasta" % config
    log: "logs/unite_unpack.log" % config
    shell:
        "unzip -d %(dbFolder)s %(dbFolder)s/sh_general_release_dynamic_%(unite_version)s.zip &> {log}" % config

rule db_makeUniteFiles:
    input: "%(dbFolder)s/sh_general_release_dynamic_all_%(unite_version)s.fasta" % config
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
        "ITSx -t . -i {input} -o %(dbFolder)s/ITSx/unite_%(unite_version)s --save_regions 5.8S --cpu {threads} --graphical F &> {log}" % config

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
    script:
        "scripts/createTax.py"


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

