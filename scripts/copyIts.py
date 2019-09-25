from Bio import SeqIO

with open(snakemake.output[0], "w") as out:
    for rec in SeqIO.parse(open(snakemake.input[0]), "fasta"):
        rId, sizeStr = rec.id.split("|")[0].strip(";").split(";")
        rec.id = "%s|ITS2;%s;" % (rId, sizeStr)
        out.write(rec.format("fasta"))
