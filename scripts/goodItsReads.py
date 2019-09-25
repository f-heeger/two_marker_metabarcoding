from Bio import SeqIO

with open(snakemake.output[0], "w") as out:
    for rec in SeqIO.parse(open(snakemake.input[0]), "fasta"):
        if rec.seq.count("N") == 0:
            out.write(rec.format("fasta"))
