#!/usr/bin/env python

import pandas as pd

"""
You should run this script if the pipeline fails in parsing the blast. It's very likely you have
a very large NUMT in your genome so all mito reads are being drown
there and the mitogenome is not being assembled independently.
"""

my_names = [
    "qseqid",
    "sseqid",
    "pident",
    "alilength",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "leng_query",
    "s_length",
]

blast_cov = pd.read_csv(
    "contigs.blastn",
    sep="\t",
    names=my_names,
)

# Get the percentage of the subject in the blast aligment
blast_cov["alilength"] * 100 / (blast_cov["s_length"])
blast_cov["%subject_in_match"] = blast_cov["alilength"] * 100 / (blast_cov["s_length"])
# blast_cov.to_csv("parsed_blast_percSubject.txt", index=False, sep="\t")


# sum percentages of subject sequence in blast match based on column query id and subject id
a = (
    blast_cov.groupby(["qseqid", "sseqid"])["%subject_in_match"]
    .sum()
    .to_frame()
    .rename(columns={"sseqid": "%subject_in_match"})
    .reset_index()
)
a.to_csv("parsed_blast_SumPercSubject.txt", index=False, sep="\t")
print(
    "If on your new created file 'parsed_blast_SumPercSubject.txt' you see 100% or close to it of"
    + "your subject in the blast match,"
    + "its very likely you have a large NUMT pulling all the mitoreads to it. \n"
    + "You will need to identify those reads and run an assembler again"
    + "(e.g. hicanu or hifiasm), then come back to MitoHiFi with the mito assembled contigs."
)
