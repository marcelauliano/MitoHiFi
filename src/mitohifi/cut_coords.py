#!/usr/bin/env python
# Credit Marcela Uliano-Silva

import sys

from Bio import SeqIO

input_fasta = sys.argv[1]
circu_info = open("circularisationCheck.txt", "r")

file = circu_info.read()
obj = file.split(",")

if obj[1] == " True":
    record = SeqIO.read(input_fasta, "fasta")
    id = record.id
    get = record[int(obj[3]) :]
    print(get.format("fasta"))
else:
    record = SeqIO.read(input_fasta, "fasta")
    id = record.id
    print(record.format("fasta"))
