from Bio import SeqIO
import sys
input_fasta = sys.argv[1]
circu_info = open('circucheck', 'r')
file = circu_info.read()
obj=file.split(",")
if obj[1] == " True":
    record=SeqIO.read(input_fasta, "fasta")
    id = record.id
    get= record[int(obj[3]):]
    print(get.format('fasta'))
