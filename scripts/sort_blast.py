import pandas as pd
import sys


my_names = ["qseqid", "sseqid", "pident", "alilength" , "mismatch", "gapopen", "qstart", "qend", 
            "sstart", "send", "evalue" , "bitscore", "leng_query", "s_length",]

blast_cov = pd.read_csv("contigs.blastn", 
                        sep="\t", names = my_names, )

#Get the percentage of the query in the blast aligment
blast_cov['alilength']*100 / (blast_cov['leng_query'])
blast_cov['%q_in_match'] = blast_cov['alilength']*100 / (blast_cov['leng_query'])

#sum percentages of query sequence in blast match based on column id
a= blast_cov.groupby('qseqid')['%q_in_match'].sum().to_frame().rename(columns={'qseqid':'%q_in_match'}).reset_index()

#get size of query and subject and drop duplicates
seqsizes = blast_cov[['qseqid', 'leng_query', 's_length']].drop_duplicates(subset='qseqid')

#merge 'a' and 'seqsizes' dataframes by 'qseqid'
result = pd.merge(a, seqsizes, on='qseqid')

# Now let's filter the blast matches
# if the lenght of the query is 5x the size of the subject (close-related mitogenome), drop it. (As its likely the match belongs to a NUMT)
five_times = (result['s_length'] * 5)
result1 = result[(result['leng_query'] < five_times)].sort_values(by='%q_in_match', ascending=False)

# if the lenght of the query is smaller than the length of the subject, drop it. Unlikely you will have a complete mitogenome.
slen=result1['s_length']
ac=result1[result1['leng_query'] > slen].sort_values(by='%q_in_match')

# if the % of the query in the blast match is smaller than 70%, drop it
ac[(ac['%q_in_match'] > 70)].sort_values(by='%q_in_match', ascending=False)
ac[(ac['%q_in_match'] > 70)].sort_values(by='%q_in_match', ascending=False).to_csv("/lustre/scratch116/vr/projects/vgp/user/mu2/mito-automation/test-pipe/again3/MitoHiFi/exampleFiles/exaploblastedited", index=False, sep="\t")

id = (ac[(ac['%q_in_match'] > 70)].sort_values(by='%q_in_match', ascending=False)['qseqid'].iloc[0])
id_series = pd.Series(id)
id_series
id_series.to_csv("/lustre/scratch116/vr/projects/vgp/user/mu2/mito-automation/test-pipe/again3/MitoHiFi/exampleFiles/potentialmito.id", index=False, header=None)
print("all done here, let's filter the fasta and them circularise it")
