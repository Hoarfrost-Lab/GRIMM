import pandas as pd
import numpy as np
import json
from Bio import SeqIO

df = pd.read_csv('uniprotkb_AND_reviewed_true_2025_05_23.tsv', header=0, sep='\t')

entries = []
uniref50 = []
uniref90 = []
uniref100 = []

for chunk in range(6):
    with open('chunk_{}.json'.format(chunk), 'r') as f:
        d = json.load(f)
        r = d['results']

        for entry in r:
            entries.append(entry['from'])
            uniref50.append(entry['to']['id'])
            uniref90.append(entry['to']['representativeMember']['uniref90Id'])
            uniref100.append(entry['to']['representativeMember']['uniref100Id'])

df2 = pd.DataFrame()
df2['Entry'] = entries
df2['UniRef50'] = uniref50
df2['UniRef90'] = uniref90
df2['UniRef100'] = uniref100

df_partial = pd.merge(df, df2, on='Entry')

uniref50 = []
seqs = []

fasta_sequences = SeqIO.parse(open('uniref50.fasta'),'fasta')
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    uniref50.append(name)
    seqs.append(sequence)

df3 = pd.DataFrame()
df3['UniRef50'] = uniref50
df3['Sequence'] = seqs

df_full = pd.merge(df_partial, df3, on='UniRef50')

def filter_ec(ecs):
    if pd.isnull(ecs):
        return np.nan

    ec_list = ecs.split(';')
    ec_list = [ec for ec in ec_list if not ec.__contains__('-')]
    
    if len(ec_list) > 1:
        return ';'.join(ec_list).replace(' ', '')
    elif len(ec_list) == 1:
        return ec_list[0]
    else:
        return np.nan

df_full['EC number'] = df_full['EC number'].map(filter_ec) 
filtered_df = df_full[df_full[['UniRef50', 'EC number']].notnull().all(1)]
filtered_df.to_csv('swissprot_full.tsv', index=False, sep='\t')
