import pandas as pd
from Bio import SeqIO

def read_fasta(fasta_path):
    seq_ids = []
    seqs = []

    for seq_record in SeqIO.parse(fasta_path, "fasta"):
        seq_ids.append(str(seq_record.id).split('|')[-1])
        seqs.append(str(seq_record.seq))

    assert(len(seq_ids) == len(seqs))
    return seq_ids, seqs

print('Loading files...')
directory = './custom/split_{}/{}'
new_directory = './custom/nuc_split_{}/{}'
file_list = ['train_reformated.csv', 'valid_reformated.csv', 'test1_reformated.csv', 'test2_reformated.csv']

embl_cds_mapping = pd.read_csv('./swissprot_full_w_cds_ids.tsv', header=0, sep='\t').fillna('').set_index('Entry').to_dict()['EMBL_CDS']

nuc_ids, nuc_seqs = read_fasta('./nuc_sequences.fasta')
seq_mapping = pd.DataFrame({'Entry' : nuc_ids, 'Sequence' : nuc_seqs}).set_index('Entry').to_dict()['Sequence']

for i in range(1, 6):
    print('Creating split {}...'.format(i))
    entries = []
    ecs = []
    seqs = []
    for filename in file_list:
        df_aa = pd.read_csv(directory.format(i, filename), header=0, sep='\t')

        for entry, ec in zip(df_aa['Entry'], df_aa['EC number']):
            entry_list = embl_cds_mapping[entry]

            if not entry_list: #empty string implies no mapping for that id
                continue

            entry_list = entry_list.split(';')
            for new_entry in entry_list:
                if new_entry in seq_mapping.keys() and new_entry not in entries:
                    entries.append(new_entry)
                    ecs.append(ec)
                    seqs.append(seq_mapping[new_entry])

        df_nuc = pd.DataFrame({'Entry' : entries, 'EC number' : ecs, 'Sequence' : seqs})
        df_nuc.to_csv(new_directory.format(i, filename), index=False, sep='\t')        
