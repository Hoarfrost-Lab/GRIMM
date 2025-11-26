import pandas as pd
import json

df = pd.read_csv('./swissprot_full.tsv', sep='\t', header=0)

embl_cds_mapping = {}
for chunk in range(6):
    with open('nuc_chunk_{}.json'.format(chunk), 'r') as f:
        d = json.load(f)
        r = d['results']

        for entry in r:
            if entry['from'] not in embl_cds_mapping:
                embl_cds_mapping[entry['from']] = [entry['to']]
            else:
                embl_cds_mapping[entry['from']].append(entry['to'])

        #some may not map
        failed = d['failedIds']
        for entry in failed:
            if entry not in embl_cds_mapping:
                embl_cds_mapping[entry] = []


l1 = []
l2 = []
for key, val in embl_cds_mapping.items():
    l1.append(key)
    l2.append(';'.join(val)) #convert to string to store

parsed_dict = {'Entry' : l1, 'EMBL_CDS' : l2}
df2 = pd.DataFrame.from_dict(parsed_dict)
df = pd.merge(df, df2, on='Entry')

print(df)
df.to_csv('swissprot_full_w_cds_ids.tsv', index=False, sep='\t')

