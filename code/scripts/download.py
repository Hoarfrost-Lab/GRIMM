#CREDIT: Soumya
#Download sequences for each embl cds id

import pandas as pd
import requests
import certifi
import contextlib
from tqdm import tqdm
import numpy as np
from multiprocessing import Pool

def download_batch(ids):
    url_root = 'https://www.ebi.ac.uk/ena/browser/api/fasta/'
    url = url_root + ','.join(ids)
    try:
        response = requests.get(url, verify=certifi.where()).text
        return response
    except requests.exceptions.RequestException as e:
        print(f"Error making request for {','.join(ids)}: {e}")
        return ""
    except IOError as e:
        print(f"Error writing to file: {e}")
        return ""

def get_seqs(ids, output_file):
    batch_size = 500  # Adjust batch size as needed
    starts = list(np.arange(0, len(ids), batch_size))
    stops = list(np.arange(batch_size, len(ids), batch_size)) + [len(ids)]

    with open(output_file, 'a') as f:
        with contextlib.closing(f):
            with Pool(processes=4) as pool:  # Adjust the number of processes as needed
                results = list(tqdm(pool.imap(download_batch, [ids[start:stop] for start, stop in zip(starts, stops)]), total=len(starts)))

            for result in results:
                f.write(result)

# Usage
df = pd.read_csv('./swissprot_full_w_cds_ids.tsv', header=0, sep='\t')
df = df.dropna(subset=['EMBL_CDS'])

_ids = pd.Series(list(';'.join(df['EMBL_CDS'].tolist()).split(';')))

ids = (
    _ids
    .dropna()  # Remove NaNs
    .astype(str)  # Convert to string
    .str.strip()  # Remove any leading/trailing whitespace
)
get_seqs(ids, output_file= 'nuc_sequences.fasta')
