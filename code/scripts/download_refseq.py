# Option 1: Using ncbi-genome-download (recommended for larger datasets)
# Install ncbi-genome-download: pip install ncbi-genome-download
# Example: Download all bacterial RefSeq genomes
# !ncbi-genome-download --genera bacteria --format fasta bacteria

# Option 2: Using BioPython and Entrez (more control, suitable for smaller datasets)
from Bio import Entrez
from Bio import SeqIO

def download_refseq_sequences(accession_ids, email, out_format="fasta"):
    Entrez.email = email  # Always set your email for NCBI
    handle = Entrez.efetch(db="nucleotide", id=accession_ids, rettype="fasta", retmode="text")
    records = SeqIO.parse(handle, "fasta")
    SeqIO.write(records, "output.fasta", out_format)
    handle.close()

# Example usage (replace with your desired accessions and email)
accession_list = ["NM_001195576.2", "NM_001317617.1"]
download_refseq_sequences(accession_list, "your_email@example.com", "gbk")
