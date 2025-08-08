import requests
import sys
import xml.etree.ElementTree as ET

# Usage: python ncbi_fetcher.py "protein_name" "organism" <api_key> [max_records]
# E.g.: python3 xtractseq_ncbi.py "prion protein" "Bos taurus" abb4f7cff84a4af777891b6f35184e703000

def search_protein(protein_name: str, organism: str, api_key: str, max_records: int = 400):
    """
    Search NCBI Protein database for given protein name and organism.
    Returns a list of protein IDs up to max_records.
    """
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        'db': 'protein',
        'term': f'{protein_name}[Protein Name] AND "{organism}"[Organism] NOT partial[All Fields]',
        'api_key': api_key,
        'retmode': 'xml',
        'retmax': str(max_records)
    }
    resp = requests.get(base_url, params=params)
    resp.raise_for_status()
    root = ET.fromstring(resp.text)
    ids = [id_elem.text for id_elem in root.findall('.//Id')]
    return ids


def fetch_fasta(ids, api_key: str, batch_size: int = 100):
    """
    Fetch FASTA sequences for given list of IDs, batching to avoid URL length limits.
    Uses POST requests to safely transmit longer ID lists.
    Returns concatenated FASTA text.
    """
    if not ids:
        return ''
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    fasta_chunks = []
    # Batch list in groups of batch_size
    for i in range(0, len(ids), batch_size):
        batch = ids[i:i+batch_size]
        # Use POST to avoid URL length limits
        data = {
            'db': 'protein',
            'id': ','.join(batch),
            'retmode': 'text',
            'rettype': 'fasta',
            'api_key': api_key
        }
        resp = requests.post(base_url, data=data)
        resp.raise_for_status()
        fasta_chunks.append(resp.text)
    return "\n".join(fasta_chunks)


def filter_by_prefix_and_organism(fasta_data: str, max_per_prefix: int, required_organism: str) -> str:
    """
    Keep only records with the exact organism name, and limit to N records per 5-letter accession prefix.
    """
    if not fasta_data:
        return ''
    records = fasta_data.strip().split('\n>')
    filtered = []
    counts = {}
#    prec = 0
    orgnomatch = 0
    for rec in records:
        header, *seq_lines = rec.split('\n')
        acc = header.split()[0]
        prefix = acc[:5]

#        if 'precursor' in header.lower():
#            prec+=1
#            continue   # skip 'precursor' sequences

        # Extract organism in brackets at end of header
        if '[' in header and ']' in header:
            org_in_brackets = header.split('[')[-1].rstrip(']')
        else:
            continue  # skip malformed headers

        if org_in_brackets != required_organism:
            orgnomatch+=1
            print(org_in_brackets)
            continue  # skip if organism doesn't match exactly

        counts.setdefault(prefix, 0)
        if counts[prefix] < max_per_prefix:
            counts[prefix] += 1
            filtered.append(('>' + rec) if not rec.startswith('>') else rec)

    print(f"Total input records: {len(records)}")
#    print(f"Excluded 'precursor': {prec}")
    print(f"Excluded wrong organism: {orgnomatch}")
    print(f"Final kept records: {len(filtered)}")


    return '\n'.join(filtered)


def main():
    if len(sys.argv) < 4 or len(sys.argv) > 5:
        print(f"Usage: python {sys.argv[0]} <protein_name> <organism> <api_key> [max_records]")
        sys.exit(1)

    protein_name = sys.argv[1]
    organism = sys.argv[2]
    api_key = sys.argv[3]
    max_records = int(sys.argv[4]) if len(sys.argv) == 5 else 400

    print(f"Searching for '{protein_name}' in '{organism}' (max {max_records} records)...")
    ids = search_protein(protein_name, organism, api_key, max_records)
    if not ids:
        print("No sequences found.")
        return

    print(f"Found {len(ids)} IDs. Fetching FASTA in batches...")
    raw_fasta = fetch_fasta(ids, api_key, batch_size=100)
    if not raw_fasta:
        print("No FASTA data returned.")
        return

    print("Filtering sequences by 5-letter prefix (max 15 each)...")
#    filtered_fasta = filter_by_prefix(raw_fasta, max_per_prefix=15)

    filtered_fasta = filter_by_prefix_and_organism(
        raw_fasta,
        max_per_prefix=15,
        required_organism=organism
    )

#    output_file = f"{protein_name.replace(' ', '_')}_{organism.replace(' ', '_')}.fasta"
#    with open(output_file, 'w') as f:
#        f.write(filtered_fasta)

    # --- Save only the first 20 sequences
    records = filtered_fasta.strip().split('\n>')
    top_20 = records[:20]
    top_20_text = '\n>'.join(top_20)
    if not top_20_text.startswith('>'):
        top_20_text = '>' + top_20_text

    output_20_file = f"{protein_name.replace(' ', '_')}_{organism.replace(' ', '_')}_20.fasta"
    with open(output_20_file, 'w') as f20:
        f20.write(top_20_text)

    print(f"Filtered FASTA saved to {output_20_file}")


if __name__ == '__main__':
    main()
