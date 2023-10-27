import csv
import time
from Bio import Entrez, SeqIO

# Read gene names from a CSV file
gene_names = []

with open("gene_names.csv", "r") as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        gene_names.append(row["Gene_Name"])

# From Biopython: Always tell Entrez who you are!
Entrez.email = "izavel.leesy01@gmail.com"

# Define the genome accession number for E. coli (RefSeq ID)
ecoli_accession = "NC_000913.3"

for gene_name in gene_names:
    # Use the Entrez.esearch function to search for the specific gene in the E. coli genome
    search_term = f"{ecoli_accession}[Accession] AND {gene_name}[Gene Name]"
    search_handle = Entrez.esearch(db="gene", term=search_term)
    search_results = Entrez.read(search_handle)
    search_handle.close()

    if search_results["Count"] == "0":
        print(f"No genes found for the name '{gene_name}' in E. coli genome ({ecoli_accession})")
    else:
        # Retrieve the GenBank record for the gene
        gene_id = search_results["IdList"][0]
        handle = Entrez.efetch(db="gene", id=gene_id, rettype="gb", retmode="text")
        genbank_record = SeqIO.read(handle, "genbank")
        handle.close()

        # Initialize an empty variable to store the DNA sequence
        dna_sequences = ""

        # Search for the "Origin" feature in the GenBank record
        for feature in genbank_record.features:
            if feature.type == "Origin":
                # Extract the DNA sequence from the feature location
                origin_location = feature.location
                dna_sequence = origin_location.extract(genbank_record.seq)
                break  # Stop searching after finding the first "Origin" feature

        if dna_sequence:
            # Append the gene name and DNA sequence to the list
            dna_sequences.append([gene_name, str(dna_sequence)])
        else:
            print(f"No 'Origin' feature found for gene '{gene_name}'")

    # Sleep for 1/3 of a second to respect Entrez's rate limit (3 queries per second)
    time.sleep(1/3)

# Save the mRNA transcript sequences to a CSV file
with open("dna_sequences.csv", "w", newline="") as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(["Gene_Name", "dna_sequence"])
    csvwriter.writerows(dna_sequences)
