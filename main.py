import csv
import time
from Bio import Entrez, SeqIO
import xml.etree.ElementTree as ET
from Bio.Seq import Seq

# Read gene names from a CSV file
gene_names = []

with open("gene_names.csv", "r") as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        gene_names.append(row["Gene_Name"])
print(gene_names)

# From Biopython: Always tell Entrez who you are!
Entrez.email = "izavel.leesy01@gmail.com"

# Define the genome accession number for E. coli (RefSeq ID)
ecoli_accession = "NC_000913.3"

# Read the E. coli genome sequence from a FASTA file
with open("sequence.fasta", "r") as fasta_file:
    ecoli_genome = "".join(line.strip() for line in fasta_file if not line.startswith(">"))

# Remove whitespace from the E. coli genome sequence
ecoli_genome = ecoli_genome.replace(" ", "")

# Initialize a list to store the DNA nucleotide sequences
nucleotide_sequences = []

for gene_name in gene_names:
    # Use the Entrez.esearch function to search for the specific gene in the E. coli genome
    search_term = f"{ecoli_accession}[Accession] AND {gene_name}[Gene Name]"
    search_handle = Entrez.esearch(db="gene", term=search_term)
    search_results = Entrez.read(search_handle)
    search_handle.close()

    if search_results["Count"] == "0":
        print(f"No genes found for the name '{gene_name}' in E. coli genome ({ecoli_accession})")
    else:
        gene_id = search_results["IdList"][0]
        print(gene_id)
        # Retrieve the DNA nucleotide sequence for the gene only
        handle = Entrez.efetch(db="gene", id=gene_id, rettype="docsum")
        # Read the XML data
        xml_data = handle.read()
        handle.close()
        # Retrieve GenBank record for the gene
        handle = Entrez.efetch(db="gene", id=gene_id, rettype="gb", retmode="text")
        if handle:
            genbank_record = SeqIO.read(handle, "genbank")
            handle.close()
        else:
            print("No GenBank record found for the given Gene ID:", gene_id)
        # Parse XML data
        root = ET.fromstring(xml_data)

        # Retrieve the ChrStart and ChrStop elements
        chr_1 = int(root.find(".//ChrStart").text)
        chr_2 = int(root.find(".//ChrStop").text) + 1

        if chr_2 > chr_1:
            chr_start = chr_1
            chr_stop = chr_2
            Take_Complement = False
        elif chr_1 > chr_2:
            chr_start = chr_2
            chr_stop = chr_1
            Take_Complement = False

        print(f"This is chr_start for '{gene_name}': ", str(chr_start))
        print(f"This is chr_stop for '{gene_name}': ", str(chr_stop))

        # Check strand information
        for feature in genbank_record.features:
            if "strand" in feature.qualifiers:
                strand_info = feature.qualifiers["strand"][0]

        if strand_info:
            # Determine if sequence is in forward or reverse strand
            if strand_info == "1" or strand_info == "+":
                strand = "Forward"
            elif strand_info == "-1" or strand_info == "-":
                strand = "Reverse"
            else:
                strand = "Unknown"

        print(strand)

        # Process the sequence from the genome file
        sliced_sequence = ecoli_genome[chr_start:chr_stop]

        if Take_Complement and strand == "Forward":
            sliced_sequence = Seq(sliced_sequence)
            nucleotide_sequence = str(sliced_sequence.complement())
        elif Take_Complement and strand == "Reverse":
            sliced_sequence = Seq(sliced_sequence)
            nucleotide_sequence = str(sliced_sequence.reverse_complement())
        else:
            nucleotide_sequence = sliced_sequence

        nucleotide_sequences.append([gene_name, str(nucleotide_sequence)])

    # Sleep for 1/3 of a second to respect Entrez's rate limit (3 queries per second)
    time.sleep(1/3)

# Save the DNA nucleotide sequences to a CSV file
with open("nucleotide_sequences.csv", "w", newline="") as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(["Gene_Name", "sliced_sequence"])
    csvwriter.writerows(nucleotide_sequences)

# cspC sequence can only be retrieved if reverse complement instead of complement - but why?
# acs works with just complement. check strand?