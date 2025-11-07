# Mapping gene names to UniProt IDs
gene_to_uniprot = {
    "TP53": "P04637",
    "BRCA1": "P38398",
    "BRCA2": "P51587",
    "CFTR": "P13569",
    "EGFR": "P00533",
    "HBB": "P68871",
    "MYH7": "P12883",
    "BRAF": "P15056",
    "KRAS": "P01116",
    "SARS-CoV-2_Spike": "P0DTC2",
    "PTEN": "P60484",
    "APC": "P25054"
}

# Original mutation list
mutations = [
    "TP53 R273H",
    "TP53 R175H",
    "BRCA1 185delAG",
    "BRCA1 C61G",
    "BRCA1 5382insC",
    "BRCA2 6174delT",
    "CFTR F508del",
    "EGFR L858R",
    "HBB E6V",
    "MYH7 R403Q",
    "BRAF V600E",
    "KRAS G12D",
    "SARS-CoV-2_Spike N501Y",
    "PTEN R130Q",
    "APC R1450*"
]

# Convert gene names to UniProt ID format
formatted_mutations = []
for m in mutations:
    gene, mutation = m.split()
    uniprot_id = gene_to_uniprot.get(gene, gene)  # fallback to gene if not found
    formatted_mutations.append(f"{uniprot_id} {mutation}")

# Save to file
file_path = "./mutations_list.txt"
with open(file_path, "w") as f:
    for m in formatted_mutations:
        f.write(m + "\n")

print(f"File saved to: {file_path}")
