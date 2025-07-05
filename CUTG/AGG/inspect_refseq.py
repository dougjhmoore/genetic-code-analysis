import pandas as pd

df = pd.read_csv(r"AGG\refseq_codon_species.tsv", sep="\t")

print("Shape:", df.shape)
print("\nTop 10 by total codons:")
print(df.nlargest(10, "#Codons")[["Taxid", "Organelle", "#CDS", "#Codons"]])

