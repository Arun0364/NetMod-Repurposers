import pandas as pd
import mygene

#load GWAS data
df = pd.read_csv("gwas_data.tsv", sep="\t")

#convert p-values
df["P-VALUE"] = pd.to_numeric(df["P-VALUE"], errors="coerce")

#filter significant SNPs
df_sig = df[df["P-VALUE"] < 5e-8].copy()

print("Significant SNPs:", len(df_sig))

#clean mapped genes
df_sig["MAPPED_GENE"] = df_sig["MAPPED_GENE"].str.replace(" - ", ",")
df_sig["MAPPED_GENE"] = df_sig["MAPPED_GENE"].str.split(",")

#split rows
df_genes = df_sig.explode("MAPPED_GENE")

#clean whitespace
df_genes["MAPPED_GENE"] = df_genes["MAPPED_GENE"].str.strip()

#remove invalid entries
df_genes = df_genes[df_genes["MAPPED_GENE"].notna()]
df_genes = df_genes[df_genes["MAPPED_GENE"] != "NR"]

#unique genes
gene_list = df_genes["MAPPED_GENE"].unique()

print("Unique genes:", len(gene_list))

#standardize gene identifiers
mg = mygene.MyGeneInfo()

result = mg.querymany(
    list(gene_list),
    scopes="symbol",
    fields="entrezgene,ensembl.gene,symbol",
    species="human"
)

mapped_df = pd.DataFrame(result)
mapped_df = mapped_df[~mapped_df["notfound"].fillna(False)]

final_genes = mapped_df[["symbol", "entrezgene", "ensembl"]].drop_duplicates()

print("Final standardized genes:", len(final_genes))

#save final genes list for future reference
final_genes.to_csv("gwas_disease_genes.csv", index=False)
