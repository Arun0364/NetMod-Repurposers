import pandas as pd
import mygene

#load DisGeNET file
df = pd.read_csv("disgenet_data.tsv", sep="\t")

print(df.columns)

#filter genes by association score
df_filtered = df[df["score"] >= 0.3].copy()

print("Associated genes:", len(df_filtered))

#extract gene symbols
gene_list = df_filtered["geneSymbol"].dropna().unique()

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

#save data for future reference
final_genes.to_csv("disgenet_disease_genes.csv", index=False)
