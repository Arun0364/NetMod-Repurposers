import pandas as pd
import requests
import networkx as nx

# ==========================
# PARAMETERS
# ==========================

GENE_FILE ='use your path name chat'
#main change (from symbol to gene)
GENE_COLUMN = "gene"

SPECIES_ID = 9606
CONFIDENCE_SCORE = 700

# ==========================
# STEP 1: LOAD GENE LIST
# ==========================

genes_df = pd.read_csv(GENE_FILE)

gene_list = genes_df[GENE_COLUMN].dropna().unique().tolist()

print(f"Number of genes loaded: {len(gene_list)}")

# ==========================
# STEP 2: QUERY STRING API
# ==========================

print("Querying STRING database...")

string_api_url = "https://string-db.org/api"
output_format = "tsv"
method = "network"

request_url = "/".join([string_api_url, output_format, method])

genes_str = "%0d".join(gene_list)

params = {
    "identifiers": genes_str,
    "species": SPECIES_ID,
    "required_score": CONFIDENCE_SCORE
}

response = requests.post(request_url, data=params)

if response.status_code != 200:
    raise Exception("STRING API request failed")

with open("string_interactions.tsv", "w") as f:
    f.write(response.text)

print("STRING interactions saved")

# ==========================
# STEP 3: LOAD INTERACTIONS
# ==========================

ppi_df = pd.read_csv("string_interactions.tsv", sep="\t")

print(f"Interactions retrieved: {len(ppi_df)}")

# ==========================
# STEP 4: BUILD NETWORK
# ==========================

G = nx.from_pandas_edgelist(
    ppi_df,
    source="preferredName_A",
    target="preferredName_B"
)

print("Total nodes:", G.number_of_nodes())
print("Total edges:", G.number_of_edges())

# ==========================
# STEP 5: EXTRACT DISEASE SUBNETWORK
# ==========================

disease_genes = set(gene_list)

sub_nodes = set()

for gene in disease_genes:
    
    if gene in G:
        
        sub_nodes.add(gene)
        
        neighbors = list(G.neighbors(gene))
        
        sub_nodes.update(neighbors)

disease_network = G.subgraph(sub_nodes)

print("Disease subnetwork nodes:", disease_network.number_of_nodes())
print("Disease subnetwork edges:", disease_network.number_of_edges())

# ==========================
# STEP 6: SAVE NETWORK
# ==========================

nx.write_edgelist(disease_network, "disease_subnetwork.edgelist")

print("Disease subnetwork saved.")
