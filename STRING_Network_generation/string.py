import pandas as pd
import requests
import networkx as nx

# ==========================
# PARAMETERS (CHAT PLSSSS EDIT THESE)
# ==========================

GENE_FILE = "disease_genes.csv"   # file containing gene symbols
GENE_COLUMN = "symbol"            # column name containing gene symbols
SPECIES_ID = 9606                 # human
CONFIDENCE_SCORE = 700            # STRING high-confidence interactions

# ==========================
# STEP 1: LOAD GENE LIST
# ==========================

genes_df = pd.read_csv(GENE_FILE)

gene_list = genes_df[GENE_COLUMN].dropna().unique().tolist()

print(f"Number of genes loaded: {len(gene_list)}")

# ==========================
# STEP 2: QUERY STRING API
# ==========================

string_api_url = "https://string-db.org/api"
output_format = "tsv"
method = "network"

genes_str = "%0d".join(gene_list)

request_url = "/".join([string_api_url, output_format, method])

params = {
    "identifiers": genes_str,
    "species": SPECIES_ID,
    "required_score": CONFIDENCE_SCORE
}

print("Querying STRING database...")

response = requests.post(request_url, data=params)

if response.status_code != 200:
    raise Exception("STRING API request failed")

# Save raw interaction data
with open("string_interactions.tsv", "w") as f:
    f.write(response.text)

print("STRING interactions saved.")

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

print(f"Total nodes in network: {G.number_of_nodes()}")
print(f"Total edges in network: {G.number_of_edges()}")

# ==========================
# STEP 5: EXTRACT DISEASE SUBNETWORK
# ==========================

disease_genes = set(gene_list)

sub_nodes = [
    node for node in G.nodes()
    if node in disease_genes
]

disease_network = G.subgraph(sub_nodes)

print(f"Disease subnetwork nodes: {disease_network.number_of_nodes()}")
print(f"Disease subnetwork edges: {disease_network.number_of_edges()}")

# ==========================
# STEP 6: SAVE NETWORK
# ==========================

nx.write_edgelist(disease_network, "disease_subnetwork.edgelist")

print("Disease subnetwork saved.")
