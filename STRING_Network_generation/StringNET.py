import pandas as pd
import requests
import networkx as nx
import matplotlib.pyplot as plt

# ==========================
# PARAMETERS
# ==========================
DISEASE = "DISEASE_NAME"  # name of disease for files
GENE_FILE = "(FILENAME)"  # input gene list file
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

interaction_file = f"{DISEASE}_string_interactions.tsv"
with open(interaction_file, "w") as f:
    f.write(response.text)
print(f"STRING interactions saved to {interaction_file}")

# ==========================
# STEP 3: LOAD INTERACTIONS
# ==========================
ppi_df = pd.read_csv(interaction_file, sep="\t")
print(f"Interactions retrieved: {len(ppi_df)}")

# ==========================
# STEP 4: BUILD NETWORK
# ==========================
G = nx.from_pandas_edgelist(
    ppi_df,
    source="preferredName_A",
    target="preferredName_B"
)
print(f"Total nodes: {G.number_of_nodes()}")
print(f"Total edges: {G.number_of_edges()}")

# ==========================
# STEP 5: EXTRACT DISEASE SUBNETWORK
# ==========================
sub_nodes = set()
for gene in gene_list:
    if gene in G:
        sub_nodes.add(gene)
        sub_nodes.update(G.neighbors(gene))

disease_network = G.subgraph(sub_nodes).copy()
print(f"Disease subnetwork nodes: {disease_network.number_of_nodes()}")
print(f"Disease subnetwork edges: {disease_network.number_of_edges()}")

# ==========================
# STEP 6: SAVE NETWORK
# ==========================
edgelist_file = f"{DISEASE}_disease_subnetwork.edgelist"
nx.write_edgelist(disease_network, edgelist_file)
print(f"Disease subnetwork saved to {edgelist_file}")

# ==========================
# STEP 7: VISUALIZE NETWORK WITH ADJUSTED NODE/EDGE SIZE
# ==========================
plt.figure(figsize=(12,12))

# layout
if disease_network.number_of_nodes() < 50:
    pos = nx.spring_layout(disease_network, seed=42, k=0.5)  # spread out small networks
else:
    pos = nx.kamada_kawai_layout(disease_network)  # better for larger networks

# separate connected and isolated nodes
connected_nodes = [n for n in disease_network.nodes() if disease_network.degree(n) > 0]
isolated_nodes = [n for n in disease_network.nodes() if disease_network.degree(n) == 0]

# adjust node size and edge width based on network size
node_size = max(3000 // disease_network.number_of_nodes(), 50)
edge_width = max(200 // disease_network.number_of_edges(), 0.5)

# draw edges first
nx.draw_networkx_edges(disease_network, pos, alpha=0.7, width=edge_width)

# draw connected nodes
nx.draw_networkx_nodes(disease_network, pos, nodelist=connected_nodes, node_size=node_size, node_color="skyblue")

# draw isolated nodes on top
nx.draw_networkx_nodes(disease_network, pos, nodelist=isolated_nodes, node_size=node_size, node_color="lightcoral")

# draw labels
nx.draw_networkx_labels(disease_network, pos, font_size=8)

plt.title(f"{DISEASE} PPI Subnetwork")
plt.axis("off")
plt.tight_layout()

# save figure
fig_file = f"{DISEASE}_disease_subnetwork.png"
plt.savefig(fig_file, dpi=300)
plt.show()
print(f"Network figure saved to {fig_file}")
