# pip install python-louvain for community detection
import pandas as pd
import requests
import networkx as nx
import community as community_louvain
import matplotlib.pyplot as plt
from io import StringIO


DISEASE = "Arthritis"
GENE_FILE = "/Users/arun/Desktop/02604/Project/NetMod-Repurposers/Filtered_Genes/Copy of top_300_genes_arthritis.csv"
GENE_COLUMN = "gene"

SPECIES_ID = 9606

# Computing sensitivity parameters
CONFIDENCE_SCORES = [700] #if not working then [700, 900]]
RESOLUTIONS = [1.0] #if not working then [0.5, 1.0, 1.5]

# Output prefix for files
OUTPUT_PREFIX = DISEASE.replace(" ", "_")

# Loading gene list
genes_df = pd.read_csv(GENE_FILE)
gene_list = genes_df[GENE_COLUMN].dropna().unique().tolist()

print(f"Loaded {len(gene_list)} genes")

# Fetching STRING PPI network for gene list and confidence score.
def get_network(gene_list, conf_score):
    string_api_url = "https://string-db.org/api"
    output_format = "tsv"
    method = "network"
    request_url = "/".join([string_api_url, output_format, method])

    genes_str = "%0d".join(gene_list)

    params = {
        "identifiers": genes_str,
        "species": SPECIES_ID,
        "required_score": conf_score
    }

    response = requests.post(request_url, data=params)

    if response.status_code != 200:
        raise Exception(f"STRING API request failed with status code {response.status_code}")

    ppi_df = pd.read_csv(StringIO(response.text), sep="\t")

    if ppi_df.empty:
        raise Exception("No interactions returned from STRING.")
    
    #Checking available edge weight column
    edge_attr_col = None
    if "score" in ppi_df.columns:
        edge_attr_col = "score"
    elif "combined_score" in ppi_df.columns:
        edge_attr_col = "combined_score"

    if edge_attr_col:
        G = nx.from_pandas_edgelist(
            ppi_df,
            source="preferredName_A",
            target="preferredName_B",
            edge_attr=edge_attr_col
    )
    else:
        G = nx.from_pandas_edgelist(
            ppi_df,
            source="preferredName_A",
            target="preferredName_B"
    ) 

    return G, ppi_df

# Get largest connected component as subgraph
def get_largest_component(G):
    if G.number_of_nodes() == 0:
        raise Exception("Graph is empty.")
    
    largest_cc = max(nx.connected_components(G), key=len)
    G_main = G.subgraph(largest_cc).copy()
    return G_main

# Computing centrality measures (Degree, Betweenness, Closeness)
def compute_centralities(G):
    degree = nx.degree_centrality(G)
    betweenness = nx.betweenness_centrality(G)
    closeness = nx.closeness_centrality(G)

    centrality_df = pd.DataFrame({
        "gene": list(G.nodes()),
        "degree": [degree[g] for g in G.nodes()],
        "betweenness": [betweenness[g] for g in G.nodes()],
        "closeness": [closeness[g] for g in G.nodes()]
    })

    # Combined score = Average of three centralities
    centrality_df["combined_score"] = (
        centrality_df["degree"] + 
        centrality_df["betweenness"] +
        centrality_df["closeness"]
    ) / 3

    centrality_df = centrality_df.sort_values("combined_score", ascending=False).reset_index(drop=True)
    return centrality_df

# Running Louvain community detection
def detect_communities(G, resolution=1.0):
    partition = community_louvain.best_partition(G, resolution=resolution)
    return partition

# Building per-community summary
def summarize_modules(gene_metrics_df):
    module_summary = []

    for comm in sorted(gene_metrics_df["community"].unique()):
        sub = gene_metrics_df[gene_metrics_df["community"] == comm].copy()
        sub = sub.sort_values("combined_score", ascending=False)

        top_genes = sub["gene"].head(5).tolist()

        module_summary.append({
            "community": comm,
            "num_genes": len(sub),
            "top_5_genes": ", ".join(top_genes)
        })

    module_summary_df = pd.DataFrame(module_summary).sort_values("num_genes", ascending=False).reset_index(drop=True)
    return module_summary_df

# Plotting network with communities colored, sized by centrality
def plot_network(G, partition, degree_dict, out_file, title):
    plt.figure(figsize=(14, 12))
    pos = nx.spring_layout(G, seed=42)

    node_colors = [partition[node] for node in G.nodes()]
    node_sizes = [3000 * degree_dict[node] + 100 for node in G.nodes()]

    nx.draw_networkx_edges(G, pos, alpha=0.3, edge_color="gray")
    nx.draw_networkx_nodes(
        G,
        pos,
        node_color=node_colors,
        node_size=node_sizes,
        cmap=plt.cm.tab20
    )
    nx.draw_networkx_labels(G, pos, font_size=8)
    
    plt.title(title)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close()

results = []

for conf in CONFIDENCE_SCORES:
   
    print(f"\n--- Confidence Score: {conf} ---")

    # Step 1: Fetching STRING network
    G_full, ppi_df = get_network(gene_list, conf)

    print(f"Full graph -> Nodes: {G_full.number_of_nodes()}, Edges: {G_full.number_of_edges()}")

    # Saving raw PPI table
    ppi_outfile = f"{OUTPUT_PREFIX}_STRING_PPI_conf{conf}.csv"
    ppi_df.to_csv(ppi_outfile, index=False)
    print(f"STRING PPI data saved to {ppi_outfile}")

    # Step 2: Using largest connected component
    G = get_largest_component(G_full)
    print(f"Largest connected component -> Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")

    # Step 3: Computing centrality measures
    centrality_df = compute_centralities(G)

    centrality_outfile = f"{OUTPUT_PREFIX}_centrality_scores_conf{conf}.csv"
    centrality_df.to_csv(centrality_outfile, index=False)
    print(f"Saved centrality scores: {centrality_outfile}")

    # Step 4: Getting top hubs
    top_hubs = centrality_df["gene"].head(10).tolist()
    print("Top 10 hubs:", top_hubs)
   
    for res in RESOLUTIONS:
        print(f"\n--- Louvain resolution: {res} ---")

        # Step 5: Comunity detection
        partition = detect_communities(G, resolution=res)
        num_communities = len(set(partition.values()))
        print(f"Detected communities: {num_communities}")

        community_df = pd.DataFrame({
            "gene": list(partition.keys()),
            "community": list(partition.values())
        })

        # Step 6: Merging centrality with community
        gene_metrics_df = community_df.merge(centrality_df, on="gene", how="left")

        metrics_outfile = f"{OUTPUT_PREFIX}_gene_metrics_conf{conf}_res{res}.csv"
        gene_metrics_df.to_csv(metrics_outfile, index=False)
        print(f"Saved gene metrics with communities: {metrics_outfile}")

        # Step 7: Module summary
        module_summary_df = summarize_modules(gene_metrics_df)

        module_outfile = f"{OUTPUT_PREFIX}_module_summary_conf{conf}_res{res}.csv"
        module_summary_df.to_csv(module_outfile, index=False)
        print(f"Saved module summary: {module_outfile}")

        # Step 8: Plotting network
        degree_dict = nx.degree_centrality(G)
        plot_outfile = f"{OUTPUT_PREFIX}_network_conf{conf}_res{res}.png"
        plot_title = f"{DISEASE} PPI Network | STRING confidence={conf}, Louvain resolution={res}"
        plot_network(G, partition, degree_dict, plot_outfile, plot_title)
        print(f"Saved network plot: {plot_outfile}")

        # Step 9: Saving summary row
        results.append({
            "confidence": conf,
            "resolution": res,
            "nodes_full_graph": G_full.number_of_nodes(),
            "edges_full_graph": G_full.number_of_edges(),
            "nodes_largest_component": G.number_of_nodes(),
            "edges_largest_component": G.number_of_edges(),
            "num_communities": num_communities,
            "top_10_hubs": ", ".join(top_hubs)
        })

# Step 10: Saving overall sensitivity summary
results_df = pd.DataFrame(results)
summary_outfile = f"{OUTPUT_PREFIX}_sensitivity_analysis_summary.csv"
results_df.to_csv(summary_outfile, index=False)

print(f"\nSaved overall sensitivity analysis summary: {summary_outfile}")