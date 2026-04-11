import pandas as pd
import networkx as nx
import random
import numpy as np
import scipy.stats as stats


#load disease network
disease_network = nx.read_edgelist('/Users/mrunmayeewankhede/Downloads/CardioVascularDisease_disease_subnetwork (1).edgelist')

#load centrality scores
df = pd.read_csv("/Users/mrunmayeewankhede/Downloads/Cardiovascular_Disease_centrality_scores_conf700.csv")

#sort hubs
df = df.sort_values(by="degree", ascending=False)

top_hubs = df.head(10)["gene"].tolist()

#degree values for the top hubs
degree_dict = dict(disease_network.degree())

real_avg = np.mean([degree_dict[g] for g in top_hubs])

print("Real avg degree:", real_avg)

#randomization
N = 10000
all_nodes = list(disease_network.nodes())

random_avgs = []

for _ in range(N):
    rand_genes = random.sample(all_nodes, 10)
    avg = np.mean([degree_dict[g] for g in rand_genes])
    random_avgs.append(avg)

#p-value
count = sum(1 for x in random_avgs if x >= real_avg)
p_value = count / N

print("P-value:", p_value)

#z-score of your real hubs against the null distribution
z_score = (real_avg - np.mean(random_avgs)) / np.std(random_avgs)

#analytical p-value from the normal distribution
analytical_p = stats.norm.sf(z_score)  #one-tailed

print(f"Z-score: {z_score:.4f}")
print(f"Analytical p-value: {analytical_p:.2e}")

"""print("Top hubs:", top_hubs[:5])
print("Sample network nodes:", list(disease_network.nodes())[:5])
missing = [g for g in top_hubs if g not in disease_network]

print("Missing genes:", missing)"""
