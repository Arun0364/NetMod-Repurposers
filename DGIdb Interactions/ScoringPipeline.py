import pandas as pd
import numpy as np
from pathlib import Path


def minmax(series):
    series = series.fillna(0)
    if series.nunique() <= 1:
        return pd.Series(np.zeros(len(series)), index=series.index)
    return (series - series.min()) / (series.max() - series.min())


def compute_module_coverage(drug_gene_df):
    # size of each module
    module_sizes = (
        drug_gene_df.groupby("community")["gene"]
        .nunique()
        .reset_index(name="module_size")
    )

    # targets per (drug, community)
    drug_module_hits = (
        drug_gene_df.groupby(["drug", "community"])["gene"]
        .nunique()
        .reset_index(name="targets_in_module")
        .merge(module_sizes, on="community", how="left")
    )

    drug_module_hits["module_coverage"] = (
        drug_module_hits["targets_in_module"] / drug_module_hits["module_size"]
    )

    module_cov = (
        drug_module_hits.groupby("drug")
        .agg(
            n_communities_hit=("community", "nunique"),
            max_module_coverage=("module_coverage", "max"),
            mean_module_coverage=("module_coverage", "mean")
        )
        .reset_index()
    )

    return module_cov, drug_module_hits


def score_drugs_from_interactions(drug_gene_file, output_dir, output_prefix):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) Load your already-enriched drug-gene table
    df = pd.read_csv(drug_gene_file)

    required = [
        "drug", "gene", "community",
        "degree", "betweenness", "closeness",
        "combined_score", "is_hub"
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in {drug_gene_file}: {missing}")

    # 2) Drop duplicates at (drug, gene) level
    df = df.drop_duplicates(subset=["drug", "gene"])

    # 3) Compute module coverage
    module_cov, drug_module_hits = compute_module_coverage(df)

    # 4) Aggregate per-drug features
    drug_scores = (
        df.groupby("drug")
        .agg(
            n_targets_in_network=("gene", "nunique"),
            n_hub_targets=("is_hub", "sum"),
            sum_centrality=("combined_score", "sum"),
            mean_centrality=("combined_score", "mean"),
            max_centrality=("combined_score", "max"),
            mean_degree=("degree", "mean"),
            mean_betweenness=("betweenness", "mean"),
            mean_closeness=("closeness", "mean"),
            mean_interaction_score=("interaction_score", "mean")
        )
        .reset_index()
    )

    # 5) Merge module coverage metrics
    drug_scores = drug_scores.merge(module_cov, on="drug", how="left").fillna(0)

    # 6) Normalize core features
    norm_cols = [
        "n_targets_in_network",
        "n_hub_targets",
        "sum_centrality",
        "mean_centrality",
        "n_communities_hit",
        "max_module_coverage",
        "mean_module_coverage",
        "mean_interaction_score"
    ]

    for col in norm_cols:
        drug_scores[f"{col}_norm"] = minmax(drug_scores[col])

    # 7) Components reflecting your proposal
    drug_scores["hub_component"] = drug_scores["n_hub_targets_norm"]

    drug_scores["centrality_component"] = (
        0.5 * drug_scores["sum_centrality_norm"] +
        0.3 * drug_scores["mean_centrality_norm"] +
        0.2 * drug_scores["n_targets_in_network_norm"]
    )

    drug_scores["module_component"] = (
        0.5 * drug_scores["n_communities_hit_norm"] +
        0.3 * drug_scores["max_module_coverage_norm"] +
        0.2 * drug_scores["mean_module_coverage_norm"]
    )

    # Optional small evidence component
    drug_scores["evidence_component"] = drug_scores["mean_interaction_score_norm"]

    # 8) Final composite score
    drug_scores["drug_score"] = (
        0.35 * drug_scores["hub_component"] +
        0.40 * drug_scores["centrality_component"] +
        0.20 * drug_scores["module_component"] +
        0.05 * drug_scores["evidence_component"]
    )

    drug_scores = drug_scores.sort_values("drug_score", ascending=False).reset_index(drop=True)

    # 9) Save outputs
    scores_file = output_dir / f"{output_prefix}_drug_scores_scoringOnly.csv"
    module_hits_file = output_dir / f"{output_prefix}_drug_module_hits_scoringOnly.csv"

    drug_scores.to_csv(scores_file, index=False)
    drug_module_hits.to_csv(module_hits_file, index=False)

    print(f"Saved ranked drug scores: {scores_file}")
    print(f"Saved drug-module hits: {module_hits_file}")
    print(f"Total ranked drugs: {len(drug_scores)}")


if __name__ == "__main__":
    score_drugs_from_interactions(
        drug_gene_file="/Users/vineetpaliwal/NetMod-Repurposers/DGIdb Interactions/output/Cardiovascular_Disease_drug_gene_interactions.csv",
        output_dir="/Users/vineetpaliwal/NetMod-Repurposers/DGIdb Interactions/output",
        output_prefix="Cardiovascular_Disease"
    )