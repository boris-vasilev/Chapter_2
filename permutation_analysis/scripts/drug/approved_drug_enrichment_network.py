#!/usr/bin/env python
"""Plot a network linking coloc (H4) genes to approved-drug-target hits.

For one (disease, gene_set, condition) triple — the same triple you'd pass to
approved_drug_enrichment_drilldown.py — builds a three-shell graph:

  - inner shell: "coloc" nodes (H4 GWAS genes), gold squares;
  - middle shell: "hit" nodes (gene-set members with approved/Phase II+ drug
    for the focal disease), green gradient by best clinical stage;
  - outer shell: "drug" nodes (deduplicated by Open Targets drugId), light
    gradient by best stage for any target.

Edges:
  - coloc→hit: Perturb-seq DE relationships gated on H4 membership, coloured
    by direction (red = +logFC, blue = -logFC), weighted by |logFC|.
  - hit→drug: thin grey lines, one per (target, drug) pair, showing which
    drug acts on which target.

The visual emphasises how multiple already-drugged genes converge on shared
coloc genes — the qualitative narrative for recovering convergent disease
biology — and lets you see the actual drug repertoire at one glance.

Outputs PNG by default; pass --graphml to also export a Cytoscape-compatible
graph file for interactive exploration.
"""

import argparse
from pathlib import Path

import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import networkx as nx
import polars as pl
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

APPROVED_STAGES = [
    "APPROVAL",
    "PREAPPROVAL",
    "PHASE_4",
    "PHASE_3",
    "PHASE_2_3",
    "PHASE_2",
]
# Strongest evidence first.
STAGE_PRIORITY = [
    "APPROVAL",
    "PHASE_4",
    "PREAPPROVAL",
    "PHASE_3",
    "PHASE_2_3",
    "PHASE_2",
]


def expand_subtree(disease_ids: list[str], ontology: pl.DataFrame) -> set[str]:
    matches = ontology.filter(pl.col("id").is_in(disease_ids))
    subtree = set(disease_ids)
    for d in matches.get_column("descendants").to_list():
        if d:
            subtree.update(d)
    return subtree


def load_gene_set(
    gene_sets_dir: Path, studies: list[str], gene_set: str, condition: str
) -> set[str]:
    genes: set[str] = set()
    for study in studies:
        for file in gene_sets_dir.glob(f"{study}_*_{condition}_{gene_set}.txt"):
            with open(file) as fh:
                for line in fh:
                    g = line.strip()
                    if g:
                        genes.add(g)
    return genes


def pick_max_stage(stages: list[str]) -> str:
    for s in STAGE_PRIORITY:
        if s in stages:
            return s
    return ""


def stage_rank(stage: str) -> int:
    """0 = strongest evidence (APPROVAL), larger = weaker."""
    try:
        return STAGE_PRIORITY.index(stage)
    except ValueError:
        return len(STAGE_PRIORITY)


def stage_color(stage: str):
    """Green gradient: darker = stronger clinical evidence."""
    if not stage:
        return "lightgray"
    try:
        rank = STAGE_PRIORITY.index(stage)
    except ValueError:
        return "lightgray"
    intensity = 1.0 - (rank / max(len(STAGE_PRIORITY) - 1, 1)) * 0.65
    return plt.cm.Greens(intensity)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-sets-dir", required=True)
    parser.add_argument("--studies", nargs="+", required=True)
    parser.add_argument("--drug-targets", required=True)
    parser.add_argument("--disease-ontology", required=True)
    parser.add_argument("--disease-id", nargs="+", required=True)
    parser.add_argument("--de-results", required=True)
    parser.add_argument(
        "--gene-set",
        required=True,
        choices=["regulators_H4", "targets_H4"],
        help="H4_genes is excluded (self-loop network is uninformative).",
    )
    parser.add_argument(
        "--condition", required=True, help="e.g. Rest, Stim8hr, Stim48hr"
    )
    parser.add_argument(
        "--title", default="", help="Optional disease label, e.g. 'Asthma (ASTCO)'"
    )
    parser.add_argument("--output", required=True, help="PNG output path")
    parser.add_argument(
        "--graphml",
        default=None,
        help="Optional GraphML export for Cytoscape / Gephi",
    )
    parser.add_argument(
        "--max-drugs-per-target",
        type=int,
        default=0,
        help=(
            "Cap drug nodes per hit gene (top K by best clinical stage). "
            "0 = no cap. Useful for very-drugged targets like NR3C1."
        ),
    )
    parser.add_argument(
        "--drug-stages",
        nargs="+",
        default=["APPROVAL", "PREAPPROVAL", "PHASE_4"],
        choices=APPROVED_STAGES,
        help=(
            "Clinical stages a drug must reach to appear as an outer-ring node. "
            "Default = Phase 4+ (APPROVAL/PREAPPROVAL/PHASE_4) for a clean "
            "'validated drug chain' figure. Hits with no drug at this threshold "
            "are dropped from the plot."
        ),
    )
    args = parser.parse_args()

    # 1. Disease subtree
    ontology = (
        pl.scan_parquet(args.disease_ontology).select(["id", "descendants"]).collect()
    )
    subtree = expand_subtree(args.disease_id, ontology)

    # 2. Approved drug evidence for the focal disease (incl. drug identity).
    drug_evidence = (
        pl.scan_parquet(args.drug_targets)
        .select(["targetId", "diseaseId", "clinicalStage", "drugId", "drugFromSource"])
        .filter(pl.col("diseaseId").is_in(subtree))
        .filter(pl.col("clinicalStage").is_in(APPROVED_STAGES))
        .collect()
    )
    approved = set(drug_evidence.get_column("targetId").unique().to_list())

    # 3. Gene set + coloc set
    gene_sets_dir = Path(args.gene_sets_dir)
    chosen = load_gene_set(gene_sets_dir, args.studies, args.gene_set, args.condition)
    h4_genes = load_gene_set(gene_sets_dir, args.studies, "H4_genes", args.condition)
    hits = chosen & approved
    if not hits:
        raise SystemExit("No hits — nothing to plot.")
    print(
        f"|chosen set|={len(chosen)}  |approved|={len(approved)}  "
        f"|hits|={len(hits)}  |coloc|={len(h4_genes)}"
    )

    # 4. Per-hit max clinical stage
    hit_max_stage: dict[str, str] = {}
    for ensg in hits:
        stages = (
            drug_evidence.filter(pl.col("targetId") == ensg)
            .get_column("clinicalStage")
            .unique()
            .to_list()
        )
        hit_max_stage[ensg] = pick_max_stage(stages)

    # 4b. Aggregate drugs per hit by drugId (deduplicates formulations/casing).
    # drug_targets_pairs[drug_id] = {"label", "best_stage", "targets": set(ensg)}
    drug_targets_pairs: dict[str, dict] = {}
    # For per-target capping: target_ensg -> list of (drug_id, stage)
    per_target_pairs: dict[str, list[tuple[str, str]]] = {ensg: [] for ensg in hits}

    drug_stage_filter = set(args.drug_stages)
    hit_drug_rows = drug_evidence.filter(pl.col("targetId").is_in(hits))
    for r in hit_drug_rows.iter_rows(named=True):
        target = r["targetId"]
        drug_id = r["drugId"]
        drug_name = r["drugFromSource"]
        stage = r["clinicalStage"]
        if drug_id is None or stage is None:
            continue
        if stage not in drug_stage_filter:
            continue
        per_target_pairs[target].append((drug_id, stage))
        label = (drug_name or drug_id).strip().lower()
        if drug_id not in drug_targets_pairs:
            drug_targets_pairs[drug_id] = {
                "label": label,
                "best_stage": stage,
                "targets": {target},
            }
        else:
            info = drug_targets_pairs[drug_id]
            info["targets"].add(target)
            if stage_rank(stage) < stage_rank(info["best_stage"]):
                info["best_stage"] = stage

    # Apply per-target cap (only the K best-stage drugs survive per target).
    if args.max_drugs_per_target > 0:
        keep_per_target: dict[str, set[str]] = {}
        for target, pairs in per_target_pairs.items():
            best_per_drug: dict[str, str] = {}
            for did, s in pairs:
                if did not in best_per_drug or stage_rank(s) < stage_rank(
                    best_per_drug[did]
                ):
                    best_per_drug[did] = s
            ranked = sorted(best_per_drug.items(), key=lambda kv: stage_rank(kv[1]))
            keep_per_target[target] = {
                d for d, _ in ranked[: args.max_drugs_per_target]
            }
        filtered: dict[str, dict] = {}
        for drug_id, info in drug_targets_pairs.items():
            kept_targets = {
                t for t in info["targets"] if drug_id in keep_per_target.get(t, set())
            }
            if kept_targets:
                filtered[drug_id] = {**info, "targets": kept_targets}
        drug_targets_pairs = filtered
    print(f"|drug nodes|={len(drug_targets_pairs)} (after dedup by drugId)")

    # 5. DE edges for this condition, gated on H4 membership
    de = pl.read_csv(args.de_results).filter(
        pl.col("culture_condition") == args.condition
    )
    if args.gene_set == "targets_H4":
        edges_df = de.filter(
            pl.col("trans_ensg").is_in(hits) & pl.col("cis_ensg").is_in(h4_genes)
        ).select(["cis_gene", "trans_gene", "trans_ensg", "logFC", "padj"])
        coloc_sym_col, hit_sym_col, hit_ensg_col = (
            "cis_gene",
            "trans_gene",
            "trans_ensg",
        )
    else:  # regulators_H4
        edges_df = de.filter(
            pl.col("cis_ensg").is_in(hits) & pl.col("trans_ensg").is_in(h4_genes)
        ).select(["cis_gene", "trans_gene", "cis_ensg", "logFC", "padj"])
        coloc_sym_col, hit_sym_col, hit_ensg_col = "trans_gene", "cis_gene", "cis_ensg"

    if edges_df.height == 0:
        raise SystemExit(
            "No DE edges link hits to coloc genes — nothing to plot. "
            "Check that --condition matches the gene-set file."
        )

    # 6. Build the graph (collapse multi-edges to strongest |logFC|).
    G = nx.Graph()
    for r in edges_df.iter_rows(named=True):
        coloc = r[coloc_sym_col]
        hit = r[hit_sym_col]
        hit_ensg = r[hit_ensg_col]
        lfc = r["logFC"]
        if coloc is None or hit is None:
            continue
        G.add_node(coloc, layer="coloc")
        G.add_node(
            hit,
            layer="hit",
            ensg=hit_ensg,
            max_stage=hit_max_stage.get(hit_ensg, ""),
        )
        if G.has_edge(coloc, hit):
            if abs(lfc) > abs(G[coloc][hit]["logfc"]):
                G[coloc][hit]["logfc"] = lfc
        else:
            G.add_edge(coloc, hit, logfc=lfc)

    # Hit-ensg -> hit-symbol lookup, used to wire drug nodes to their targets.
    ensg_to_symbol: dict[str, str] = {
        d["ensg"]: n
        for n, d in G.nodes(data=True)
        if d.get("layer") == "hit" and "ensg" in d
    }
    # Add drug nodes (outer shell) and hit -> drug edges.
    for drug_id, info in drug_targets_pairs.items():
        drug_node = f"drug:{drug_id}"
        G.add_node(
            drug_node,
            layer="drug",
            label=info["label"],
            max_stage=info["best_stage"],
        )
        for target_ensg in info["targets"]:
            target_sym = ensg_to_symbol.get(target_ensg)
            if target_sym is None:
                continue
            G.add_edge(drug_node, target_sym, kind="hit_drug")

    coloc_nodes = [n for n, d in G.nodes(data=True) if d.get("layer") == "coloc"]
    hit_nodes = [n for n, d in G.nodes(data=True) if d.get("layer") == "hit"]
    drug_nodes = [n for n, d in G.nodes(data=True) if d.get("layer") == "drug"]
    print(
        f"|nodes|={G.number_of_nodes()} (coloc={len(coloc_nodes)}, "
        f"hit={len(hit_nodes)}, drug={len(drug_nodes)})  "
        f"|edges|={G.number_of_edges()}"
    )

    # 7. Layout — three concentric shells with explicit radii so labels at
    # larger font sizes don't collide. Place coloc in the innermost circle
    # (or at the origin if there is only one), hits in the middle ring, drugs
    # in the outer ring.
    import math

    def _ring(nodes: list[str], radius: float) -> dict:
        if not nodes:
            return {}
        if len(nodes) == 1 and radius < 0.05:
            return {nodes[0]: (0.0, 0.0)}
        if len(nodes) == 1:
            return {nodes[0]: (radius, 0.0)}
        return {
            node: (
                radius * math.cos(2 * math.pi * i / len(nodes)),
                radius * math.sin(2 * math.pi * i / len(nodes)),
            )
            for i, node in enumerate(sorted(nodes))
        }

    R_COLOC, R_HIT, R_DRUG = 0.35, 1.0, 1.55
    pos: dict = {}
    if len(coloc_nodes) <= 1:
        pos.update(_ring(coloc_nodes, 0.0))
    else:
        pos.update(_ring(coloc_nodes, R_COLOC))
    pos.update(_ring(hit_nodes, R_HIT))
    pos.update(_ring(drug_nodes, R_DRUG))

    # 8. Plot.
    fig, ax = plt.subplots(figsize=(18, 15))
    ax.set_xlim(-2.0, 2.0)
    ax.set_ylim(-2.0, 2.0)

    # Split edges by kind so we don't KeyError on hit-drug edges (no logfc).
    coloc_hit_edges = [
        (u, v) for u, v, d in G.edges(data=True) if d.get("kind") != "hit_drug"
    ]
    hit_drug_edges = [
        (u, v) for u, v, d in G.edges(data=True) if d.get("kind") == "hit_drug"
    ]
    pos_edges = [(u, v) for u, v in coloc_hit_edges if G[u][v].get("logfc", 0) > 0]
    neg_edges = [(u, v) for u, v in coloc_hit_edges if G[u][v].get("logfc", 0) < 0]
    pos_widths = [min(abs(G[u][v]["logfc"]) * 1.2, 4.0) for u, v in pos_edges]
    neg_widths = [min(abs(G[u][v]["logfc"]) * 1.2, 4.0) for u, v in neg_edges]

    # Background hit -> drug edges first.
    nx.draw_networkx_edges(
        G,
        pos,
        edgelist=hit_drug_edges,
        edge_color="#222222",
        width=1.8,
        alpha=0.9,
        ax=ax,
    )
    # coloc -> hit signed edges over the top.
    nx.draw_networkx_edges(
        G,
        pos,
        edgelist=pos_edges,
        edge_color="#d62728",
        width=pos_widths,
        alpha=0.6,
        ax=ax,
    )
    nx.draw_networkx_edges(
        G,
        pos,
        edgelist=neg_edges,
        edge_color="#1f77b4",
        width=neg_widths,
        alpha=0.6,
        ax=ax,
    )

    # Coloc (inner) nodes.
    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=coloc_nodes,
        node_color="#ffd700",
        edgecolors="black",
        linewidths=2.5,
        node_size=3500,
        node_shape="s",
        ax=ax,
    )
    # Hit (middle) nodes coloured by best clinical stage for the focal disease.
    hit_colors = [stage_color(G.nodes[n].get("max_stage", "")) for n in hit_nodes]
    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=hit_nodes,
        node_color=hit_colors,
        edgecolors="black",
        linewidths=1.2,
        node_size=1500,
        ax=ax,
    )
    # Drug (outer) nodes — diamonds, smaller.
    drug_colors = [stage_color(G.nodes[n].get("max_stage", "")) for n in drug_nodes]
    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=drug_nodes,
        node_color=drug_colors,
        edgecolors="#555555",
        linewidths=0.8,
        node_size=900,
        node_shape="D",
        ax=ax,
    )

    # Labels: gene symbols for inner/middle (bold), drug names for outer.
    # White-stroke path effects keep black text readable against any node
    # fill colour (the dark green APPROVAL fill is the worst offender).
    stroke = [
        path_effects.Stroke(linewidth=7.5, foreground="white"),
        path_effects.Normal(),
    ]
    drug_stroke = [
        path_effects.Stroke(linewidth=5.5, foreground="white"),
        path_effects.Normal(),
    ]
    coloc_label_texts = nx.draw_networkx_labels(
        G,
        pos,
        labels={n: n for n in coloc_nodes},
        font_size=30,
        font_weight="bold",
        ax=ax,
    )
    for txt in coloc_label_texts.values():
        txt.set_path_effects(stroke)
    hit_label_texts = nx.draw_networkx_labels(
        G,
        pos,
        labels={n: n for n in hit_nodes},
        font_size=22,
        font_weight="bold",
        ax=ax,
    )
    for txt in hit_label_texts.values():
        txt.set_path_effects(stroke)

    def _trim_drug(name: str, max_len: int = 16) -> str:
        # Drop trailing parenthetical (e.g. "(removab)", "(mt103)") then truncate.
        n = name.split(" (", 1)[0].strip()
        if len(n) > max_len:
            n = n[: max_len - 1] + "…"
        return n.upper()

    drug_labels = {n: _trim_drug(G.nodes[n].get("label", "")) for n in drug_nodes}
    drug_label_texts = nx.draw_networkx_labels(
        G,
        pos,
        labels=drug_labels,
        font_size=24,
        font_color="#222222",
        ax=ax,
    )
    for txt in drug_label_texts.values():
        txt.set_path_effects(drug_stroke)

    legend_handles = [
        Patch(facecolor="#ffd700", edgecolor="black", label="coloc gene"),
        Patch(
            facecolor=stage_color("APPROVAL"),
            edgecolor="black",
            label="drug target",
        ),
        Line2D(
            [0],
            [0],
            marker="D",
            color="w",
            markerfacecolor=stage_color("APPROVAL"),
            markeredgecolor="#000000",
            markersize=14,
            label="drug",
        ),
    ]
    ax.legend(
        handles=legend_handles,
        loc="upper right",
        bbox_to_anchor=(1.0, 1.0),
        frameon=True,
        framealpha=0.92,
        edgecolor="#cccccc",
        fontsize=22,
    )

    title_prefix = f"{args.title} | " if args.title else ""
    ax.set_title(
        f"{title_prefix}{args.gene_set} | {args.condition}\n"
        f"{len(coloc_nodes)} coloc gene(s) → {len(hit_nodes)} drug-targeted hit(s) "
        f"→ {len(drug_nodes)} drug(s)",
        fontsize=22,
    )
    ax.axis("off")
    plt.tight_layout()
    fig.savefig(args.output, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {args.output}")

    if args.graphml:
        nx.write_graphml(G, args.graphml)
        print(f"Wrote {args.graphml}")


if __name__ == "__main__":
    main()
