"""
Supplementary Code 2
Cross-architecture regime classification and linear dynamical nulls for
"Mechanistic null models bound evidence for consciousness in aneural systems".

This code adds two adversarial analyses without changing the locked prospective
Stentor protocol or its confirmatory analysis.
"""
from pathlib import Path
import importlib.util
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from numpy.random import default_rng
from scipy.stats import wilcoxon
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, accuracy_score

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("supp1", HERE / "Supplementary_Code_1.py")
supp1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(supp1)

SEED = 20261016


def bootstrap_mean_ci(x, nboot=2000, seed=SEED):
    x = np.asarray(x, float)
    rng = default_rng(seed)
    means = np.empty(nboot)
    for i in range(nboot):
        means[i] = rng.choice(x, size=len(x), replace=True).mean()
    return tuple(np.quantile(means, [0.025, 0.975]))


def make_topology_reservoir(seed, topology, N=60, rho=1.05):
    rng = default_rng(seed)
    if topology == "ER":
        G = nx.gnp_random_graph(N, 0.10, seed=seed, directed=True)
        A = nx.to_numpy_array(G, dtype=float)
    elif topology == "WS":
        G = nx.watts_strogatz_graph(N, 6, 0.20, seed=seed)
        A = nx.to_numpy_array(G, dtype=float)
    elif topology == "BA":
        G = nx.barabasi_albert_graph(N, 3, seed=seed)
        A = nx.to_numpy_array(G, dtype=float)
    else:
        raise ValueError(topology)
    mean_degree = max(A.sum() / N, 1.0)
    W = rng.normal(0, 1 / np.sqrt(mean_degree), (N, N)) * A
    W *= rho / (supp1.spectral_radius(W) + 1e-12)
    Win = rng.normal(0, 0.6, N)
    return W, Win


def topology_dataset(n_networks=60):
    rows = []
    for ti, topology in enumerate(["ER", "WS", "BA"]):
        for network in range(n_networks):
            seed = 1000 * (ti + 1) + network
            W, Win = make_topology_reservoir(seed, topology)
            for suppressed, gain in [(0, 1.0), (1, 0.45)]:
                dim, lz, comp, energy, divergence = supp1.perturb_metrics(W, Win, gain, network)
                rows.append((topology, network, suppressed, dim, lz, comp,
                             energy, divergence))
    return pd.DataFrame(rows, columns=[
        "topology", "network", "suppressed", "effective_dimension",
        "lz_complexity", "compressed_complexity", "energy",
        "divergence"
    ])


def cross_topology_classifier(df):
    # Energy is deliberately excluded so transfer is not a trivial amplitude detector.
    features = ["effective_dimension", "lz_complexity", "compressed_complexity",
                "divergence"]
    train = df[df.topology == "ER"]
    clf = make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=5000,
                                                              random_state=SEED))
    clf.fit(train[features], train["suppressed"])
    rows = []
    pred_rows = []
    for topology in ["ER", "WS", "BA"]:
        d = df[df.topology == topology].copy()
        prob = clf.predict_proba(d[features])[:, 1]
        pred = (prob >= 0.5).astype(int)
        auc = roc_auc_score(d["suppressed"], prob)
        acc = accuracy_score(d["suppressed"], pred)
        temp = d[["network", "suppressed"]].copy()
        temp["prob_suppressed"] = prob
        piv = temp.pivot(index="network", columns="suppressed", values="prob_suppressed")
        pair_rank = float(np.mean(piv[1] > piv[0]))
        margin = (piv[1] - piv[0]).to_numpy()
        ci = bootstrap_mean_ci(margin, seed=SEED + ["ER", "WS", "BA"].index(topology))
        rows.append((topology, auc, acc, pair_rank, float(margin.mean()), ci[0], ci[1]))
        d["prob_suppressed"] = prob
        pred_rows.append(d)
    summary = pd.DataFrame(rows, columns=[
        "topology", "auc", "accuracy_at_0.5", "paired_ranking_fraction",
        "mean_within_network_probability_margin", "margin_ci_low", "margin_ci_high"
    ])
    return summary, pd.concat(pred_rows, ignore_index=True)


def run_linear(W, Win, inputs, gain=1.0, leak=0.35, x0=None):
    N = W.shape[0]
    x = np.zeros(N) if x0 is None else x0.copy()
    states = np.zeros((len(inputs), N))
    for t, u in enumerate(inputs):
        # No nonlinearity and no adaptation state.
        x = (1 - leak) * x + leak * (gain * (W @ x) + Win * u)
        states[t] = x
    return states


def linear_metrics(W, Win, gain=1.0):
    u = np.zeros(160)
    u[20] = 1.0
    S = run_linear(W, Win, u, gain=gain)
    post = S[20:140]
    return (
        supp1.participation_ratio(post),
        supp1.normalized_lz(post),
        supp1.compression_complexity(post),
        float(np.mean(np.sum(post ** 2, axis=1))),
    )


def linear_null_ensemble(n_networks=100):
    rows = []
    for network in range(n_networks):
        # rho < 1 ensures the baseline linear dynamics are stable.
        W, Win = supp1.make_reservoir(network, N=60, p=0.10, rho=0.95)
        b = linear_metrics(W, Win, 1.0)
        s = linear_metrics(W, Win, 0.45)
        rows.append((network, *b, *s))
    raw = pd.DataFrame(rows, columns=[
        "network", "dimension_baseline", "lz_baseline", "compression_baseline", "energy_baseline",
        "dimension_suppressed", "lz_suppressed", "compression_suppressed", "energy_suppressed"
    ])
    summaries = []
    for metric in ["dimension", "lz", "compression", "energy"]:
        b = raw[f"{metric}_baseline"].to_numpy()
        s = raw[f"{metric}_suppressed"].to_numpy()
        delta = s - b
        ci = bootstrap_mean_ci(delta, seed=SEED + len(summaries) + 10)
        dz = delta.mean() / delta.std(ddof=1)
        w = wilcoxon(delta)
        summaries.append((metric, b.mean(), s.mean(), delta.mean(), ci[0], ci[1],
                          dz, float(w.pvalue), float(np.mean(delta < 0))))
    summary = pd.DataFrame(summaries, columns=[
        "metric", "baseline_mean", "suppressed_mean", "paired_difference_suppressed_minus_baseline",
        "difference_ci_low", "difference_ci_high", "paired_dz", "wilcoxon_p",
        "fraction_decreasing"
    ])
    return summary, raw


def effective_marker_count(k=6, rho=None):
    if rho is None:
        rho = np.linspace(0, 0.95, 96)
    rho = np.asarray(rho)
    return k / (1 + (k - 1) * rho)


def make_figure(linear_summary, predictions, outpath):
    fig = plt.figure(figsize=(11, 8.3))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.30)

    # A: linear null paired means
    ax = fig.add_subplot(gs[0, 0])
    metrics = ["dimension", "lz", "compression"]
    x = np.arange(len(metrics))
    b = [linear_summary.loc[linear_summary.metric == m, "baseline_mean"].iloc[0] for m in metrics]
    s = [linear_summary.loc[linear_summary.metric == m, "suppressed_mean"].iloc[0] for m in metrics]
    width = 0.35
    ax.bar(x - width/2, b, width, label="Baseline")
    ax.bar(x + width/2, s, width, label="Suppressed")
    ax.set_xticks(x, ["Effective\ndimension", "LZ\ncomplexity", "Compressed\ncomplexity"])
    ax.set_ylabel("Mean metric value")
    ax.set_title("a  Linear systems reproduce complexity collapse", loc="left")
    ax.legend(frameon=False)

    # B: classifier predictions by topology
    ax = fig.add_subplot(gs[0, 1])
    positions = []
    labels = []
    data = []
    pos = 1
    for topo in ["ER", "WS", "BA"]:
        for suppressed, state in [(0, "Base"), (1, "Supp")]:
            vals = predictions[(predictions.topology == topo) &
                               (predictions.suppressed == suppressed)].prob_suppressed.to_numpy()
            data.append(vals); positions.append(pos); labels.append(f"{topo}\n{state}")
            pos += 1
        pos += 0.5
    ax.boxplot(data, positions=positions, widths=0.7, showfliers=False)
    ax.axhline(0.5, linewidth=1, linestyle="--")
    ax.set_xticks(positions, labels, fontsize=8)
    ax.set_ylim(-0.03, 1.03)
    ax.set_ylabel("Predicted suppressed-state probability")
    ax.set_title("b  ER-trained marker bundle transfers across topology", loc="left")

    # C: within-network probability margin
    ax = fig.add_subplot(gs[1, 0])
    margins = []
    for topo in ["ER", "WS", "BA"]:
        d = predictions[predictions.topology == topo]
        piv = d.pivot(index="network", columns="suppressed", values="prob_suppressed")
        margins.append((piv[1] - piv[0]).to_numpy())
    ax.boxplot(margins, labels=["ER", "WS", "BA"], showfliers=False)
    ax.axhline(0, linewidth=1, linestyle="--")
    ax.set_ylabel("Within-network probability margin\n(suppressed - baseline)")
    ax.set_title("c  Every paired network is ranked correctly", loc="left")

    # D: analytic effective marker number
    ax = fig.add_subplot(gs[1, 1])
    rho = np.linspace(0, 0.95, 96)
    ax.plot(rho, effective_marker_count(6, rho))
    ax.axhline(1, linewidth=1, linestyle="--")
    ax.set_xlabel("Pairwise marker correlation, ρ")
    ax.set_ylabel("Effective independent-marker count")
    ax.set_title("d  Correlation collapses nominal marker multiplicity", loc="left")
    ax.set_xlim(0, 0.95)
    ax.set_ylim(0.8, 6.2)

    fig.suptitle("Extended Data Figure 3 | Consciousness-linked metrics can identify a generic dynamical regime",
                 y=0.98, fontsize=12)
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    topo = topology_dataset(30)
    clf_summary, predictions = cross_topology_classifier(topo)
    linear_summary, linear_raw = linear_null_ensemble(100)

    clf_summary.to_csv(HERE / "ExtendedData_Table3_cross_topology_classifier.csv", index=False)
    predictions.to_csv(HERE / "ExtendedData_Table3b_classifier_predictions.csv", index=False)
    linear_summary.to_csv(HERE / "ExtendedData_Table4_linear_nulls.csv", index=False)
    linear_raw.to_csv(HERE / "ExtendedData_Table4b_linear_null_raw.csv", index=False)
    make_figure(linear_summary, predictions, HERE / "ExtendedData_Figure3_regime_classifier.png")

    print("Cross-topology classifier (energy excluded):")
    print(clf_summary.to_string(index=False))
    print("\nLinear null ensemble:")
    print(linear_summary.to_string(index=False))
    print("\nEffective marker count for k=6:")
    for r in [0, 0.2, 0.5, 0.8, 0.95]:
        print(r, effective_marker_count(6, [r])[0])


if __name__ == "__main__":
    main()
