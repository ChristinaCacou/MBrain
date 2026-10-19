"""
Supplementary Code 1
Aneural adaptive dynamics recapitulate signatures often linked to consciousness

Reproducible simulations for the computational analyses in the manuscript.
Python 3.11+; requires numpy, pandas, scipy, networkx, matplotlib.
"""
import numpy as np
import pandas as pd
import networkx as nx
from numpy.random import default_rng
from scipy.stats import wilcoxon
import zlib


def spectral_radius(W):
    return float(np.max(np.abs(np.linalg.eigvals(W))))


def make_reservoir(seed, N=60, p=0.10, rho=1.05):
    rng = default_rng(seed)
    mask = rng.random((N, N)) < p
    W = rng.normal(0, 1/np.sqrt(p*N), (N, N)) * mask
    W *= rho / (spectral_radius(W) + 1e-12)
    Win = rng.normal(0, 0.6, N)
    return W, Win


def run_reservoir(W, Win, inputs, gain=1.0, leak=0.35,
                  adapt_beta=0.15, adapt_tau=0.02, x0=None, a0=None):
    N = W.shape[0]
    x = np.zeros(N) if x0 is None else x0.copy()
    a = np.zeros(N) if a0 is None else a0.copy()
    states = np.zeros((len(inputs), N))
    for t, u in enumerate(inputs):
        h = gain * (W @ x) + Win * u - adapt_beta * a
        x = (1-leak) * x + leak * np.tanh(h)
        a = (1-adapt_tau) * a + adapt_tau * x
        states[t] = x
    return states


def participation_ratio(states):
    X = states - states.mean(0, keepdims=True)
    C = np.cov(X, rowvar=False)
    ev = np.clip(np.linalg.eigvalsh(C), 0, None)
    return float(ev.sum()**2 / (np.sum(ev**2) + 1e-12))


def lz76_count(bits):
    s = ''.join('1' if b else '0' for b in bits)
    n = len(s)
    if n == 0:
        return 0
    i = 0; c = 1; l = 1; k = 1; k_max = 1
    while True:
        if i+k >= n:
            c += 1; break
        if s[i+k-1] == s[l+k-1]:
            k += 1
            if l+k > n:
                c += 1; break
        else:
            if k > k_max:
                k_max = k
            i += 1
            if i == l:
                c += 1
                l += k_max
                if l+1 > n:
                    break
                i = 0; k = 1; k_max = 1
            else:
                k = 1
    return c


def normalized_lz(states):
    threshold = np.median(states, axis=0, keepdims=True)
    bits = (states > threshold).astype(np.uint8).ravel()
    n = len(bits)
    return float(lz76_count(bits) * np.log2(max(n, 2)) / n)


def compression_complexity(states):
    threshold = np.median(states, axis=0, keepdims=True)
    packed = np.packbits((states > threshold).astype(np.uint8).ravel())
    return len(zlib.compress(packed.tobytes(), level=9)) / len(packed.tobytes())


def memory_capacity(W, Win, gain=1.0, seed=0, T=1400, wash=200, maxlag=15):
    rng = default_rng(seed + 990)
    u = rng.normal(size=T)
    X = run_reservoir(W, Win, u, gain=gain)
    start = wash + maxlag
    idx = np.arange(start, T)
    Xd = X[idx]
    cut = int(len(idx) * 0.65)
    Xtr, Xte = Xd[:cut], Xd[cut:]
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-6
    Xtr = (Xtr - mu) / sd
    Xte = (Xte - mu) / sd
    Xtr = np.c_[np.ones(len(Xtr)), Xtr]
    Xte = np.c_[np.ones(len(Xte)), Xte]
    lam = 1e-3
    A = Xtr.T @ Xtr + lam * np.eye(Xtr.shape[1])
    capacity = 0.0
    for lag in range(1, maxlag+1):
        y = u[idx-lag]
        ytr, yte = y[:cut], y[cut:]
        w = np.linalg.solve(A, Xtr.T @ ytr)
        pred = Xte @ w
        r2 = 1 - np.sum((yte-pred)**2) / (np.sum((yte-yte.mean())**2) + 1e-12)
        capacity += max(0.0, r2)
    return float(capacity)


def perturb_metrics(W, Win, gain=1.0, seed=0):
    u = np.zeros(160); u[20] = 1.0
    S = run_reservoir(W, Win, u, gain=gain)
    post = S[20:140]
    pr = participation_ratio(post)
    lz = normalized_lz(post)
    comp = compression_complexity(post)
    energy = float(np.mean(np.sum(post**2, axis=1)))
    x0 = S[20].copy()
    rng = default_rng(seed + 1234)
    delta = rng.normal(size=len(x0)); delta = 1e-7 * delta / np.linalg.norm(delta)
    Z = np.zeros(45)
    S1 = run_reservoir(W, Win, Z, gain=gain, x0=x0)
    S2 = run_reservoir(W, Win, Z, gain=gain, x0=x0+delta)
    dist = np.linalg.norm(S2-S1, axis=1) + 1e-15
    divergence = float(np.polyfit(np.arange(1,16), np.log(dist[:15]), 1)[0])
    return pr, lz, comp, energy, divergence


def simulate_habituation(seed=12, ncells=1200, n_train=30,
                         p_hab_mean=0.085, p_rec_between=0.002,
                         p_rec_rest=0.045, rest_steps=20,
                         p_resp_R=0.96, p_resp_H=0.08):
    rng = default_rng(seed)
    concentration = 30
    p_hab = rng.beta(p_hab_mean*concentration,
                     (1-p_hab_mean)*concentration, size=ncells)
    state = np.zeros(ncells, dtype=np.int8)
    curve = []
    responses = np.zeros((ncells, n_train+4), dtype=np.int8)
    for t in range(n_train):
        prob = np.where(state==0, p_resp_R, p_resp_H)
        r = rng.random(ncells) < prob
        responses[:,t] = r
        curve.append(r.mean())
        state[(state==0) & (rng.random(ncells) < p_hab)] = 1
        state[(state==1) & (rng.random(ncells) < p_rec_between)] = 0
    for _ in range(rest_steps):
        state[(state==1) & (rng.random(ncells) < p_rec_rest)] = 0
    post_rest = []
    for j in range(3):
        prob = np.where(state==0, p_resp_R, p_resp_H)
        r = rng.random(ncells) < prob
        post_rest.append(r.mean())
        responses[:,n_train+j] = r
        state[(state==0) & (rng.random(ncells) < p_hab)] = 1
    r = rng.random(ncells) < 0.97
    post_rest.append(r.mean())
    responses[:,n_train+3] = r
    return np.asarray(curve), np.asarray(post_rest), responses


def virtual_trap_closure(gain, n=20000, seed=22, touch_amp=0.62,
                         decay=0.86, ntouches=2):
    rng = default_rng(seed)
    threshold = np.clip(rng.normal(1.0, 0.09, n), 0.75, 1.25)
    v = np.zeros(n)
    for _ in range(ntouches):
        v = decay*v + gain*touch_amp
    return float(np.mean(v >= threshold))


def random_grid_maze(seed, n=8, remove_frac=0.15):
    rng = default_rng(seed)
    G = nx.grid_2d_graph(n,n)
    edges = list(G.edges()); rng.shuffle(edges)
    target = int(remove_frac * len(edges)); removed = 0
    for e in edges:
        if removed >= target: break
        G.remove_edge(*e)
        if nx.has_path(G,(0,0),(n-1,n-1)):
            removed += 1
        else:
            G.add_edge(*e)
    for u,v in G.edges():
        G[u][v]['length'] = float(rng.uniform(0.8,1.2))
    return G


def physarum_adapt(G, source=(0,0), sink=(7,7), seed=0,
                    steps=160, dt=0.12, gamma=1.0):
    rng = default_rng(seed+10000)
    nodes = list(G.nodes()); idx = {node:i for i,node in enumerate(nodes)}
    edges = list(G.edges()); n = len(nodes)
    L = np.array([G[u][v]['length'] for u,v in edges], float)
    D = np.exp(rng.normal(0, 0.08, len(edges)))
    s, t = idx[source], idx[sink]
    keep = np.array([i for i in range(n) if i != t])
    b = np.zeros(n); b[s] = 1; b[t] = -1
    for _ in range(steps):
        A = np.zeros((n,n)); c = D/L
        for k,(u,v) in enumerate(edges):
            i,j = idx[u],idx[v]; ck = c[k]
            A[i,i] += ck; A[j,j] += ck; A[i,j] -= ck; A[j,i] -= ck
        p = np.zeros(n)
        p[keep] = np.linalg.lstsq(A[np.ix_(keep,keep)], b[keep], rcond=1e-10)[0]
        Q = np.array([D[k]/L[k]*(p[idx[u]]-p[idx[v]]) for k,(u,v) in enumerate(edges)])
        D = np.maximum((1-dt)*D + dt*np.abs(Q)**gamma, 1e-6)
    H = nx.Graph()
    for k,(u,v) in enumerate(edges):
        H.add_edge(u,v,cost=L[k]/(D[k]+1e-12))
    path = nx.shortest_path(H, source, sink, weight='cost')
    selected = sum(G[path[i]][path[i+1]]['length'] for i in range(len(path)-1))
    optimum = nx.shortest_path_length(G, source, sink, weight='length')
    return float(selected/optimum)


def main():
    # Recurrent-network ensemble
    rows = []
    for seed in range(100):
        W, Win = make_reservoir(seed)
        for condition, gain in [('baseline',1.0), ('suppressed',0.45)]:
            pr,lz,comp,energy,div = perturb_metrics(W,Win,gain,seed)
            mc = memory_capacity(W,Win,gain,seed)
            rows.append((seed,condition,pr,lz,comp,energy,div,mc))
    df = pd.DataFrame(rows, columns=['seed','condition','effective_dimension','lz_complexity',
                                     'compressed_complexity','energy','divergence','memory_capacity'])
    print(df.groupby('condition').mean(numeric_only=True))
    for metric in ['effective_dimension','lz_complexity','compressed_complexity','energy','divergence','memory_capacity']:
        x = df.pivot(index='seed', columns='condition', values=metric)
        print(metric, wilcoxon(x['baseline'],x['suppressed']).pvalue)

    # Habituation
    curve, post, _ = simulate_habituation()
    print('Habituation: trial 1, trial 30, post-rest, strong stimulus:',
          curve[0], curve[-1], post[0], post[-1])

    # Excitable threshold / anaesthetic-like suppression
    print('Trap closure baseline/suppressed:', virtual_trap_closure(1.0), virtual_trap_closure(0.35))

    # Physarum-like adaptive routing
    ratios = []
    for seed in range(120):
        ratios.append(physarum_adapt(random_grid_maze(seed), seed=seed))
    ratios = np.asarray(ratios)
    print('Maze ratio mean/median/max:', ratios.mean(), np.median(ratios), ratios.max())


# ---- Added robustness analyses for strengthened manuscript ----
from scipy.stats import qmc


def perturb_metrics_general(W, Win, gain=1.0, leak=0.35,
                            adapt_beta=0.15, adapt_tau=0.02,
                            post_len=100):
    """Fast perturbational metrics used for the space-filling robustness sweep."""
    u = np.zeros(20 + post_len + 10)
    u[20] = 1.0
    S = run_reservoir(W, Win, u, gain=gain, leak=leak,
                      adapt_beta=adapt_beta, adapt_tau=adapt_tau)
    post = S[20:20+post_len]
    pr = participation_ratio(post)
    comp = compression_complexity(post)
    energy = float(np.mean(np.sum(post**2, axis=1)))
    return pr, comp, energy


def sobol_robustness(n_power=9, seed=20261016):
    """Space-filling sweep over seven dynamical parameters (2**n_power systems)."""
    Z = qmc.Sobol(d=7, scramble=True, seed=seed).random_base2(m=n_power)
    rows = []
    for i, z in enumerate(Z):
        p = 0.05 + 0.20*z[0]
        rho = 0.75 + 0.60*z[1]
        leak = 0.15 + 0.50*z[2]
        beta = 0.35*z[3]
        tau = np.exp(np.log(0.005) + z[4]*(np.log(0.08)-np.log(0.005)))
        g0 = 0.80 + 0.45*z[5]
        suppression_ratio = 0.30 + 0.40*z[6]
        g1 = g0*suppression_ratio
        W, Win = make_reservoir(200000+i, N=40, p=p, rho=rho)
        b = perturb_metrics_general(W, Win, g0, leak, beta, tau)
        s = perturb_metrics_general(W, Win, g1, leak, beta, tau)
        rows.append((i, p, rho, leak, beta, tau, g0, suppression_ratio,
                     b[0], s[0], b[1], s[1], b[2], s[2]))
    return pd.DataFrame(rows, columns=[
        'sample','edge_density','spectral_radius','leak','adapt_beta','adapt_tau',
        'baseline_gain','suppression_ratio','dimension_baseline','dimension_suppressed',
        'compression_baseline','compression_suppressed','energy_baseline','energy_suppressed'])


def output_scaling_control(scale=0.45, n_networks=100):
    """Post-hoc amplitude attenuation control; complexity metrics are scale-invariant."""
    rows = []
    for seed in range(n_networks):
        W, Win = make_reservoir(seed)
        u = np.zeros(160); u[20] = 1.0
        S = run_reservoir(W, Win, u, gain=1.0)
        post = S[20:140]
        scaled = scale*post
        rows.append((seed,
                     participation_ratio(post), participation_ratio(scaled),
                     compression_complexity(post), compression_complexity(scaled),
                     float(np.mean(np.sum(post**2, axis=1))),
                     float(np.mean(np.sum(scaled**2, axis=1)))))
    return pd.DataFrame(rows, columns=['seed','dimension','dimension_scaled',
                                        'compression','compression_scaled',
                                        'energy','energy_scaled'])


def robustness_summary():
    """Return the principal robustness and amplitude-control summary used in the manuscript."""
    rob = sobol_robustness(n_power=9, seed=20261016)
    dim = float(np.mean(rob.dimension_suppressed < rob.dimension_baseline))
    comp = float(np.mean(rob.compression_suppressed < rob.compression_baseline))
    both = float(np.mean((rob.dimension_suppressed < rob.dimension_baseline) &
                         (rob.compression_suppressed < rob.compression_baseline)))
    ctrl = output_scaling_control(scale=0.45, n_networks=100)
    return {
        'sobol_n': len(rob),
        'dimension_decrease_fraction': dim,
        'compression_decrease_fraction': comp,
        'both_decrease_fraction': both,
        'amplitude_control_dimension_max_abs_change': float(np.max(np.abs(ctrl.dimension-ctrl.dimension_scaled))),
        'amplitude_control_compression_max_abs_change': float(np.max(np.abs(ctrl.compression-ctrl.compression_scaled))),
        'amplitude_control_mean_energy_ratio': float(np.mean(ctrl.energy_scaled/ctrl.energy)),
    }


if __name__ == '__main__':
    main()
    print('Robustness summary:', robustness_summary())
