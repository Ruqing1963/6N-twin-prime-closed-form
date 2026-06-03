# 6N Twin-Prime Closed-Form Survival (Part VII)

A closed form for the conditional twin-gap distribution on the 6N ± 1 skeleton —
resolving the open problem of Part VI.

**Background.** Part VI's bridge theorem reduced the gap preference to the
right-centre survival via an ω-independent constant: r(d|ω) = P(N+d twin|ω)/C(d).
It left one problem: a closed form for the **absolute** survival P(N+d twin|ω)
from per-prime data, since a naive independent product overstates it.

**The closed form (S₁₀, 23,988,173 twin centres, no fitted parameters).**

```
    P(N+d twin | N twin, ω) = K · ∏_{q>3} f_q(d, N)
```

with the **pure CRT survival factor** (dead(q) = {±6⁻¹ mod q}):

```
    q | N :  f_q = 1 if (d mod q) ∉ dead(q) else 0
    q ∤ N :  f_q = #{ r ∈ A_q : (r+d) ∉ dead(q) } / #A_q
             A_q = admissible nonzero residues { r ≠ 0 : r ∉ dead(q) }
```

and **K** an ω-independent tail constant (primes beyond the working set {5..47}),
fixed once from the ω-merged ratio.

**Result.** Error ≤ 0.5% for ω ≤ 5, and ≤ 2.5% at the sparsest stratum ω=6, for
both 6ΔN=42 and 210. The two tail constants (K=0.1033 for 42, K=0.1052 for 210)
are nearly equal — consistent with a single **gap-independent** K. No correction
term is needed for the Part V cross-prime hedge: it is captured automatically by
selecting the q|N or q∤N branch of f_q per prime according to N's factorisation.

**End to end.** Composed with the Part VI bridge constant:

```
    r(d|ω) = (K / C(d)) · ⟨ ∏_{q>3} f_q(d, N) ⟩_ω
```

So r(d|ω) — the rise of 42 to 1.55, the collapse of 210 to 0.41 — is now a
*computed* quantity from the centre's factorisation, not a described one.

> **Precision boundary (reported, not smoothed).** The only departures above 0.5%
> are at ω=6 (−2.1% for 42, +2.5% for 210). This is where the uniform-residue
> approximation in the q∤N branch of f_q is least accurate: real twin-centre
> residues carry the Part I enrichment bias, not exact uniformity over A_q. We
> identify the ω=6 residual as this enrichment correction and leave quantifying
> it as the natural refinement — it is not absorbed into a fitted term.
>
> **Scope.** Experimental / computational number theory. No claim about the
> infinitude of twin primes or any prime k-tuple conjecture.

Part I: doi:10.5281/zenodo.20470367 · II: doi:10.5281/zenodo.20477664 ·
III: doi:10.5281/zenodo.20498668 · IV: doi:10.5281/zenodo.20500465 ·
V: doi:10.5281/zenodo.20510700 · VI: doi:10.5281/zenodo.20517990

---

## Layout

```
.
├── README.md
├── LICENSE                 (MIT)
├── CITATION.cff
├── data/
│   └── closed_S10_data.csv  gap, omega, measured_P, closed_form, err_pct  (S10)
├── code/
│   ├── closed_form.py       evaluates P = K·∏ f_q (pure CRT) vs measured survival
│   │                        per omega; prints K and err%; emits closed_S{K}_data.csv
│   └── make_closed_fig.py   builds the 3-panel closed-form-vs-measured figure
├── figures/                fig_paper7_closed.{pdf,png}
└── paper/                  Chen_6N_Paper7.{tex,pdf} + figure
```

## Reproducing

Requirements: Python 3.8+, `numpy`, `matplotlib`.

```bash
pip install numpy matplotlib

# 1. Evaluate the closed form. Default S10 (~11 min scan). Prints K and per-omega
#    error; writes closed_S{K}_data.csv.
python code/closed_form.py            # S10
MAXK=9 python code/closed_form.py     # S9 (faster, for validation)

# 2. Build the figure (reads ../data/closed_S10_data.csv).
cd code && python make_closed_fig.py
```

### Conventions (same as Parts II–VI)

- Twin centre: N with 6N−1, 6N+1 both prime. Centre-step d; physical gap 6d.
- dead(q) = {±6⁻¹ mod q}: a centre is q-safe iff its residue ∉ dead(q).
- Right centre M = N+d; q|N forces M ≡ d (mod q). The factor f_q is the pure CRT
  conditional safety of M given the left centre's divisibility by q.
- K is fixed once (ω-merged), then held constant across all ω.
- Working primes q ∈ {5,7,…,47}; K carries the tail q > 47.
- Engine: complete segmented-sieve factorisation + deterministic interval-sieve
  primality; S₁₀ twin count 23,988,173 matches Part I. "N+d is a twin centre"
  tested by binary search on the sorted twin array.

## License

MIT — see `LICENSE`.
