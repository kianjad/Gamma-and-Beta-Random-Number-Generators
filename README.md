# Gamma and Beta Random Number Generators
**Jayden Gould, Harrison Funk, Elijah DiFuria, Kian Jadbabaei**  
*STAT 194CS — Monte Carlo Methods | Winter 2026*

---

## Overview

This project implements a suite of acceptance-rejection sampling algorithms for generating Gamma and Beta random variables in R, based on the 1974 technical report by U. Dieter and J.H. Ahrens:

> *"Acceptance-Rejection Techniques for Sampling from the Gamma and Beta Distributions"*  
> Stanford University, Department of Statistics — Technical Report No. 83

The goal was to go beyond theoretical understanding and build working, verified random number generators that mirror the architecture of foundational statistical software like R's built-in `rgamma()`.

---

## Algorithms Implemented

The core insight from Dieter & Ahrens (1974) is that no single acceptance-rejection envelope can cover the entire domain of the Gamma shape parameter `a ∈ (0, ∞)`. Our solution is a **master routing function** that dynamically dispatches to the appropriate algorithm based on `a`.

| Algorithm | Parameter Range | Envelope Type |
|-----------|----------------|---------------|
| **GS** | `0 < a ≤ 1` | Piecewise power/exponential (Lemma 4) |
| **GC** | `1 < a ≤ 2.53` | Cauchy distribution (Lemma 1) |
| **GO** | `a > 2.53` | Normal-Exponential split with polynomial squeeze bounds (Lemma 3) |
| **Beta (Gamma-Ratio)** | Any `a, b > 0` | Built on top of the Gamma generator |

The threshold `a = 2.53` is not an arbitrary tuning choice — it is the exact value where the cubic discriminant in Lemma 3 becomes positive, below which GO's polynomial squeeze bounds are mathematically invalid.

---

## Key Features

- **Master wrapper function** `generate_gamma_complete(n, a)` routes automatically to GS, GC, or GO
- **Algorithm GO** implements the "Quick Acceptance" squeeze bounds from Lemma 3, bypassing expensive log/exp evaluations in the vast majority of iterations
- **Beta generator** `generate_beta_from_gamma(n, a, b)` uses the Gamma-ratio method: if `x ~ Gamma(a)` and `y ~ Gamma(b)`, then `x/(x+y) ~ Beta(a,b)`
- All acceptance tests evaluated on **log scale** to prevent floating-point overflow
- **Empirical verification** via Kolmogorov-Smirnov tests and histogram overlays for all four distributions

---

## Verification Results

All modules were tested with `n = 10,000` samples and evaluated against theoretical CDFs using `ks.test()`.

| Distribution | Algorithm | K-S p-value |
|---|---|---|
| Gamma(a=0.5) | GS | 0.615 |
| Gamma(a=2.0) | GC | 0.629 |
| Gamma(a=10.0) | GO | 0.568 |
| Beta(2.5, 5.0) | Gamma-Ratio | 0.792 |

All p-values are well above the 0.05 threshold — we fail to reject the null hypothesis that the empirical samples match the theoretical distributions.

---

## File Structure

```
├── Gamma_and_Beta_Random_Number_Generators.Rmd   # Full report with code and analysis
└── Gamma-and-Beta-Random-Number-Generators.pdf   # Rendered PDF output
```

---

## How to Run

Open `Gamma_and_Beta_Random_Number_Generators.Rmd` in RStudio and knit to PDF or HTML. No external packages are required — only base R.

```r
# Generate 10,000 Gamma samples using the master wrapper
samples <- generate_gamma_complete(n = 10000, a = 5.0)

# Generate 10,000 Beta samples via the Gamma-ratio method
beta_samples <- generate_beta_from_gamma(n = 10000, a = 2.5, b = 5.0)
```

---

## What We Learned

Implementing these algorithms taught us that software functions like `rgamma()` are not single equations — they are sophisticated routing systems backed by algebraic proofs. The `a = 2.53` boundary isn't a programmer's guess; it's the solution to a cubic discriminant. Understanding *why* the boundary exists was what made the implementation actually correct.

---

## Reference

Dieter, U. and Ahrens, J.H. (1974). *Acceptance-Rejection Techniques for Sampling from the Gamma and Beta Distributions.* Technical Report No. 83, Department of Statistics, Stanford University.
