# J/ψ Spin Asymmetry Analysis in Ultra-Peripheral Collisions

**STAR Experiment** | p+p √s = 510 GeV (2017) and p+Au √s_NN = 200 GeV (2015)

---

## Overview

This repository contains the full analysis chain for measuring the J/ψ single-transverse-spin asymmetry (A_N) in ultra-peripheral collisions (UPC) at STAR. The analysis targets exclusive photo-production of J/ψ via γp → J/ψ → e⁺e⁻, where the proton serves as the polarized target. The measured asymmetry provides sensitivity to the gluon GPD E_g, which is related to the orbital angular momentum carried by gluons inside the proton.

The raw asymmetry is constructed as:

```
ε = ( √(N_L+ · N_R−) − √(N_L− · N_R+) ) / ( √(N_L+ · N_R−) + √(N_L− · N_R+) )
```

where N_{L/R ±} are J/ψ candidate counts in the left/right halves of the detector for spin-up(+)/spin-down(−) beam bunches.

---

## Datasets

| System | Year | √s_NN | Production Tag | Stream |
|--------|------|-------|----------------|--------|
| p+p    | 2017 | 510 GeV  | Production_pp500_2017 / P20ic | st_rp |
| p+Au   | 2015 | 200 GeV  | Production_pAu200_2015 / P16id | st_rp |

### Triggers

The primary trigger is **HTTPJPsi**: requires back-to-back high towers in BEMC (or EEMC) phi sextants, with BBC veto to suppress hadronic backgrounds.

- **p+p** (Trigger IDs: 570209, 570219, 570229): BEMC active, 2–6 TOF hits required
- **p+Au** (Trigger IDs: 500730, 500750, 500710): ZDC veto to reject nuclear breakup

---

## Analysis Chain

The analysis runs in three sequential stages:

```
MuDST files
    │
    ▼
[1] uDstSkim / uDstSkim_pA     ← skim MuDST → reduced sDST
    │
    ▼
[2] sDstAna / sDstAna_pA       ← fill histograms & trees from sDST
    │
    ▼
[3] hadd output files
    │
    ▼
[4] __draw_mee_fast.C           ← produce final plots & asymmetry results
```

---

## Step 1: MuDST Skimming

Skims the full MuDST files, applying trigger selection and basic event/track quality cuts to produce a reduced output suitable for the analysis.

**For p+p:**
```bash
star-submit uDstSkim.xml
```

**For p+Au:**
```bash
star-submit uDstSkim_pA.xml
```

### Event Quality Cuts Applied at Skim Stage
- |Z_vtx| < 100 cm (vertex reconstruction quality)
- Back-to-back trigger sextant requirement
- Number of BEMC clusters < 4

### Track Selection Cuts
- nHitsTPC ≥ 15
- nHitsFit_dedx ≥ 11
- High-tower energy > 0.5 GeV
- χ²_ee < 10
- Exclude χ²_AB < 10 (for π, K, p hypotheses — used for electron PID)

---

## Step 2: Histogram and Tree Filling

Reads the skimmed sDST files and fills analysis histograms and trees, including invariant mass spectra, rapidity distributions, and spin-sorted (Left/Right, Up/Down) candidate counts.

**For p+p:**
```bash
star-submit sDstAna.xml
```

**For p+Au:**
```bash
star-submit sDstAna_pA.xml
```

### Pair Selection
Trigger combinations handled:
- **p+p**: BEMC+BEMC, BEMC+EEMC, EEMC+EEMC
- **p+Au**: BEMC+BEMC only (electrons registered in barrel region only)

J/ψ signal window: **2.8 < M_ee < 3.2 GeV/c²**

---

## Step 3: Merge Output Files

After all grid jobs complete, merge the output ROOT files:

```bash
hadd merged_pp.root path/to/pp_output/*.root
hadd merged_pAu.root path/to/pAu_output/*.root
```

---

## Step 4: Final Plots and Asymmetry Extraction

Run the final plotting macro in ROOT to produce the invariant mass fits and asymmetry results:

```bash
root -l -b -q '__draw_mee_fast.C'
```

This macro:
- Subtracts combinatorial background (same-sign pairs for p+p; STARLight-based fit for p+Au)
- Fits the J/ψ peak with a Crystal Ball function convoluted with a 2nd-order polynomial
- Integrates signal counts in 2.8 < M_ee < 3.2 GeV/c² for each Left/Right and spin-Up/Down configuration
- Applies beam polarization correction (mean values: 59.8% for p+p, 60.1% for p+Au)
- Applies geometric correction via √(⟨cosφ⟩_L · ⟨cosφ⟩_R)
- Computes the final asymmetry:

```
A_N = ε / ( P · √( ⟨cosφ⟩_L · ⟨cosφ⟩_R ) )
```

### Background Estimates
- γAu contribution: ~12%
- γγ contribution: ~7%
- After background rejection: J/ψ signal purity ~86%

---

## Results

A_N results are plotted as a function of the photon-proton center-of-mass energy W_γp, computed as:

```
W²_γp = 2 · E_p · M_J/ψ · e^(−y)
```

Output figures include:
- Invariant mass spectra for all four spin/geometry configurations (Figs. 12, 13)
- A_N vs W_γp for p+p (Fig. 16)
- A_N vs W_γp for p+Au (Fig. 17)
- Combined A_N comparison (Fig. 18)

---

## Dependencies

- STAR software environment (STAR SL or equivalent)
- ROOT (for final macros)
- STARLight (Monte Carlo for background shape estimation)
- `star-submit` for grid job submission

---

## References

- J.P. Lansberg, L. Massacrier, L. Szymanowski, J. Wagner, *Single-Transverse-Spin Asymmetries in Exclusive Photo-production of J/ψ in Ultra-Peripheral Collisions*, Phys. Lett. B **793**, 33–40 (2019)

---

## Authors

Alexander Jentsch, Evgeny Shulga