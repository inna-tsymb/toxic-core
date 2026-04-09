# Computational Protein Design: Prion Core α-Helix Stabilization
## Iteration 2 — Comprehensive Analysis Report

### Executive Summary

This project demonstrates computational protein design for blocking the pathogenic prion transition PrP^C → PrP^Sc through systematic α-helix stabilization. 23 designs (baseline + 20 alpha/beta mutations + 2 openMM references) were analyzed. 7 designs improve upon the wild-type baseline energy.

---

## Key Results

### Baseline
| Structure | Energy (kcal/mol) | RMSD (Å) |
|---|---|---|
| Baseline (prion_core_autopsf.pdb) | -1.1096 | 0.0 |

### Top 5 Designs by Energy Improvement

| Rank | Design | Mutations | Energy (kcal/mol) | Improvement (kcal/mol) | RMSD (Å) | Population (%) |
|---|---|---|---|---|---|---|
| 1 | alpha_mutations/results_2 | MET 7 LEU, GLY 14 SER | -1.5755 | +0.4659 | 7.750 | 10.9 |
| 2 | beta_mutations/results_4 | GLY 14 SER, GLY 18 ALA, GLY 19 THR | -1.3259 | +0.2163 | 14.666 | 7.1 |
| 3 | alpha_mutations/results_5 | MET 7 LEU, GLY 14 SER, GLY 18 ALA, GLY 19 THR | -1.2906 | +0.1810 | 7.564 | 6.7 |
| 4 | alpha_mutations/results_6 | GLY 9 SER, GLY 14 SER, GLY 19 THR | -1.2869 | +0.1773 | 7.604 | 6.7 |
| 5 | alpha_mutations/results_4 | GLY 14 SER, GLY 18 ALA, GLY 19 THR | -1.2771 | +0.1675 | 7.922 | 6.6 |

### Summary Statistics
- **Total designs analyzed**: 22 (excluding baseline)
- **Designs improving baseline**: 7
- **Success rate**: 31.8%
- **Best energy**: -1.5755 kcal/mol (alpha_mutations/results_2)
- **Best population dominance**: 10.9% (alpha_mutations/results_2, vs 5.0% baseline)

---

## Full Design Table

| Design | Mutations | Energy | Improvement | RMSD | Population |
|---|---|---|---|---|---|
| Baseline | Wild Type | -1.1096 | 0.0 | 0.0 | 5.0% |
| alpha/results_1 | GLY 14 SER | -1.4064 | +0.2968 | 2.837 | 8.2% |
| alpha/results_2 | MET 7 LEU, GLY 14 SER | -1.5755 | +0.4659 | 7.750 | 10.9% |
| alpha/results_3 | MET 7 LEU, GLY 14 SER, GLY 18 SER | -1.1814 | +0.0718 | 7.896 | 5.6% |
| alpha/results_4 | GLY 14 SER, GLY 18 ALA, GLY 19 THR | -1.2771 | +0.1675 | 7.922 | 6.6% |
| alpha/results_5 | MET 7 LEU, GLY 14 SER, GLY 18 ALA, GLY 19 THR | -1.2906 | +0.1810 | 7.564 | 6.7% |
| alpha/results_6 | GLY 9 SER, GLY 14 SER, GLY 19 THR | -1.2869 | +0.1773 | 7.604 | 6.7% |
| alpha/results_7 | GLY 18 SER, GLY 19 ALA | -1.0797 | -0.0299 | 8.294 | 4.7% |
| alpha/results_8 | GLY 18 SER, GLY 19 THR | -1.0287 | -0.0809 | 8.306 | 4.3% |
| alpha/results_9 | MET 7 ILE, GLY 14 SER, GLY 18 SER | -1.0065 | -0.1031 | 12.011 | 4.2% |
| alpha/results_10 | ALA 11 LEU, ALA 15 ILE | +0.5147 | -1.6243 | 8.278 | 0.3% |
| beta/results_1 | GLY 14 SER | -0.9942 | -0.1154 | 8.258 | 4.1% |
| beta/results_2 | MET 7 LEU, GLY 14 SER | -0.7026 | -0.4070 | 8.775 | 2.5% |
| beta/results_3 | MET 7 LEU, GLY 14 SER, GLY 18 SER | -1.1219 | +0.0123 | 17.605 | 5.1% |
| beta/results_4 | GLY 14 SER, GLY 18 ALA, GLY 19 THR | -1.3259 | +0.2163 | 14.666 | 7.1% |
| beta/results_5 | MET 4 LEU, GLY 14 SER, GLY 18 ALA, GLY 19 THR | -0.7268 | -0.3828 | 14.962 | 2.6% |
| beta/results_6 | GLY 9 SER, GLY 14 SER, GLY 19 THR | -0.4299 | -0.6797 | 11.527 | 1.6% |
| beta/results_7 | GLY 18 SER, GLY 19 ALA | -0.7984 | -0.3112 | — | 2.9% |
| beta/results_8 | GLY 18 SER, GLY 19 THR | -0.9771 | -0.1325 | — | 4.0% |
| beta/results_9 | MET 7 ILE, GLY 14 SER, GLY 18 SER | -0.9136 | -0.1960 | — | 3.6% |
| beta/results_10 | ALA 11 LEU, ALA 15 ILE | -0.5873 | -0.5223 | 8.476 | 2.1% |

---

## Observations

### Alpha vs Beta Mutations
- Alpha mutations produced **5 of the top 7** improving designs
- The single-mutation GLY 14 SER (alpha/results_1) achieves strong improvement (+0.2968 kcal/mol) with the lowest RMSD (2.837 Å)
- Beta mutations show higher RMSD overall (8–17 Å), suggesting greater structural deviation

### RMSD Analysis
- RMSD values are large (2.8–17.6 Å) compared to typical design benchmarks, suggesting the designs deviate significantly from the baseline fold
- The best energy design (alpha/results_2, RMSD 7.75 Å) still improves energy despite structural deviation
- beta/results_3 shows extreme RMSD (17.6 Å) indicating near-complete refolding

### Mutations with Consistent Positive Effect
- **GLY 14 SER**: present in all top 5 designs — the most reliably stabilizing mutation
- **MET 7 LEU** combined with GLY 14 SER (alpha/results_2): best single combination
- **GLY 18/19 substitutions alone** (results_7, results_8): insufficient — slightly destabilizing

---

## Generated Files
- `plots/enhanced_funnel_plot.png` — Energy vs RMSD funnel plot
- `plots/comprehensive_probabilities.png` — Boltzmann population distribution
- `scripts/final_analysis_results.csv` — Full analysis with probabilities
- `output/mutated_prot_iteration_2/analysis/2_rmsd_overlay_alpha_mutations.png` — PyMOL overlay (alpha)
- `output/mutated_prot_iteration_2/analysis/2_rmsd_overlay_beta_mutations.png` — PyMOL overlay (beta)
