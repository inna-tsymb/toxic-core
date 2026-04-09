# Prion Core α-Helix Stabilization: Report Recommendations

## 1. Introduction

### Background
The prion protein core region (residues 106–126, sequence `KTNMKHMAGAAAAGAVVGGLG`) is prone to misfolding from the healthy α-helical PrP^C conformation to the pathogenic β-sheet PrP^Sc state. Key structural vulnerabilities:
- **Position 7 (MET)**: oxidation susceptibility
- **Position 9, 14, 18, 19 (GLY)**: excessive backbone flexibility promoting β-sheet formation
- **Positions 11, 15 (ALA)**: weak helix propensity

### Objective
Systematically mutate flexibility hotspots to thermodynamically stabilize the α-helical conformation and block the PrP^C → PrP^Sc transition.

---

## 2. Methods

### Structure Preparation
- Baseline: `prion_core_autopsf.pdb`, energy = -1.1096 kcal/mol (Damietta SE scoring)
- VMD CHARMM typing for hydrogen addition

### Mutation Strategy
Guided by ProteinMPNN favorability heatmaps. Two parallel tracks:
- **Alpha mutations**: target GLY flexibility hotspots with helix-favoring substitutions (SER, ALA, THR)
- **Beta mutations**: same positions, scored against β-sheet baseline

### Analysis Pipeline
1. `1_extract_damietta_energies.py` — extract energies from Damietta output
2. `2_enhanced_pymol_rmsd.pml` — calculate RMSD per design vs baseline
3. `3_enhanced_boltzmann_analysis.py` — Boltzmann population + funnel plot

---

## 3. Results

### Energy Stabilization
- **7/22 designs (31.8%)** improve upon baseline
- **Best design**: alpha/results_2 (MET 7 LEU + GLY 14 SER), -1.5755 kcal/mol (+0.4659 improvement)
- **Most destabilizing**: alpha/results_10 (ALA 11 LEU + ALA 15 ILE), +0.5147 kcal/mol

### Boltzmann Population (300K)
- Baseline population: ~5.0%
- Best design population: 10.9% (2.2× improvement)
- Bulky hydrophobic substitutions at ALA 11/15 strongly destabilize

### Structural Integrity (RMSD)
- RMSD range: 2.8 – 17.6 Å
- Lowest RMSD: alpha/results_1 (GLY 14 SER only), 2.837 Å
- High RMSDs suggest Damietta is sampling conformations beyond the local α-helical basin

### Key Finding
**GLY 14 SER** is the most reliably stabilizing single mutation — present in all top 5 designs. MET 7 LEU adds further stabilization when combined with GLY 14 SER.

---

## 4. Discussion

### Design Strategy Evaluation
- Alpha-track mutations consistently outperform beta-track — the α-helical context scores better
- Multi-mutation designs (3–4 mutations) do not always outperform simpler 1–2 mutation designs
- ALA substitutions (results_10) are strongly counterproductive — likely disrupt packing

### Limitations
- Large RMSD values (7–17 Å) indicate designs may not preserve the target fold
- Damietta energy is a single-conformation score, not a full free energy
- Only 3 beta designs missing RMSD (beta/results_7–9) — run PyMOL script to complete

### Next Steps
1. Run `2_enhanced_pymol_rmsd.pml` to fill in missing RMSD for beta/results_7–9
2. Filter designs by RMSD < 5 Å for structural viability before final selection
3. Consider MD simulation on top candidates to validate thermodynamic stability
4. Experimental validation (CD spectroscopy, ThT aggregation assay)

---

## 5. Conclusions

- Computational design successfully identifies stabilizing mutations for prion core α-helix
- GLY 14 SER is the most robust single mutation
- MET 7 LEU + GLY 14 SER (alpha/results_2) is the recommended candidate for further study
- 31.8% success rate is reasonable for a first combinatorial iteration; iteration 3 should focus on GLY 14 SER combinations with lower RMSD
