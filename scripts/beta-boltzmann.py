#!/usr/bin/env python3
"""
Boltzmann Probability Analysis for Prion Core Design Project
experimental data from Damietta & PyMOL
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def calculate_boltzmann_probabilities(energies, temperature=300):
    kB = 1.987e-3  # Boltzmann constant in kcal/mol/K
    kB_T = kB * temperature
    energies = np.array(energies)
    E_min = np.min(energies)
    E_relative = energies - E_min
    weights = np.exp(-E_relative / kB_T)
    Z = np.sum(weights)
    probabilities = weights / Z
    return {
        'probabilities': probabilities,
        'percentages': probabilities * 100
    }

def create_funnel_plot(rmsd_values, energies, design_labels, output_file='prion_funnel_plot_final.png'):
    plt.figure(figsize=(12, 9))
    
    # Використовуємо кольорову мапу для відображення енергії
    scatter = plt.scatter(rmsd_values, energies, c=energies, cmap='coolwarm', 
                          s=150, zorder=5, edgecolors='black', linewidth=1)
    
    # Підписуємо найважливіші точки, щоб не перевантажувати графік
    key_labels = ['WT', 'Result 3', 'dam_4', 'Result 4', 'Arom-Shield 3', 'Des 11']
    for i, (rmsd, energy, label) in enumerate(zip(rmsd_values, energies, design_labels)):
        if label in key_labels:
            plt.annotate(f'{label}\n({energy:.2f})', (rmsd, energy), 
                        xytext=(10, 5), textcoords='offset points', 
                        fontsize=10, ha='left', weight='bold')
    
    plt.xlabel('RMSD to Baseline Structure (Å)', fontsize=14)
    plt.ylabel('Alpha-Helix Total Energy (kcal/mol)', fontsize=14) 
    plt.title('Prion Core Folding Funnel: α-Helix Stabilization\n(All 18 Experimental Designs)', fontsize=16)
    plt.colorbar(scatter, label='Total Energy (kcal/mol)')
    plt.grid(True, alpha=0.3)
    
    # Лінія Дикого типу
    plt.axhline(y=energies[0], color='red', linestyle='--', alpha=0.5, label='WT Baseline Energy')
    
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Збережено графік: {output_file}")

def create_probability_bar_chart(probabilities, design_labels, output_file='prion_boltzmann_probabilities_final.png'):
    percentages = np.array(probabilities) * 100
    
    # Виводимо тільки Топ-6 найкращих, інші мають 0.0%
    top_indices = np.argsort(percentages)[::-1][:6]
    top_labels = [design_labels[i] for i in top_indices]
    top_pcts = [percentages[i] for i in top_indices]
    
    plt.figure(figsize=(10, 6))
    colors = plt.cm.viridis(np.linspace(0, 0.8, len(top_labels)))
    bars = plt.bar(top_labels, top_pcts, color=colors, alpha=0.8, edgecolor='black', linewidth=1)
    
    for bar, pct in zip(bars, top_pcts):
        height = bar.get_height()
        if pct > 0.1:
            plt.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                    f'{pct:.1f}%', ha='center', va='bottom', fontsize=12, weight='bold')
    
    plt.xlabel('Design State', fontsize=14)
    plt.ylabel('Population Percentage at 300K (%)', fontsize=14)
    plt.title('Boltzmann Population Distribution (Top Contenders)\nCompetition between Alpha-Helices', fontsize=16)
    plt.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Збережено графік: {output_file}")

# ДАНІ З DAMIETTA ТА PYMOL
design_labels = ['WT', 'Result 3', 'Result 2', 'dam_4', 'Result 4', 'Result 0', 'Result 1', 
                 'Arom-Shield 3', 'Arom-Acid 1', 'Arom-Acid 2', 'Arom-Shield 4', 'Result 8', 
                 'Des 8', 'Des 9', 'Trojan', 'Des 10', 'Des 12', 'Des 11']

energies = [-1.10, -7.64, -7.01, -7.45, -6.42, -6.76, -6.60, 
            -5.78, -4.56, -4.66, -4.06, -3.60, 
            -2.98, -2.54, -0.38, -1.29, -0.60, -0.63]

rmsd_values = [0.000, 0.116, 0.116, 0.192, 0.116, 0.116, 0.156, 
               3.140, 2.927, 2.932, 2.922, 0.273, 
               0.154, 0.153, 0.110, 2.711, 2.707, 3.254]

def main():
    print("="*60)
    print("PRION CORE α-HELIX STABILIZATION ANALYSIS (ACTUAL DATA)")
    print("="*60)
    
    results = calculate_boltzmann_probabilities(energies)
    create_funnel_plot(rmsd_values, energies, design_labels)
    create_probability_bar_chart(results['probabilities'], design_labels)
    
    print("\nТвій найкращий дизайн: Result 3")
    print(f"Енергетичний розрив (vs WT): {energies[0] - energies[1]:.2f} kcal/mol")
    print(f"Шанс утворення: {results['percentages'][1]:.1f}%")

if __name__ == "__main__":
    main()