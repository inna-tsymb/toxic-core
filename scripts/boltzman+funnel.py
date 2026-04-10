#!/usr/bin/env python3
"""
Updated Boltzmann Analysis for Complete Ensemble (Alpha & Beta)
Visually separates stable Alpha designs from collapsed Beta designs.
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def load_energy_data(csv_file=None):
    if csv_file is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        csv_file = os.path.join(script_dir, 'energy_analysis.csv')
    try:
        return pd.read_csv(csv_file)
    except FileNotFoundError:
        print(f"Error: {csv_file} not found.")
        return None

def calculate_boltzmann_probabilities(energies, temperature=300):
    kB = 1.987e-3  
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

def create_enhanced_funnel_plot(df, output_file=None):
    if output_file is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        plots_dir = os.path.join(script_dir, '..', 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        output_file = os.path.join(plots_dir, 'comprehensive_funnel_plot.png')
        
    plt.figure(figsize=(12, 9))
    
    # Створюємо списки для легенди, щоб не дублювати
    added_alpha = False
    added_beta = False
    
    for idx, row in df.iterrows():
        rmsd = row['rmsd']
        energy = row['total_energy']
        design_id = str(row['design_id'])
        
        if pd.isna(rmsd):
            continue

        # Логіка розфарбовування (Альфа - сині, Бета - червоні, Бейслайн - зелена зірка)
        if 'Baseline' in design_id or 'wt_' in design_id.lower() or design_id == 'prion_core_autopsf_openMM.pdb':
            plt.scatter(rmsd, energy, c='green', s=300, marker='*', zorder=5, 
                        edgecolors='black', linewidth=1.5, label='WT Baseline')
            plt.annotate('Baseline', (rmsd, energy), xytext=(10, 10), textcoords='offset points', weight='bold')
            
        elif 'alpha' in design_id.lower():
            label = 'Alpha Designs (Stable)' if not added_alpha else ""
            plt.scatter(rmsd, energy, c='dodgerblue', s=120, marker='o', zorder=4, 
                        edgecolors='black', linewidth=1, alpha=0.8, label=label)
            added_alpha = True
            if row['energy_improvement'] > 0.05:
                plt.annotate(design_id.split('/')[-1].replace('results_', 'A'), 
                             (rmsd, energy), xytext=(8, -12), textcoords='offset points', fontsize=8)
                
        elif 'beta' in design_id.lower():
            label = 'Beta Designs (Collapsed)' if not added_beta else ""
            plt.scatter(rmsd, energy, c='crimson', s=120, marker='X', zorder=3, 
                        edgecolors='black', linewidth=1, alpha=0.7, label=label)
            added_beta = True
            plt.annotate(design_id.split('/')[-1].replace('results_', 'B'), 
                         (rmsd, energy), xytext=(8, 8), textcoords='offset points', fontsize=8, color='darkred')

    plt.xlabel('RMSD to WT Structure (Å) - Shows Structural Deviation', fontsize=14)
    plt.ylabel('Total Energy (kcal/mol)', fontsize=14)
    plt.title('Prion Core Folding Landscape\nAlpha Stabilization vs. Beta Collapse', fontsize=16, weight='bold')
    plt.grid(True, alpha=0.3, linestyle='--')
    plt.legend(loc='upper left', fontsize=12, shadow=True)
    
    # Текстова врізка з поясненням
    plt.text(0.65, 0.95, "Notice the severe structural collapse\n(High RMSD) of Beta designs,\nproving negative design success.", 
             transform=plt.gca().transAxes, fontsize=11, bbox=dict(boxstyle='round', facecolor='whitesmoke', alpha=0.8))

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Funnel plot saved as {output_file}")

def create_comprehensive_probability_plot(df, probabilities, output_file=None):
    if output_file is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        plots_dir = os.path.join(script_dir, '..', 'plots')
        output_file = os.path.join(plots_dir, 'ensemble_probabilities.png')
        
    plt.figure(figsize=(14, 7))
    percentages = probabilities * 100
    
    # Кольори для стовпців
    colors = []
    labels = []
    for design_id in df['design_id']:
        design_id_str = str(design_id)
        if 'Baseline' in design_id_str or 'wt_' in design_id_str.lower():
            colors.append('green')
            labels.append('WT Baseline')
        elif 'alpha' in design_id_str.lower():
            colors.append('dodgerblue')
            labels.append(design_id_str.split('/')[-1].replace('results_', 'Alpha '))
        else:
            colors.append('crimson')
            labels.append(design_id_str.split('/')[-1].replace('results_', 'Beta '))

    bars = plt.bar(range(len(labels)), percentages, color=colors, alpha=0.8, edgecolor='black')
    
    for bar, pct in zip(bars, percentages):
        height = bar.get_height()
        if pct > 0.1:  # Показуємо цифри тільки якщо ймовірність більше 0.1%
            plt.text(bar.get_x() + bar.get_width()/2., height + 0.5, f'{pct:.1f}%', 
                     ha='center', va='bottom', fontsize=10, weight='bold')

    plt.ylabel('Population Percentage at 300K (%)', fontsize=12)
    plt.title('Boltzmann Population Distribution of the Entire Ensemble', fontsize=16, weight='bold')
    plt.grid(True, alpha=0.3, axis='y')
    plt.xticks(range(len(labels)), labels, rotation=45, ha='right', fontsize=10)
    
    # Додаємо кастомну легенду
    import matplotlib.patches as mpatches
    alpha_patch = mpatches.Patch(color='dodgerblue', label='Alpha Designs')
    beta_patch = mpatches.Patch(color='crimson', label='Beta Designs')
    wt_patch = mpatches.Patch(color='green', label='Wild Type')
    plt.legend(handles=[alpha_patch, beta_patch, wt_patch], loc='upper right', fontsize=12)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Probability plot saved as {output_file}")

def main():
    print("Generating Comprehensive Plots for Alpha & Beta Ensemble...")
    df = load_energy_data()
    if df is None: return
    
    # Drop rows with NaN energy to prevent math errors
    df = df.dropna(subset=['total_energy']).copy()
    df.reset_index(drop=True, inplace=True)
    
    boltzmann_results = calculate_boltzmann_probabilities(df['total_energy'].values)
    
    create_enhanced_funnel_plot(df)
    create_comprehensive_probability_plot(df, boltzmann_results['probabilities'])
    
    print("\n✅ Done! Check the 'plots' folder for the new images.")

if __name__ == "__main__":
    main()