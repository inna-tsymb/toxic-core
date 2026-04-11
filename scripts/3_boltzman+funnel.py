#!/usr/bin/env python3
"""
Updated Boltzmann Analysis for Complete Ensemble
Correctly distinguishes between WT Alpha and WT Beta baselines!
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as mpatches

def load_energy_data(csv_file=None):
    if csv_file is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.abspath(os.path.join(script_dir, '..'))
        csv_file = os.path.join(project_root, 'output', 'mutated_prot_iteration_2', 'analysis', '1_energy_analysis+rmsd.csv')
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
        analysis_dir = os.path.join(script_dir, '..', 'output', 'mutated_prot_iteration_2', 'analysis')
        os.makedirs(analysis_dir, exist_ok=True)
        output_file = os.path.join(analysis_dir, 'comprehensive_funnel_plot.png')
        
    plt.figure(figsize=(12, 9))
    
    added_alpha = False
    added_beta = False
    
    for idx, row in df.iterrows():
        rmsd = row['rmsd']
        energy = row['total_energy']
        design_id = str(row['design_id'])
        
        if pd.isna(rmsd):
            continue

        # 1. WT Alpha Baseline
        if 'prion_core' in design_id:
            plt.scatter(rmsd, energy, c='forestgreen', s=350, marker='*', zorder=5, 
                        edgecolors='black', linewidth=1.5, label='WT Alpha')
            plt.annotate('WT Alpha', (rmsd, energy), xytext=(-20, 15), textcoords='offset points', weight='bold', color='forestgreen')
            
        # 2. WT Beta Baseline
        elif 'prion_beta' in design_id:
            plt.scatter(rmsd, energy, c='darkorange', s=350, marker='*', zorder=5, 
                        edgecolors='black', linewidth=1.5, label='WT Beta')
            plt.annotate('WT Beta', (rmsd, energy), xytext=(10, 10), textcoords='offset points', weight='bold', color='darkorange')
            
        # 3. Alpha Designs
        elif 'alpha' in design_id.lower():
            label = 'Alpha Designs' if not added_alpha else ""
            plt.scatter(rmsd, energy, c='dodgerblue', s=120, marker='o', zorder=4, 
                        edgecolors='black', linewidth=1, alpha=0.8, label=label)
            added_alpha = True
            if row['energy_improvement'] > 0.05:
                plt.annotate(design_id.split('/')[-1].replace('results_', 'A'), 
                             (rmsd, energy), xytext=(8, -12), textcoords='offset points', fontsize=8)
                
        # 4. Beta Designs
        elif 'beta' in design_id.lower():
            label = 'Beta Designs (Collapsed)' if not added_beta else ""
            plt.scatter(rmsd, energy, c='crimson', s=120, marker='X', zorder=3, 
                        edgecolors='black', linewidth=1, alpha=0.7, label=label)
            added_beta = True
            plt.annotate(design_id.split('/')[-1].replace('results_', 'B'), 
                         (rmsd, energy), xytext=(8, 8), textcoords='offset points', fontsize=8)

    plt.xlabel('RMSD to Respective WT Structure (Å)', fontsize=14)
    plt.ylabel('Total Energy (kcal/mol)', fontsize=14)
    plt.title('Prion Core Folding Landscape\nResolution of Metastability', fontsize=16, weight='bold')
    plt.grid(True, alpha=0.3, linestyle='--')
    plt.legend(loc='upper left', fontsize=12, shadow=True)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Funnel plot saved as {output_file}")

def create_comprehensive_probability_plot(df, probabilities, output_file=None):
    if output_file is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        analysis_dir = os.path.join(script_dir, '..', 'output', 'mutated_prot_iteration_2', 'analysis')
        output_file = os.path.join(analysis_dir, 'ensemble_probabilities.png')
        
    plt.figure(figsize=(14, 7))
    percentages = probabilities * 100
    
    colors = []
    labels = []
    for design_id in df['design_id']:
        design_id_str = str(design_id)
        if 'prion_core' in design_id_str:
            colors.append('forestgreen')
            labels.append('WT Alpha')
        elif 'prion_beta' in design_id_str:
            colors.append('darkorange')
            labels.append('WT Beta')
        elif 'alpha' in design_id_str.lower():
            colors.append('dodgerblue')
            labels.append(design_id_str.split('/')[-1].replace('results_', 'A'))
        else:
            colors.append('crimson')
            labels.append(design_id_str.split('/')[-1].replace('results_', 'B'))

    bars = plt.bar(range(len(labels)), percentages, color=colors, alpha=0.8, edgecolor='black')
    
    for bar, pct in zip(bars, percentages):
        height = bar.get_height()
        if pct > 0.1:  
            plt.text(bar.get_x() + bar.get_width()/2., height + 0.5, f'{pct:.1f}%', 
                     ha='center', va='bottom', fontsize=10, weight='bold')

    plt.ylabel('Population Percentage at 300K (%)', fontsize=12)
    plt.title('Boltzmann Population Distribution', fontsize=16, weight='bold')
    plt.grid(True, alpha=0.3, axis='y')
    plt.xticks(range(len(labels)), labels, rotation=45, ha='right', fontsize=10)
    
    alpha_patch = mpatches.Patch(color='dodgerblue', label='Alpha Designs')
    beta_patch = mpatches.Patch(color='crimson', label='Beta Designs')
    wt_a_patch = mpatches.Patch(color='forestgreen', label='WT Alpha')
    wt_b_patch = mpatches.Patch(color='darkorange', label='WT Beta')
    plt.legend(handles=[alpha_patch, beta_patch, wt_a_patch, wt_b_patch], loc='upper right', fontsize=12)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Probability plot saved as {output_file}")

def main():
    print("Generating Comprehensive Plots for Alpha & Beta Ensemble...")
    df = load_energy_data()
    if df is None: return
    
    df = df.dropna(subset=['total_energy']).copy()
    df.reset_index(drop=True, inplace=True)
    
    boltzmann_results = calculate_boltzmann_probabilities(df['total_energy'].values)
    
    create_enhanced_funnel_plot(df)
    create_comprehensive_probability_plot(df, boltzmann_results['probabilities'])
    
    print("\n✅ Done! Check the 'plots' folder for the new images.")

if __name__ == "__main__":
    main()
