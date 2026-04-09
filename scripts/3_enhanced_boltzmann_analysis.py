#!/usr/bin/env python3
"""
Updated Boltzmann Analysis for 10 Successful Mutations
Complete analysis with visualization generation
"""

import sys
import os

# Ensure we're using the virtual environment's Python
venv_python = os.path.join(os.path.dirname(__file__), '..', '.venv', 'bin', 'python')
venv_python_real = os.path.realpath(venv_python)
current_real = os.path.realpath(sys.executable)
if current_real != venv_python_real and os.path.exists(venv_python):
    print(f"Re-executing with virtual environment Python: {venv_python}")
    os.execv(venv_python, [venv_python] + sys.argv)

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def load_energy_data(csv_file=None):
    """
    Load energy and RMSD data from CSV file
    """
    if csv_file is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        csv_file = os.path.join(script_dir, 'energy_analysis.csv')
    try:
        df = pd.read_csv(csv_file)
        return df
    except FileNotFoundError:
        print(f"Error: {csv_file} not found. Please run the energy extraction script first.")
        return None

def calculate_boltzmann_probabilities(energies, temperature=300):
    """
    Calculate Boltzmann probabilities for ensemble of states
    """
    kB = 1.987e-3  # Boltzmann constant in kcal/mol/K
    kB_T = kB * temperature
    
    energies = np.array(energies)
    
    # Numerical stability: subtract minimum energy
    E_min = np.min(energies)
    E_relative = energies - E_min
    
    # Calculate Boltzmann weights
    weights = np.exp(-E_relative / kB_T)
    
    # Partition function
    Z = np.sum(weights)
    
    # Probabilities
    probabilities = weights / Z
    
    return {
        'energies': energies,
        'relative_energies': E_relative,
        'weights': weights,
        'partition_function': Z,
        'probabilities': probabilities,
        'percentages': probabilities * 100
    }

def create_enhanced_funnel_plot(df, output_file=None):
    """
    Create enhanced funnel plot for all designs
    """
    if output_file is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        plots_dir = os.path.join(script_dir, '..', 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        output_file = os.path.join(plots_dir, 'enhanced_funnel_plot.png')
    plt.figure(figsize=(12, 9))
    
    # Extract data
    rmsd_values = df['rmsd'].astype(float)
    energies = df['total_energy']
    design_labels = df['design_id']
    
    # Create color map for different designs
    n_designs = len(df) - 1
    design_colors = plt.cm.viridis(np.linspace(0, 1, max(n_designs, 1)))
    
    # Plot baseline separately
    baseline_idx = df[df['design_id'] == 'Baseline'].index[0]
    plt.scatter(0.0, energies[baseline_idx], 
               c='red', s=150, marker='s', label='Baseline', 
               zorder=5, edgecolors='black', linewidth=2)
    
    # Plot all designs with known RMSD values
    missing_rmsd = 0
    color_index = 0
    for idx, row in df.iterrows():
        if row['design_id'] == 'Baseline':
            continue

        rmsd = row['rmsd']
        if pd.isna(rmsd):
            missing_rmsd += 1
            continue

        plt.scatter(rmsd, row['total_energy'], 
                   c=[design_colors[color_index]], s=120, marker='o', 
                   zorder=4, edgecolors='black', linewidth=1,
                   alpha=0.8)
        color_index += 1
            
        # Add labels for significant designs
        if row['energy_improvement'] > 0.05:  # Significant improvement
            plt.annotate(f"{row['design_id']}\n({row['total_energy']:.3f})", 
                       (rmsd, row['total_energy']), 
                       xytext=(10, 10), textcoords='offset points', 
                       fontsize=9, ha='left')

    if missing_rmsd > 0:
        plt.text(0.02, 0.06, f'⚠️ {missing_rmsd} designs have missing RMSD values.',
                 transform=plt.gca().transAxes, fontsize=11,
                 bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.4))

    # Define x-axis range
    valid_rmsds = rmsd_values.dropna()
    if len(valid_rmsds) > 0:
        max_rmsd = np.nanmax(valid_rmsds) * 1.1
    else:
        max_rmsd = 1.0
    plt.xlim(-0.5, max_rmsd)

    # Add baseline energy line
    baseline_energy = energies[baseline_idx]
    plt.axhline(y=baseline_energy, color='red', linestyle='--', 
               alpha=0.7, label='Baseline Energy')
    
    # Formatting
    plt.xlabel('RMSD to Baseline Structure (Å)', fontsize=14)
    plt.ylabel('Total Energy (kcal/mol)', fontsize=14)
    plt.title('Enhanced Prion Core Folding Funnel\nα-Helix Stabilization with 10 Successful Mutations', 
              fontsize=16, weight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend(loc='best')
    
    # Add success annotation
    improved_count = len(df[df['energy_improvement'] > 0]) - 1  # Exclude baseline
    if improved_count > 0:
        plt.text(0.02, 0.98, f'✓ {improved_count} designs improve upon baseline', 
                transform=plt.gca().transAxes, fontsize=12, 
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Enhanced funnel plot saved as {output_file}")

def create_comprehensive_probability_plot(probabilities, design_labels, improvements, 
                                        output_file=None):
    """
    Create comprehensive probability distribution plot
    """
    if output_file is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        plots_dir = os.path.join(script_dir, '..', 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        output_file = os.path.join(plots_dir, 'comprehensive_probabilities.png')
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 12))
    
    percentages = np.array(probabilities) * 100
    
    # Top plot: Probability distribution
    colors = ['red'] + ['steelblue'] * (len(percentages) - 1)
    bars1 = ax1.bar(range(len(design_labels)), percentages, color=colors, alpha=0.8,
                    edgecolor='black', linewidth=1)
    
    # Add percentage labels on bars
    for i, (bar, pct) in enumerate(zip(bars1, percentages)):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                f'{pct:.1f}%', ha='center', va='bottom', 
                fontsize=10, weight='bold')
    
    ax1.set_ylabel('Population Percentage at 300K (%)', fontsize=12)
    ax1.set_title('Boltzmann Population Distribution', fontsize=14, weight='bold')
    ax1.grid(True, alpha=0.3, axis='y')
    ax1.set_xticks(range(len(design_labels)))
    ax1.set_xticklabels(design_labels, rotation=45, ha='right')
    
    # Bottom plot: Energy improvements
    improvement_colors = ['gray'] + ['green' if imp > 0 else 'orange' for imp in improvements[1:]]
    bars2 = ax2.bar(range(len(design_labels)), improvements, color=improvement_colors, 
                    alpha=0.8, edgecolor='black', linewidth=1)
    
    # Add improvement labels
    for i, (bar, imp) in enumerate(zip(bars2, improvements)):
        if imp != 0:  # Skip baseline
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., 
                    height + (0.01 if height > 0 else -0.02),
                    f'{imp:+.3f}', ha='center', 
                    va='bottom' if height > 0 else 'top',
                    fontsize=9, weight='bold')
    
    ax2.set_ylabel('Energy Improvement (kcal/mol)', fontsize=12)
    ax2.set_xlabel('Design ID', fontsize=12)
    ax2.set_title('Energy Improvements vs Baseline', fontsize=14, weight='bold')
    ax2.grid(True, alpha=0.3, axis='y')
    ax2.axhline(y=0, color='black', linestyle='-', alpha=0.5)
    ax2.set_xticks(range(len(design_labels)))
    ax2.set_xticklabels(design_labels, rotation=45, ha='right')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Comprehensive probability plot saved as {output_file}")

def generate_final_report(df, boltzmann_results):
    """
    Generate comprehensive analysis report
    """
    print("\n" + "="*70)
    print("COMPREHENSIVE PRION CORE MUTATION ANALYSIS")
    print("="*70)
    
    # Basic statistics
    total_designs = len(df) - 1  # Exclude baseline
    improved_designs = len(df[df['energy_improvement'] > 0]) - 1  # Exclude baseline
    
    print(f"Total designs analyzed: {total_designs}")
    print(f"Designs improving baseline: {improved_designs}")
    print(f"Success rate: {(improved_designs/total_designs)*100:.1f}%")
    
    # Best design analysis
    best_idx = df['energy_improvement'].idxmax()
    best_design = df.loc[best_idx]
    
    print(f"\n🏆 BEST DESIGN:")
    print(f"Design: {best_design['design_id']}")
    print(f"Mutations: {best_design['position_info']}")
    print(f"Energy: {best_design['total_energy']:.4f} kcal/mol")
    print(f"Improvement: {best_design['energy_improvement']:+.4f} kcal/mol")
    print(f"RMSD: {best_design['rmsd']:.3f} Å")
    
    # Population analysis
    probabilities = boltzmann_results['probabilities']
    best_population = probabilities[best_idx] * 100
    baseline_population = probabilities[0] * 100
    
    print(f"Population dominance: {best_population:.1f}%")
    print(f"Population improvement: {best_population/baseline_population:.1f}x increase")
    
    # Summary table
    print(f"\n" + "="*70)
    print("TOP 5 DESIGNS")
    print("="*70)
    top_designs = df.nlargest(6, 'energy_improvement')[1:]  # Exclude baseline, get top 5
    for i, (_, design) in enumerate(top_designs.iterrows(), 1):
        pop_pct = probabilities[design.name] * 100
        print(f"{i}. {design['design_id']:<10} | Energy: {design['total_energy']:>8.4f} | "
              f"Improvement: {design['energy_improvement']:>+7.4f} | Population: {pop_pct:>5.1f}%")

def main():
    """
    Main analysis workflow
    """
    print("Enhanced Boltzmann Analysis for 10 Successful Mutations")
    print("=" * 60)
    
    # Load data
    df = load_energy_data()
    if df is None:
        print("Cannot proceed without energy data. Please:")
        print("1. Run extract_damietta_energies.py first")
        print("2. Calculate RMSD values using PyMOL")
        print("3. Update energy_analysis.csv with RMSD data")
        return
    
    print(f"Loaded data for {len(df)} designs (including baseline)")
    
    # Check for RMSD data
    if df['rmsd'].isna().any():
        print("⚠️  Warning: Some RMSD values are missing. Please complete PyMOL analysis.")
        print("Those designs will be omitted from the funnel plot until RMSD is added.")
    
    # Calculate Boltzmann probabilities
    energies = df['total_energy'].values
    boltzmann_results = calculate_boltzmann_probabilities(energies)
    
    # Add probabilities to dataframe
    df['population_percent'] = boltzmann_results['percentages']
    
    # Generate visualizations
    create_enhanced_funnel_plot(df)
    create_comprehensive_probability_plot(
        boltzmann_results['probabilities'], 
        df['design_id'].values,
        df['energy_improvement'].values
    )
    
    # Generate final report
    generate_final_report(df, boltzmann_results)
    
    # Save enhanced results
    script_dir = os.path.dirname(os.path.abspath(__file__))
    final_csv_path = os.path.join(script_dir, 'final_analysis_results.csv')
    df.to_csv(final_csv_path, index=False)
    print(f"\n📊 Complete results saved to: {final_csv_path}")

if __name__ == "__main__":
    main()
