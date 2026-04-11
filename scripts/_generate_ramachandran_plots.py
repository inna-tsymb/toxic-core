#!/usr/bin/env python3
"""
Ramachandran Plot Generator for Prion Core Designs
Generates publication-quality plots showing backbone geometry validation
"""

import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, calc_dihedral
from adjustText import adjust_text
import os
import glob

def calculate_phi_psi(structure):
    """
    Calculate phi and psi angles for all residues in a structure
    """
    phi_psi = []
    residues = list(structure.get_residues())
    
    for i in range(1, len(residues) - 1):
        try:
            # Get atoms for phi calculation: C(i-1), N(i), CA(i), C(i)
            c_prev = residues[i-1]['C']
            n_curr = residues[i]['N'] 
            ca_curr = residues[i]['CA']
            c_curr = residues[i]['C']
            
            # Get atoms for psi calculation: N(i), CA(i), C(i), N(i+1)
            n_next = residues[i+1]['N']
            
            # Calculate dihedral angles
            phi = calc_dihedral(c_prev.get_vector(), n_curr.get_vector(),
                              ca_curr.get_vector(), c_curr.get_vector())
            psi = calc_dihedral(n_curr.get_vector(), ca_curr.get_vector(),
                              c_curr.get_vector(), n_next.get_vector())
            
            # Convert to degrees
            phi = np.degrees(phi)
            psi = np.degrees(psi)
            
            residue_name = residues[i].get_resname()
            residue_num = residues[i].get_id()[1]
            
            phi_psi.append({
                'residue': f"{residue_name}{residue_num}",
                'phi': phi,
                'psi': psi,
                'res_name': residue_name
            })
            
        except KeyError:
            # Skip if atoms are missing
            continue
            
    return phi_psi

def create_ramachandran_plot(phi_psi_data, title, filename, highlight_outliers=True):
    """
    Create a beautiful Ramachandran plot
    """
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Extract phi and psi values
    phi_values = [d['phi'] for d in phi_psi_data]
    psi_values = [d['psi'] for d in phi_psi_data]
    
    # Color code by residue type
    colors = []
    for data in phi_psi_data:
        res = data['res_name']
        if res == 'GLY':
            colors.append('red')       # Glycine - most flexible
        elif res in ['PRO']:
            colors.append('orange')    # Proline - constrained
        elif res in ['GLU', 'LYS']:
            colors.append('blue')      # Charged residues
        elif res in ['PHE', 'LEU', 'ILE']:
            colors.append('green')     # Hydrophobic
        else:
            colors.append('purple')    # Others
    
    # Create scatter plot
    scatter = ax.scatter(phi_values, psi_values, c=colors, s=100, alpha=0.8, 
                        edgecolors='black', linewidth=1)
    
    # Add Ramachandran regions (approximate)
    # Alpha-helix region
    alpha_phi = np.linspace(-80, -40, 100)
    alpha_psi_upper = -70 + 30 * np.sin((alpha_phi + 60) * np.pi / 40)
    alpha_psi_lower = -70 - 30 * np.sin((alpha_phi + 60) * np.pi / 40)
    ax.fill_between(alpha_phi, alpha_psi_lower, alpha_psi_upper, 
                   alpha=0.2, color='lightblue', label='α-helix region')
    
    # Beta-sheet region  
    beta_phi = np.linspace(-150, -90, 100)
    beta_psi_upper = 120 + 30 * np.sin((beta_phi + 120) * np.pi / 60)
    beta_psi_lower = 120 - 30 * np.sin((beta_phi + 120) * np.pi / 60)
    ax.fill_between(beta_phi, beta_psi_lower, beta_psi_upper,
                   alpha=0.2, color='lightgreen', label='β-sheet region')
    
    # Add non-overlapping labels
    texts = [ax.text(data['phi'], data['psi'], data['residue'], fontsize=7)
             for data in phi_psi_data]
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))
    
    # Formatting
    ax.set_xlabel('φ (Phi) angle (degrees)', fontsize=12)
    ax.set_ylabel('ψ (Psi) angle (degrees)', fontsize=12)
    ax.set_title(f'Ramachandran Plot: {title}', fontsize=14, weight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    
    # Add legend
    legend_elements = [
        plt.scatter([], [], c='red', s=100, label='Glycine'),
        plt.scatter([], [], c='blue', s=100, label='Charged (E,K)'),
        plt.scatter([], [], c='green', s=100, label='Hydrophobic'),
        plt.scatter([], [], c='purple', s=100, label='Other')
    ]
    ax.legend(handles=legend_elements, loc='upper right')
    
    # Add allowed/disallowed region legend
    ax.legend(loc='upper left', bbox_to_anchor=(0, 0.2))
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"Ramachandran plot saved: {filename}")
    
    return fig, ax

def create_comparison_plot(all_data, output_file='ramachandran_comparison.png', suptitle=None):
    """
    Create a grid comparison plot showing baseline vs all designs.
    Auto-sizes the grid to fit however many entries are in all_data.
    """
    n = len(all_data)
    ncols = 4
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 5, nrows * 4))
    axes = np.array(axes).flatten()

    for i, (name, data) in enumerate(all_data.items()):
        ax = axes[i]

        phi_values = [d['phi'] for d in data]
        psi_values = [d['psi'] for d in data]

        if 'baseline' in name.lower() or 'wt' in name.lower():
            color = 'crimson'
            alpha = 0.85
        else:
            color = 'steelblue'
            alpha = 0.6

        ax.scatter(phi_values, psi_values, c=color, s=50, alpha=alpha, edgecolors='none')

        ax.set_xlim(-180, 180)
        ax.set_ylim(-180, 180)
        ax.set_title(name, fontsize=9, weight='bold', wrap=True)
        ax.grid(True, alpha=0.3)
        ax.axhline(0, color='gray', lw=0.5)
        ax.axvline(0, color='gray', lw=0.5)

        row, col = divmod(i, ncols)
        if row == nrows - 1:
            ax.set_xlabel('φ (°)', fontsize=8)
        if col == 0:
            ax.set_ylabel('ψ (°)', fontsize=8)

    # Hide unused subplots
    for j in range(n, len(axes)):
        axes[j].set_visible(False)

    title = suptitle or 'Ramachandran Plot Comparison: Baseline vs Designed Variants'
    plt.suptitle(title, fontsize=14, weight='bold')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Comparison plot saved: {output_file}")

    return fig

def analyze_outliers(phi_psi_data):
    """
    Identify and analyze Ramachandran outliers
    """
    outliers = []
    
    for data in phi_psi_data:
        phi, psi = data['phi'], data['psi']
        res_name = data['res_name']
        
        is_outlier = False
        
        # Check if outside allowed regions (simplified)
        # Alpha region: phi ~-60, psi ~-45
        # Beta region: phi ~-120, psi ~+120
        # Left-handed alpha: phi ~+60, psi ~+45
        
        alpha_dist = np.sqrt((phi + 60)**2 + (psi + 45)**2)
        beta_dist = np.sqrt((phi + 120)**2 + (psi - 120)**2)
        left_dist = np.sqrt((phi - 60)**2 + (psi - 45)**2)
        
        min_dist = min(alpha_dist, beta_dist, left_dist)
        
        # For glycine, allow more flexibility
        threshold = 100 if res_name == 'GLY' else 60
        
        if min_dist > threshold:
            is_outlier = True
            outliers.append(data)
    
    return outliers

def main():
    """
    Main function to generate all Ramachandran plots for mutated_prot_iteration_2.
    Outputs to output/mutated_prot_iteration_2/analysis/ramachandran_plots/
    """
    import csv as csv_mod

    print("Generating Ramachandran Plots for Prion Core Designs (iteration 2)")
    print("=" * 60)

    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_root = os.path.join(project_root, 'output', 'mutated_prot_iteration_2')
    output_dir = os.path.join(output_root, 'analysis', 'ramachandran_plots')
    os.makedirs(output_dir, exist_ok=True)

    # Read mutation labels from energy CSV
    csv_path = os.path.join(output_root, 'analysis', '1_energy_analysis+rmsd.csv')
    label_map = {}  # design_id -> position_info
    with open(csv_path) as f:
        for row in csv_mod.DictReader(f):
            label_map[row['design_id']] = row['position_info']

    parser = PDBParser(QUIET=True)

    # Load baseline (WT, alpha structure)
    baseline_dir = os.path.join(output_root, 'prion_core_autopsf_openMM.pdb')
    baseline_pdbs = glob.glob(os.path.join(baseline_dir, '*.pdb'))
    if not baseline_pdbs:
        raise RuntimeError(f"Baseline PDB not found in {baseline_dir}")
    baseline_pdb = baseline_pdbs[0]
    print(f"Loading baseline: {baseline_pdb}")
    baseline_structure = parser.get_structure('baseline', baseline_pdb)
    baseline_phi_psi = calculate_phi_psi(baseline_structure)
    create_ramachandran_plot(
        baseline_phi_psi, 'Baseline (WT)',
        os.path.join(output_dir, 'ramachandran_baseline.png')
    )

    # Process each mutation group
    for group in ['alpha_mutations', 'beta_mutations']:
        print(f"\nProcessing {group}...")
        search_dir = os.path.join(output_root, group)

        # Collect result dirs sorted numerically
        result_dirs = sorted(
            [p for p in glob.glob(os.path.join(search_dir, 'results_*'))
             if os.path.isdir(p)],
            key=lambda p: int(os.path.basename(p).split('_')[1])
        )

        all_data = {'Baseline (WT)': baseline_phi_psi}

        for result_dir in result_dirs:
            result_name = os.path.basename(result_dir)
            design_id = f"{group}/{result_name}"
            position_info = label_map.get(design_id, result_name)
            short_title = f"{result_name}\n{position_info}"

            pdb_files = glob.glob(os.path.join(result_dir, '*.pdb'))
            if not pdb_files:
                print(f"  WARNING: no PDB in {result_dir}, skipping")
                continue
            pdb_file = pdb_files[0]

            print(f"  {design_id}: {position_info}")
            structure = parser.get_structure(result_name, pdb_file)
            phi_psi_data = calculate_phi_psi(structure)

            if not phi_psi_data:
                print(f"    WARNING: no dihedral angles calculated")
                continue

            all_data[short_title] = phi_psi_data
            plot_file = os.path.join(output_dir, f'ramachandran_{group}_{result_name}.png')
            create_ramachandran_plot(phi_psi_data, short_title, plot_file)

            outliers = analyze_outliers(phi_psi_data)
            print(f"    {len(outliers)} outlier(s)")

        # Grid comparison for this group
        group_label = group.replace('_', ' ').title()
        create_comparison_plot(
            all_data,
            output_file=os.path.join(output_dir, f'ramachandran_comparison_{group}.png'),
            suptitle=f'Ramachandran Comparison — {group_label} vs Baseline'
        )

    print("\n" + "=" * 60)
    print(f"Done. Plots written to:\n  {output_dir}")

if __name__ == "__main__":
    main()
