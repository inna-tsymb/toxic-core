#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, calc_dihedral
from adjustText import adjust_text
import os
import glob

def calculate_phi_psi(structure):
    phi_psi = []
    residues = list(structure.get_residues())
    for i in range(1, len(residues) - 1):
        try:
            c_prev = residues[i-1]['C']
            n_curr = residues[i]['N']
            ca_curr = residues[i]['CA']
            c_curr = residues[i]['C']
            n_next = residues[i+1]['N']
            phi = np.degrees(calc_dihedral(c_prev.get_vector(), n_curr.get_vector(),
                                           ca_curr.get_vector(), c_curr.get_vector()))
            psi = np.degrees(calc_dihedral(n_curr.get_vector(), ca_curr.get_vector(),
                                           c_curr.get_vector(), n_next.get_vector()))
            phi_psi.append({'residue': f"{residues[i].get_resname()}{residues[i].get_id()[1]}",
                            'phi': phi, 'psi': psi, 'res_name': residues[i].get_resname()})
        except KeyError:
            continue
    return phi_psi

def create_ramachandran_plot(phi_psi_data, title, filename):
    fig, ax = plt.subplots(figsize=(10, 8))
    phi_values = [d['phi'] for d in phi_psi_data]
    psi_values = [d['psi'] for d in phi_psi_data]
    color_map = {'GLY': 'red', 'PRO': 'orange'}
    charged = {'GLU', 'LYS'}
    hydrophobic = {'PHE', 'LEU', 'ILE'}
    colors = [color_map.get(d['res_name'], 'blue' if d['res_name'] in charged
              else 'green' if d['res_name'] in hydrophobic else 'purple')
              for d in phi_psi_data]
    ax.scatter(phi_values, psi_values, c=colors, s=100, alpha=0.8, edgecolors='black', linewidth=1)
    alpha_phi = np.linspace(-80, -40, 100)
    ax.fill_between(alpha_phi, -70 - 30*np.sin((alpha_phi+60)*np.pi/40),
                    -70 + 30*np.sin((alpha_phi+60)*np.pi/40), alpha=0.2, color='lightblue', label='α-helix region')
    beta_phi = np.linspace(-150, -90, 100)
    ax.fill_between(beta_phi, 120 - 30*np.sin((beta_phi+120)*np.pi/60),
                    120 + 30*np.sin((beta_phi+120)*np.pi/60), alpha=0.2, color='lightgreen', label='β-sheet region')
    texts = [ax.text(d['phi'], d['psi'], d['residue'], fontsize=7) for d in phi_psi_data]
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))
    ax.set_xlabel('φ (Phi) angle (degrees)', fontsize=12)
    ax.set_ylabel('ψ (Psi) angle (degrees)', fontsize=12)
    ax.set_title(f'Ramachandran Plot: {title}', fontsize=14, weight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.legend(loc='upper left', bbox_to_anchor=(0, 0.2))
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Ramachandran plot saved: {filename}")

def create_comparison_plot(all_data, output_file='ramachandran_comparison.png', suptitle=None):
    n = len(all_data)
    ncols = 4
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 5, nrows * 4))
    axes = np.array(axes).flatten()
    for i, (name, data) in enumerate(all_data.items()):
        ax = axes[i]
        phi_values = [d['phi'] for d in data]
        psi_values = [d['psi'] for d in data]
        color = 'crimson' if ('wt' in name.lower() or 'baseline' in name.lower()) else 'steelblue'
        ax.scatter(phi_values, psi_values, c=color, s=50, alpha=0.7, edgecolors='none')
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
    for j in range(n, len(axes)):
        axes[j].set_visible(False)
    plt.suptitle(suptitle or 'Ramachandran Comparison', fontsize=14, weight='bold')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Comparison plot saved: {output_file}")

def analyze_outliers(phi_psi_data):
    outliers = []
    for data in phi_psi_data:
        phi, psi, res_name = data['phi'], data['psi'], data['res_name']
        min_dist = min(np.sqrt((phi+60)**2+(psi+45)**2),
                       np.sqrt((phi+120)**2+(psi-120)**2),
                       np.sqrt((phi-60)**2+(psi-45)**2))
        if min_dist > (100 if res_name == 'GLY' else 60):
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
    label_map = {}
    if os.path.exists(csv_path):
        with open(csv_path) as f:
            for row in csv_mod.DictReader(f):
                label_map[row['design_id']] = row['position_info']

    parser = PDBParser(QUIET=True)

    # ==========================================
    # 1. ЗАВАНТАЖУЄМО ОБИДВІ БАЗИ (ALPHA ТА BETA)
    # ==========================================
    
    # Alpha Baseline
    alpha_dir = os.path.join(output_root, 'prion_core_autopsf_openMM.pdb')
    alpha_pdbs = glob.glob(os.path.join(alpha_dir, '*.pdb'))
    if not alpha_pdbs:
        raise RuntimeError(f"Alpha Baseline PDB not found in {alpha_dir}")
    print(f"Loading Alpha baseline: {alpha_pdbs[0]}")
    alpha_structure = parser.get_structure('alpha_wt', alpha_pdbs[0])
    alpha_phi_psi = calculate_phi_psi(alpha_structure)
    create_ramachandran_plot(alpha_phi_psi, 'WT Alpha Baseline', os.path.join(output_dir, 'ramachandran_baseline_alpha.png'))

    # Beta Baseline
    beta_dir = os.path.join(output_root, 'prion_beta_autopsf_openMM.pdb')
    beta_pdbs = glob.glob(os.path.join(beta_dir, '*.pdb'))
    if not beta_pdbs:
        raise RuntimeError(f"Beta Baseline PDB not found in {beta_dir}")
    print(f"Loading Beta baseline: {beta_pdbs[0]}")
    beta_structure = parser.get_structure('beta_wt', beta_pdbs[0])
    beta_phi_psi = calculate_phi_psi(beta_structure)
    create_ramachandran_plot(beta_phi_psi, 'WT Beta Baseline', os.path.join(output_dir, 'ramachandran_baseline_beta.png'))


    # ==========================================
    # 2. ОБРОБКА ГРУП МУТАНТІВ
    # ==========================================
    for group in ['alpha_mutations', 'beta_mutations']:
        print(f"\nProcessing {group}...")
        search_dir = os.path.join(output_root, group)

        result_dirs = sorted(
            [p for p in glob.glob(os.path.join(search_dir, 'results_*')) if os.path.isdir(p)],
            key=lambda p: int(os.path.basename(p).split('_')[1])
        )

        # Вибираємо правильну базу для порівняльної сітки (Grid Plot)
        if group == 'alpha_mutations':
            current_baseline_data = alpha_phi_psi
            baseline_title = 'WT Alpha'
        else:
            current_baseline_data = beta_phi_psi
            baseline_title = 'WT Beta'

        all_data = {baseline_title: current_baseline_data}

        for result_dir in result_dirs:
            result_name = os.path.basename(result_dir)
            design_id = f"{group}/{result_name}"
            position_info = label_map.get(design_id, result_name)
            short_title = f"{result_name}\n{position_info}"

            pdb_files = glob.glob(os.path.join(result_dir, '*.pdb'))
            if not pdb_files:
                continue
            pdb_file = pdb_files[0]

            structure = parser.get_structure(result_name, pdb_file)
            phi_psi_data = calculate_phi_psi(structure)

            if not phi_psi_data:
                continue

            all_data[short_title] = phi_psi_data
            plot_file = os.path.join(output_dir, f'ramachandran_{group}_{result_name}.png')
            create_ramachandran_plot(phi_psi_data, short_title, plot_file)

            outliers = analyze_outliers(phi_psi_data)
            print(f"  {design_id}: {len(outliers)} outlier(s)")

        # Grid comparison for this group using the correct baseline
        group_label = group.replace('_', ' ').title()
        create_comparison_plot(
            all_data,
            output_file=os.path.join(output_dir, f'ramachandran_comparison_{group}.png'),
            suptitle=f'Ramachandran Comparison — {group_label} vs {baseline_title}'
        )

    print("\n" + "=" * 60)
    print(f"Done. Plots written to:\n  {output_dir}")

if __name__ == "__main__":
    main()
