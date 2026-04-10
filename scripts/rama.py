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
