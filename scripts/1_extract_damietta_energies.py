#!/usr/bin/env python3
"""
Extract Energy Data from Damietta Output Files
Processes your mutation results for Boltzmann analysis
"""

import pandas as pd
import json
import os
import glob

def parse_result_csv(csv_file):
    """
    Parse the result.csv file from Damietta output
    """
    try:
        df = pd.read_csv(csv_file)
        return df
    except Exception as e:
        print(f"Error reading {csv_file}: {e}")
        return None

def extract_energies_from_directory(output_dir):
    """
    Extract energy data from a Damietta output tree by finding all nested result.csv files.
    """
    results = []
    csv_files = glob.glob(os.path.join(output_dir, '**', 'result.csv'), recursive=True)

    if not csv_files:
        print(f"No result.csv files found under {output_dir}")
        return results

    for csv_file in sorted(csv_files):
        relative_dir = os.path.relpath(os.path.dirname(csv_file), output_dir)
        print(f"Processing {csv_file}...")
        df = parse_result_csv(csv_file)

        if df is None:
            continue

        for index, row in df.iterrows():
            design_id = f"{relative_dir}"
            if len(df) > 1:
                design_id = f"{relative_dir}_{index+1}"

            total_energy = row.get('total', row.get('ΔGtotal', None))
            result_entry = {
                'design_id': design_id,
                'filename': row.get('Result File', f'result{index}.pdb'),
                'position_info': row.get('Position', 'Unknown'),
                'total_energy': total_energy,
                'pp_energy': row.get('pp_dG', row.get('ΔGpp', None)),
                'k_energy': row.get('k_dG', row.get('ΔGk', None)),
                'lj_energy': row.get('lj', row.get('ΔGLJ', None)),
                'solv_energy': row.get('solv', row.get('ΔGsolv', None)),
                'elec_energy': row.get('elec', row.get('ΔGelec', None))
            }
            results.append(result_entry)

    return results

def create_energy_summary(results, baseline_energy, output_file='energy_analysis.csv'):
    """
    Create comprehensive energy summary
    """
    # Add baseline as first entry
    baseline_entry = {
        'design_id': 'Baseline',
        'filename': 'prion_core_autopsf.pdb',
        'position_info': 'Wild Type',
        'total_energy': baseline_energy,
        'pp_energy': None,
        'k_energy': None,
        'lj_energy': None,
        'solv_energy': None,
        'elec_energy': None,
        'energy_improvement': 0.0,
        'rmsd': 0.0
    }
    
    # Process mutation results
    summary_data = [baseline_entry]
    
    for result in results:
        total_energy = result['total_energy']
        if total_energy is not None:
            improvement = baseline_energy - total_energy
            summary_entry = {
                **result,
                'energy_improvement': improvement,
                'rmsd': None  # To be filled from PyMOL
            }
            summary_data.append(summary_entry)
    
    # Create DataFrame and save
    df = pd.DataFrame(summary_data)
    df.to_csv(output_file, index=False)
    
    print(f"Energy summary saved to: {output_file}")
    return df

def main():
    """
    Main function to extract and organize energy data
    """
    print("Damietta Energy Data Extractor")
    print("=" * 40)
    
    # Configuration - update these paths if needed
    script_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    output_directory = os.path.join(script_root, 'output', 'mutated_prot_iteration_2')
    baseline_energy = -1.1096  # Your established baseline
    
    # Extract energy data
    print(f"Processing directory: {output_directory}")
    
    if os.path.exists(output_directory):
        results = extract_energies_from_directory(output_directory)
        
        if results:
            print(f"Found {len(results)} mutation results")
            
            # Create summary
            df = create_energy_summary(results, baseline_energy)
            
            # Display summary
            print("\n" + "=" * 60)
            print("ENERGY SUMMARY")
            print("=" * 60)
            print(df[['design_id', 'position_info', 'total_energy', 'energy_improvement']].to_string(index=False))
            
            # Identify best performers
            best_designs = df[df['energy_improvement'] > 0].sort_values('energy_improvement', ascending=False)
            
            if len(best_designs) > 0:
                print(f"\n🎉 SUCCESS: {len(best_designs)} designs improve upon baseline!")
                print(f"Best improvement: {best_designs.iloc[0]['energy_improvement']:.4f} kcal/mol")
                print(f"Best design: {best_designs.iloc[0]['design_id']} - {best_designs.iloc[0]['position_info']}")
            else:
                print("\n⚠️  No designs show energy improvement over baseline")
            
            print(f"\nNext steps:")
            print(f"1. Calculate RMSD values using PyMOL")
            print(f"2. Update energy_analysis.csv with RMSD data")
            print(f"3. Run Boltzmann analysis")
            
        else:
            print("No energy data found in the specified directory")
    else:
        print(f"Directory not found: {output_directory}")
        print("Please update the output_directory path in the script")

if __name__ == "__main__":
    main()
