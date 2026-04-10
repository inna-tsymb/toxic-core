# PyMOL Script for Batch RMSD Calculations
# Enhanced version for your 10 mutation results

python
import glob
import os

# Attempt to resolve the project root.
# In PyMOL, __file__ may not be defined, so fall back to an explicit user project path.
project_root = None
try:
    script_dir = os.path.abspath(os.path.dirname(__file__))
    candidate_root = os.path.abspath(os.path.join(script_dir, '..'))
    if os.path.isdir(os.path.join(candidate_root, 'output')):
        project_root = candidate_root
except NameError:
    project_root = None

if project_root is None:
    candidate_root = os.path.expanduser('~/AnnaProjects/toxic-core')
    if os.path.isdir(os.path.join(candidate_root, 'output')):
        project_root = candidate_root

if project_root is None:
    raise RuntimeError('Could not determine project_root. Please set the project path explicitly in the script.')

output_root = os.path.join(project_root, 'output', 'mutated_prot_iteration_2')

# Load baseline structure from the prion_core_autopsf_openMM.pdb directory
baseline_dir = os.path.join(output_root, 'prion_core_autopsf_openMM.pdb')
baseline_files = glob.glob(os.path.join(baseline_dir, '*.pdb'))
baseline_file = baseline_files[0] if baseline_files else None

if baseline_file and os.path.exists(baseline_file):
    cmd.load(baseline_file, 'baseline')
else:
    raise RuntimeError(f'Baseline PDB not found in {baseline_dir}')

# Find all final mutant PDB files under alpha_mutations and beta_mutations
result_files = {'alpha_mutations': [], 'beta_mutations': []}
for folder in ['alpha_mutations', 'beta_mutations']:
    search_dir = os.path.join(output_root, folder)
    files = [f for f in sorted(glob.glob(os.path.join(search_dir, '**', '*.pdb'), recursive=True)) if os.path.isfile(f)]
    result_files[folder] = files

# Also include WT baseline directories (design_id = bare directory name)
wt_dirs = ['prion_core_autopsf_openMM.pdb', 'prion_beta_autopsf_openMM.pdb']
wt_files = []
for wt_dir in wt_dirs:
    wt_path = os.path.join(output_root, wt_dir)
    files = [f for f in glob.glob(os.path.join(wt_path, '*.pdb')) if os.path.isfile(f)]
    if files:
        wt_files.append((wt_dir, files[0]))

import csv

csv_path = os.path.join(project_root, 'energy_analysis.csv')
with open(csv_path) as f:
    csv_rows = list(csv.DictReader(f))

print('=== RMSD CALCULATIONS ===')
print('Design | RMSD (Angstroms) | Filename | Group')
print('-' * 60)

rmsd_results = {'alpha_mutations': [], 'beta_mutations': []}

def write_csv():
    fieldnames = csv_rows[0].keys()
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(csv_rows)

global_index = 1
for folder, files in result_files.items():
    for pdb_file in files:
        design_name = f"design_{global_index}"
        filename = os.path.basename(pdb_file)
        # Build design_id as used in the CSV (e.g. alpha_mutations/results_1)
        result_dir = os.path.basename(os.path.dirname(pdb_file))
        design_id = f"{folder}/{result_dir}"

        # Find matching CSV row
        matched_row = None
        for row in csv_rows:
            if row['design_id'] == design_id:
                matched_row = row
                break

        # Skip if RMSD already computed
        if matched_row and matched_row.get('rmsd', '').strip():
            print(f"SKIP    | already computed ({matched_row['rmsd']:>6} Å) | {filename:<20} | {folder}/{result_dir}")
            global_index += 1
            continue

        try:
            cmd.load(pdb_file, design_name)
            cmd.align(f"{design_name} and name CA", "baseline and name CA")
            rmsd_value = cmd.rms_cur('baseline', design_name)

            rmsd_results[folder].append({
                'design': design_name,
                'design_id': design_id,
                'filename': filename,
                'rmsd': rmsd_value
            })

            # Update CSV row and write immediately so partial results are saved
            if matched_row:
                matched_row['rmsd'] = f"{rmsd_value:.3f}"
            write_csv()

            print(f"{design_name:<7} | {rmsd_value:<15.3f} | {filename:<20} | {folder}")
        except Exception as e:
            print(f"ERROR   | {design_id}: {e}")

        global_index += 1

# Process WT directories separately
for design_id, pdb_file in wt_files:
    matched_row = None
    for row in csv_rows:
        if row['design_id'] == design_id:
            matched_row = row
            break

    if matched_row and matched_row.get('rmsd', '').strip():
        print(f"SKIP    | already computed ({matched_row['rmsd']:>6} Å) | {os.path.basename(pdb_file):<20} | {design_id}")
        continue

    design_name = f"wt_{design_id.replace('.pdb','').replace('_','')}"
    try:
        cmd.load(pdb_file, design_name)
        cmd.align(f"{design_name} and name CA", "baseline and name CA")
        rmsd_value = cmd.rms_cur('baseline', design_name)

        if matched_row:
            matched_row['rmsd'] = f"{rmsd_value:.3f}"
        write_csv()

        print(f"{design_name:<7} | {rmsd_value:<15.3f} | {os.path.basename(pdb_file):<20} | {design_id}")
    except Exception as e:
        print(f"ERROR   | {design_id}: {e}")

print(f"RMSD values written to {csv_path}")

# Save separate PNG per group
for folder, results in rmsd_results.items():
    if not results:
        continue
    cmd.hide('everything', 'all')
    cmd.show('cartoon', 'baseline')
    cmd.color('green', 'baseline')
    for r in results:
        cmd.show('cartoon', r['design'])
    if folder == 'alpha_mutations':
        cmd.color('blue', ' or '.join(r['design'] for r in results))
    else:
        cmd.color('red', ' or '.join(r['design'] for r in results))
    cmd.orient()
    cmd.set('ray_shadows', 'off')
    cmd.bg_color('white')
    cmd.ray(1200, 900)
    out_png = os.path.join(project_root, 'output', 'mutated_prot_iteration_2', 'analysis', f'rmsd_overlay_{folder}.png')
    cmd.png(out_png, dpi=300)
    print(f"Saved {out_png}")

python end

print "\n=== RMSD ANALYSIS COMPLETE ==="
print "Separate overlay images saved for alpha_mutations and beta_mutations"
