import glob
import os
import csv

# Attempt to resolve the project root
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
    raise RuntimeError('Could not determine project_root. Please set the project path explicitly.')

output_root = os.path.join(project_root, 'output', 'mutated_prot_iteration_2')

# ==========================================
# 1. ЗАВАНТАЖУЄМО ОБИДВІ БАЗИ (ALPHA ТА BETA)
# ==========================================
alpha_baseline_dir = os.path.join(output_root, 'prion_core_autopsf_openMM.pdb')
alpha_files = glob.glob(os.path.join(alpha_baseline_dir, '*.pdb'))
if alpha_files:
    cmd.load(alpha_files[0], 'baseline_alpha')
else:
    print(f"WARNING: Alpha baseline not found in {alpha_baseline_dir}")

beta_baseline_dir = os.path.join(output_root, 'prion_beta_autopsf_openMM.pdb')
beta_files = glob.glob(os.path.join(beta_baseline_dir, '*.pdb'))
if beta_files:
    cmd.load(beta_files[0], 'baseline_beta')
else:
    print(f"WARNING: Beta baseline not found in {beta_baseline_dir}")

# Знаходимо всі мутантні PDB
result_files = {'alpha_mutations': [], 'beta_mutations': []}
for folder in ['alpha_mutations', 'beta_mutations']:
    search_dir = os.path.join(output_root, folder)
    files = [f for f in sorted(glob.glob(os.path.join(search_dir, '**', '*.pdb'), recursive=True)) if os.path.isfile(f)]
    result_files[folder] = files

# Знаходимо файли дикого типу
wt_dirs = ['prion_core_autopsf_openMM.pdb', 'prion_beta_autopsf_openMM.pdb']
wt_files = []
for wt_dir in wt_dirs:
    wt_path = os.path.join(output_root, wt_dir)
    files = [f for f in glob.glob(os.path.join(wt_path, '*.pdb')) if os.path.isfile(f)]
    if files:
        wt_files.append((wt_dir, files[0]))

csv_path = os.path.join(project_root, 'energy_analysis.csv')
with open(csv_path) as f:
    csv_rows = list(csv.DictReader(f))

print('=== RMSD CALCULATIONS (CORRECTED: ALPHA->ALPHA, BETA->BETA) ===')
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

# ==========================================
# 2. РОЗРАХУНОК ДЛЯ МУТАНТІВ ІЗ ПРАВИЛЬНОЮ БАЗОЮ
# ==========================================
for folder, files in result_files.items():
    # Визначаємо, до якої бази рівняти цю папку!
    current_baseline = 'baseline_alpha' if folder == 'alpha_mutations' else 'baseline_beta'
    
    for pdb_file in files:
        design_name = f"design_{global_index}"
        filename = os.path.basename(pdb_file)
        result_dir = os.path.basename(os.path.dirname(pdb_file))
        design_id = f"{folder}/{result_dir}"

        matched_row = next((row for row in csv_rows if row['design_id'] == design_id), None)

        try:
            cmd.load(pdb_file, design_name)
            # Рівняємо до ПРАВИЛЬНОЇ бази
            cmd.align(f"{design_name} and name CA", f"{current_baseline} and name CA")
            rmsd_value = cmd.rms_cur(current_baseline, design_name)

            rmsd_results[folder].append({
                'design': design_name,
                'design_id': design_id,
                'filename': filename,
                'rmsd': rmsd_value
            })

            if matched_row:
                matched_row['rmsd'] = f"{rmsd_value:.3f}"
            write_csv()

            print(f"{design_name:<7} | {rmsd_value:<15.3f} | {filename:<20} | {folder} (vs {current_baseline})")
        except Exception as e:
            print(f"ERROR   | {design_id}: {e}")

        global_index += 1

# ==========================================
# 3. РОЗРАХУНОК ДЛЯ ДИКОГО ТИПУ
# ==========================================
for design_id, pdb_file in wt_files:
    matched_row = next((row for row in csv_rows if row['design_id'] == design_id), None)
    
    # Визначаємо базу для дикого типу
    current_baseline = 'baseline_alpha' if 'core' in design_id else 'baseline_beta'
    design_name = f"wt_{design_id.replace('.pdb','').replace('_','')}"
    
    try:
        cmd.load(pdb_file, design_name)
        cmd.align(f"{design_name} and name CA", f"{current_baseline} and name CA")
        rmsd_value = cmd.rms_cur(current_baseline, design_name)

        if matched_row:
            matched_row['rmsd'] = f"{rmsd_value:.3f}"
        write_csv()

        print(f"{design_name:<7} | {rmsd_value:<15.3f} | {os.path.basename(pdb_file):<20} | {design_id} (vs {current_baseline})")
    except Exception as e:
        print(f"ERROR   | {design_id}: {e}")

print(f"\nRMSD values written to {csv_path}")

# ==========================================
# 4. РЕНДЕР КАРТИНОК ІЗ ПРАВИЛЬНИМИ БАЗАМИ
# ==========================================
for folder, results in rmsd_results.items():
    if not results:
        continue
    
    current_baseline = 'baseline_alpha' if folder == 'alpha_mutations' else 'baseline_beta'
    
    cmd.hide('everything', 'all')
    cmd.show('cartoon', current_baseline)
    cmd.color('green', current_baseline)
    
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

print("\n=== RMSD ANALYSIS COMPLETE ===")
print("Separate overlay images saved with correct respective baselines.")
