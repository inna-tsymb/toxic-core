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
    candidate_root = os.path.expanduser('~/Desktop/KSE_Bioinformatics/Struct_Bio/pair_proj/toxic-core')
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
result_files = []
for folder in ['alpha_mutations', 'beta_mutations']:
    search_dir = os.path.join(output_root, folder)
    result_files.extend(glob.glob(os.path.join(search_dir, '**', 'dam_scored.pdb'), recursive=True))

result_files = [f for f in sorted(result_files) if os.path.isfile(f)]

print('=== RMSD CALCULATIONS ===')
print('Design | RMSD (Angstroms) | Filename')
print('-' * 50)

rmsd_results = []

for i, pdb_file in enumerate(result_files):
    design_name = f"design_{i+1}"
    filename = os.path.basename(pdb_file)

    cmd.load(pdb_file, design_name)
    cmd.align(f"{design_name} and name CA", "baseline and name CA")
    rmsd_value = cmd.rms_cur('baseline', design_name)

    rmsd_results.append({
        'design': design_name,
        'filename': filename,
        'rmsd': rmsd_value
    })

    print(f"{design_name:<7} | {rmsd_value:<15.3f} | {filename}")

python end

# Create visualization showing best designs
color green, baseline
color blue, design_1
color red, design_2
color yellow, design_3

# Show cartoon representation
show cartoon, all
hide lines, all

# Create publication-quality image
orient
set ray_shadows, off
bg white
ray 1200, 900
png rmsd_overlay_comparison.png, dpi=300

print "\n=== RMSD ANALYSIS COMPLETE ==="
print "Results saved in rmsd_overlay_comparison.png"
print "Copy the RMSD values above into your energy_analysis.csv file"
