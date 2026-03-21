# Load your structures
load prion_core_autopsf.pdb, baseline
load damietta_result0.pdb, design1
load damietta_result2.pdb, design2
load damietta_result4.pdb, design3
load damietta_result5.pdb, design4
load damietta_result6.pdb, design5
load damietta_result7.pdb, design6
load damietta_result8.pdb, design7
load damietta_result9.pdb, design8
load damietta_result10.pdb, design9


# Now align
align design1, baseline
align design2, baseline
align design3, baseline
align design4, baseline
align design5, baseline
align design6, baseline
align design7, baseline
align design8, baseline
align design9, baseline

# Create visualization
orient

# Save visuals
ray 1200, 900
png prion_overlay.png, dpi=300