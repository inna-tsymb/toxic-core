# Завантажуємо дикий тип як еталон
load prion_beta.pdb, WT

# Завантажуємо всі твої мутанти з папки beta
load dam4_beta.pdb, dam4
load des8_beta.pdb, des8
load des9_beta.pdb, des9
load des10_beta.pdb, des10
load des11_beta.pdb, des11
load des12_beta.pdb, des12
load n_res0_beta.pdb, n_res0
load n_res1_beta.pdb, n_res1
load n_res2_beta.pdb, n_res2
load n_res3_beta.pdb, n_res3
load n_res4_beta.pdb, n_res4
load res0_beta.pdb, res0
load res1_beta.pdb, res1
load res2_beta.pdb, res2
load res3_beta.pdb, res3
load res4_beta.pdb, res4
load result8_beta.pdb, result8
load trojan_beta.pdb, trojan

# Рахуємо RMSD для кожного відносно дикого типу і виводимо результати
print "====== RMSD RESULTS FOR BETA SHEETS ======"
align dam4, WT, cycles=0
align des8, WT, cycles=0
align des9, WT, cycles=0
align des10, WT, cycles=0
align des11, WT, cycles=0
align des12, WT, cycles=0
align n_res0, WT, cycles=0
align n_res1, WT, cycles=0
align n_res2, WT, cycles=0
align n_res3, WT, cycles=0
align n_res4, WT, cycles=0
align res0, WT, cycles=0
align res1, WT, cycles=0
align res2, WT, cycles=0
align res3, WT, cycles=0
align res4, WT, cycles=0
align result8, WT, cycles=0
align trojan, WT, cycles=0
print "=========================================="