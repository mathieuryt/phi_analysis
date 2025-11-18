import matplotlib.pyplot as plt

# Lecture des données
with open("/local/home/mr282803/Documents/irfu/AE/simu/new_AE/fichier_phigen.txt", "r") as f:
    phigen = [float(l.strip()) for l in f]

with open("/local/home/mr282803/Documents/irfu/AE/simu/new_AE/nb_phi.txt", "r") as f:
    phirec = [float(l.strip()) for l in f]

Acceptance = [phirec[i] / phigen[i] for i in range(len(phirec))]

# Sauvegarde de l'acceptance
with open('/local/home/mr282803/Documents/irfu/AE/simu/new_AE/fichier_AE_fit.txt', 'w', encoding='utf-8') as fichier:
    for element in Acceptance:
        fichier.write(str(element) + '\n')

# Bornes des bins
E_gamma = [2, 3.05333, 3.334, 3.566, 3.788, 4.006, 4.23333, 4.49, 4.76933, 5.1, 5.494, 6.076, 7.07933]  # 6 bornes → 5 bins

# Calcul des largeurs des bins
widths = [E_gamma[i+1] - E_gamma[i] for i in range(len(E_gamma) - 1)]

# Tracé de l'histogramme
plt.bar(E_gamma[:-1], Acceptance, width=widths, align='edge',
        edgecolor='black', facecolor='none', alpha=0.7)

plt.xlabel("Eγ [GeV]", fontsize=14)
plt.ylabel("Acceptance-efficiency", fontsize=14)
plt.ylim(0, 0.14)
plt.xlim(1, 8)  
plt.tick_params(axis='both', labelsize=12)
plt.grid(True, axis='y', linestyle='--', alpha=0.5)

plt.show()