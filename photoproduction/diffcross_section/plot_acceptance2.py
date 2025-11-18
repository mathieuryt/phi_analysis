import matplotlib.pyplot as plt

# Lecture des données
with open("/local/home/mr282803/Documents/irfu/diff_cs_test/AE/AE_bin_E_1.txt", "r") as f:
    Acceptance = [float(l.strip()) for l in f]

with open("/local/home/mr282803/Documents/irfu/diff_cs_test/diff_bins/fichier_bornes_bins_0.txt", "r") as f:
    bin_t = [float(l.strip()) for l in f]



# Calcul des largeurs des bins
widths = [bin_t[i+1] - bin_t[i] for i in range(len(bin_t) - 1)]

# Tracé de l'histogramme
plt.bar(bin_t[:-1], Acceptance, width=widths, align='edge',
        edgecolor='black', facecolor='none', alpha=0.7)

plt.xlabel(" |t| [GeV²]", fontsize=14)
plt.ylabel("Acceptance-efficiency", fontsize=14)
plt.ylim(0, 0.08)
plt.xlim(0, 1.2)  
plt.tick_params(axis='both', labelsize=12)
plt.grid(True, axis='y', linestyle='--', alpha=0.5)

plt.show()