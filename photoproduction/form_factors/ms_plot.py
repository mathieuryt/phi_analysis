import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Données



with open("/local/home/mr282803/Documents/irfu/diff_cs_new/diff_bins/fichier_bornes_bins.txt", "r") as fichier:
    borne = [float(ligne.strip()) for ligne in fichier]

with open("/local/home/mr282803/Documents/irfu/diff_cs_new/diff_bins/fichier_E_moy.txt", "r") as fichier:
    centre_bins = [float(ligne.strip()) for ligne in fichier]

with open("/local/home/mr282803/Documents/irfu/form_factors/MS.txt", "r") as fichier:
    ms = [float(ligne.strip()) for ligne in fichier]

with open("/local/home/mr282803/Documents/irfu/form_factors/MS_err.txt", "r") as fichier:
    ms_err = [float(ligne.strip()) for ligne in fichier]

ms_rich = [0.85]
ms_err_rich = [0.05]
centre_bins_rich = [3.2]

bornes_glueX = [8.20, 9.28, 10.36, 11.44]
centre_glueX = [8.93, 9.86, 10.82]
msX = [1.105, 1.472, 1.313]
err_msX = [0.168, 0.075, 0.049]


plt.figure(figsize=(8, 6))

plt.errorbar(centre_bins, ms, yerr=ms_err, fmt='o', color='blue', ecolor='blue', markersize=8, capsize=5, capthick=2, label=f"This work (CLAS12): ϕ -> e⁺e⁻")
plt.errorbar(centre_bins_rich, ms_rich, yerr = ms_err_rich, fmt='o', color='green', ecolor = 'green', markersize=8, capsize=5, capthick=2, label=f"R. Tyson (CLAS12) in Eγ [2.5 - 4.0 GeV]: ϕ -> e⁺e⁻")
plt.errorbar(centre_glueX, msX, yerr = err_msX, fmt='o', color='red', ecolor = 'red', markersize=8, capsize=5, capthick=2, label=f"GlueX : J/ψ -> e⁺e⁻")

# Ajout des crochets pour les bornes
for i in range(len(borne) - 1):
    x_start, x_end = borne[i], borne[i+1]
    y_pos = 0
    plt.hlines(y=y_pos, xmin=x_start, xmax=x_end, color='blue', linewidth=2)
    plt.vlines(x=x_start, ymin=y_pos+0.1, ymax=y_pos, color='blue', linewidth=2)
    plt.vlines(x=x_end, ymin=y_pos+0.1, ymax=y_pos, color='blue', linewidth=2)

for i in range(len(bornes_glueX) - 1):
    x_start, x_end = bornes_glueX[i], bornes_glueX[i+1]
    y_pos = 0
    plt.hlines(y=y_pos, xmin=x_start, xmax=x_end, color='red', linewidth=2)
    plt.vlines(x=x_start, ymin=y_pos+0.1, ymax=y_pos, color='red', linewidth=2)
    plt.vlines(x=x_end, ymin=y_pos+0.1, ymax=y_pos, color='red', linewidth=2)

plt.text(
    0.5, 0.5, "Preliminary",             # Position au centre
    transform=plt.gca().transAxes,       # Coordonnées relatives (0 à 1)
    fontsize=50,
    color='gray',
    alpha=0.2,                           # Transparence (0 = invisible, 1 = opaque)
    rotation=-45,                         # Angle en degrés (ex : 45 pour diagonale nette)
    ha='center', va='center'             # Alignement horizontal/vertical
)


plt.xlabel("Eγ [GeV]", size=18)
plt.ylabel("ms [GeV]", size=18)
plt.tick_params(axis='both', labelsize=14)
plt.ylim(0, 1.75)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=14)
plt.show()






