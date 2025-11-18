import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Données



with open("/local/home/mr282803/Documents/irfu/diff_cs_new/diff_bins/fichier_bornes_bins.txt", "r") as fichier:
    borne = [float(ligne.strip()) for ligne in fichier]

with open("/local/home/mr282803/Documents/irfu/diff_cs_new/diff_bins/fichier_E_moy.txt", "r") as fichier:
    centre_bins = [float(ligne.strip()) for ligne in fichier]

with open("/local/home/mr282803/Documents/irfu/form_factors/B.txt", "r") as fichier:
    B = [float(ligne.strip()) for ligne in fichier]

with open("/local/home/mr282803/Documents/irfu/form_factors/B_err.txt", "r") as fichier:
    B_err = [float(ligne.strip()) for ligne in fichier]


centre_CLAS = [1.65, 1.85, 2.05, 2.25, 2.45, 2.65, 2.85, 3.05, 3.25, 3.45]

B_CLAS = [3.1636, 2.8647, 3.5221, 2.8456, 2.2048, 3.0617, 3.4009, 3.8577, 3.1422, 3.3680]

B_CLAS_err = [0.33, 0.15, 0.15, 0.12, 0.06, 0.34, 0.20, 0.35, 0.35, 0.47]

plt.figure(figsize=(8, 6))

plt.errorbar(centre_bins, B, yerr=B_err, fmt='o', color='blue', ecolor='blue', markersize=8, capsize=5, capthick=2, label=f"This work (CLAS12): ϕ -> e⁺e⁻")
plt.errorbar(centre_CLAS, B_CLAS, yerr=B_CLAS_err, fmt='o', color='red', ecolor='red', markersize=8, capsize=5, capthick=2, label=f"CLAS: ϕ -> KsKl")

# Ajout des crochets pour les bornes
for i in range(len(borne) - 1):
    x_start, x_end = borne[i], borne[i+1]
    y_pos = 0
    plt.hlines(y=y_pos, xmin=x_start, xmax=x_end, color='blue', linewidth=2)
    plt.vlines(x=x_start, ymin=y_pos+0.1, ymax=y_pos, color='blue', linewidth=2)
    plt.vlines(x=x_end, ymin=y_pos+0.1, ymax=y_pos, color='blue', linewidth=2)


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
plt.ylabel("B [GeV²]", size=18)
plt.tick_params(axis='both', labelsize=14)
plt.ylim(0, 6.5)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=14)
plt.show()






