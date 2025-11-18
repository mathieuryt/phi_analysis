import matplotlib.pyplot as plt
import numpy as np


with open("/local/home/mr282803/Documents/irfu/bins_Ephoton/global/fichier_bornes_bins.txt", "r") as fichier:
    borne = [float(ligne.strip()) for ligne in fichier]

with open("/local/home/mr282803/Documents/irfu/bins_Ephoton/global/fichier_E_moy.txt", "r") as fichier:
    centre_bins = [float(ligne.strip()) for ligne in fichier]

with open("/local/home/mr282803/Documents/irfu/bins_Ephoton/global/nb_phi.txt", "r") as fichier:
    N_phi = [float(ligne.strip()) for ligne in fichier]


with open("/local/home/mr282803/Documents/irfu/bins_Ephoton/global/err_nb_phi.txt", "r") as fichier:
    err_Nphi = [float(ligne.strip()) for ligne in fichier]


print("bornes des bins :", borne)
print("centre bins :", centre_bins)
print("nombre de phi :", N_phi)
print("erreur sur le nb phi :", err_Nphi)


# Création du plot
plt.figure(figsize=(8, 6))  # Taille de la figure
plt.errorbar(centre_bins, N_phi, yerr=err_Nphi, fmt='o', capsize=5, capthick=2)

# Ajout des intervalles sous forme de ligne horizontale avec crochets plus grands
for i in range(len(borne) - 1):
    x_start, x_end = borne[i], borne[i+1]
    y_pos = -50  # Position sous l'axe des x pour afficher les crochets
    plt.hlines(y=y_pos, xmin=x_start, xmax=x_end, color='black', linewidth=2)  # Ligne horizontale
    plt.vlines(x=x_start, ymin=y_pos-20, ymax=y_pos, color='black', linewidth=2)  # Crochet gauche plus grand
    plt.vlines(x=x_end, ymin=y_pos-20, ymax=y_pos, color='black', linewidth=2)  # Crochet droit plus grand

# Personnalisation du graph
plt.xlabel("E photon (GeV)", size = 18)
plt.ylabel("Number of φ", size = 18)
plt.title("2000 interesting event per bins of energy", size = 18)
plt.ylim(-100, 1000)  # Ajustement de l'échelle pour inclure les crochets
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()

# Affichage
plt.show()
