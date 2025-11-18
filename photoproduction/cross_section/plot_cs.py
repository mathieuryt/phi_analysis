import math
import numpy as np
import matplotlib.pyplot as plt


L_int = 150843.9422 #pb-1

Br = 2.954e-4 #channel e+ e-

wc = 1.0

wc_rich = 0.388




with open("/local/home/mr282803/Documents/irfu/cross_section/fichier_bornes_bins.txt", "r") as fichier:
    borne = [float(ligne.strip()) for ligne in fichier]

with open("/local/home/mr282803/Documents/irfu/cross_section/fichier_E_moy.txt", "r") as fichier:
    centre_bins = [float(ligne.strip()) for ligne in fichier]

with open("/local/home/mr282803/Documents/irfu/cross_section/nb_phi.txt", "r") as fichier:
    N_phi = [float(ligne.strip()) for ligne in fichier]

with open("/local/home/mr282803/Documents/irfu/cross_section/err_nb_phi.txt", "r") as fichier:
    err_Nphi = [float(ligne.strip()) for ligne in fichier]

with open("/local/home/mr282803/Documents/irfu/cross_section/F.txt", "r") as fichier:
    F = [float(ligne.strip()) for ligne in fichier]

with open("/local/home/mr282803/Documents/irfu/cross_section/A.txt", "r") as fichier:
    A = [float(ligne.strip()) for ligne in fichier]


print("bornes des bins :", borne)
print("centre bins :", centre_bins)
print("nombre de phi :", N_phi)
print("erreur sur le nb phi :", err_Nphi)
print("acceptance :", A)
print("Flux de photon", F)

Cs = []
err_tot = []

for i in range(len(N_phi)):
    conv = 1e-3 #conversion en nb 
    Cs.append((N_phi[i]/(L_int*F[i]*A[i]*Br*wc))*conv)  #nb



for i in range(len(N_phi)):

    conv = 1e-3 #conversion en nb
    err_tot.append((err_Nphi[i]/(L_int*F[i]*A[i]*Br*wc))*conv)  #nb

print("Cross section", Cs)


# Création du plot
plt.figure(figsize=(8, 6))  # Taille de la figure
plt.errorbar(centre_bins, Cs, yerr=err_tot, fmt='o', color='#00b7eb', ecolor='#00b7eb',markersize=7, capsize=5, capthick=2, label = "This work")

# Ajout des intervalles sous forme de ligne horizontale avec crochets plus grands
for i in range(len(borne) - 1):
    x_start, x_end = borne[i], borne[i+1]
    y_pos = 0  # Position sous l'axe des x pour afficher les crochets
    plt.hlines(y=y_pos, xmin=x_start, xmax=x_end, color='black', linewidth=2)  # Ligne horizontale
    plt.vlines(x=x_start, ymin=y_pos+10, ymax=y_pos, color='black', linewidth=2)  # Crochet gauche plus grand
    plt.vlines(x=x_end, ymin=y_pos+10, ymax=y_pos, color='black', linewidth=2)  # Crochet droit plus grand

cs_rich = [130.936, 237.01, 352.178, 383.771, 412.107, 452.538, 381.355, 382.25]
cs_rich2 = [cs_rich[i]*wc_rich for i in range(len(cs_rich))]

cs_err_rich = [26.5443, 43.6486, 38.2942, 46.8208, 58.4148, 49.0882, 60.6068,
65.9573]
cs_err_rich2 = [cs_err_rich[i]*wc_rich for i in range(len(cs_err_rich))]

E_gamma_rich = [2.78624, 3.25028, 3.74083, 4.23346, 4.73797, 5.24201, 5.74143, 6.50149]

plt.errorbar(E_gamma_rich, cs_rich2,yerr = cs_err_rich2, fmt='o', color='green', ecolor='green',markersize=7, capsize=5, capthick=2,label = "CLAS12 previous analysis")

print(cs_rich)

# Personnalisation du graph
plt.xlabel("Eγ [GeV]", fontsize=18)
plt.ylabel("Cross section [nb]", size = 18)
#plt.yscale('log')
plt.tick_params(axis='both', labelsize=14)
#plt.ylim(1e-8, 1e-1)
plt.ylim(0, 250)  # Ajustement de l'échelle pour inclure les crochets
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()



plt.savefig("/local/home/mr282803/Documents/irfu/cross_section/test.pdf", dpi = 500)



