import math
import numpy as np
import matplotlib.pyplot as plt


f2 = 1.28; # facteur car toute les simu n'ont pas reussi sur osg
Q = 1.2157*1e7
L_int = 1.51434*Q*f2 #nb-1
#L_int = 150843.9422 #pb-1

Br = 0.492 #channel e+ e-

wc = 1.0



with open("/Users/mr282803/Documents/analysisRGB/neutron/neutron_analysis/bornes_bins/bornes_t_1.txt", "r") as fichier:
    borne = [float(x) for ligne in fichier for x in ligne.split()]

with open("/Users/mr282803/Documents/analysisRGB/neutron/neutron_analysis/bornes_bins/bornes_t_moy_1.txt", "r") as fichier:
    centre_bins = [float(x) for ligne in fichier for x in ligne.split()]

with open("/Users/mr282803/Documents/analysisRGB/neutron/neutron_crosssection/N/integralN_t_1.txt", "r") as fichier:
    N_phi = [float(x) for ligne in fichier for x in ligne.split()]

with open("/Users/mr282803/Documents/analysisRGB/neutron/neutron_crosssection/N/errorN_t_1.txt", "r") as fichier:
    err_Nphi = [float(x) for ligne in fichier for x in ligne.split()]

with open("/Users/mr282803/Documents/analysisRGB/neutron/neutron_crosssection/A/integralNrec_mc_t_1.txt", "r") as fichier:
    A_rec = [float(x) for ligne in fichier for x in ligne.split()]

with open("/Users/mr282803/Documents/analysisRGB/neutron/neutron_crosssection/A/errorNrec_mc_t_1.txt", "r") as fichier:
    A_rec_err = [float(x) for ligne in fichier for x in ligne.split()]

with open("/Users/mr282803/Documents/analysisRGB/neutron/neutron_crosssection/A/integralNgen_mc_t_1.txt", "r") as fichier:
    A_gen = [float(x) for ligne in fichier for x in ligne.split()]

with open("/Users/mr282803/Documents/analysisRGB/neutron/neutron_crosssection/F/FluxQ2.txt", "r") as fichier:
    FluxQ2 = [float(x) for ligne in fichier for x in ligne.split()]


print("bornes des bins en t :", borne)
print("centre bins :", centre_bins)
print("nombre de phi :", N_phi)
print("erreur sur le nb phi :", err_Nphi)
print("MC : Nrec : ", A_rec)
print("MC : Ngen : ", A_gen)
print("Flux Q2 : ", FluxQ2)


Cs = []
err_tot = []

deltaQ2 = 1.41

for i in range(len(N_phi)):

    Cs.append(A_rec[i]/A_gen[i])  #nb



for i in range(len(N_phi)):

    err_tot.append((A_rec_err[i]/A_gen[i]))  #nb

print("Cross section", Cs)


# Création du plot
plt.figure(figsize=(8, 6))  # Taille de la figure
plt.errorbar(centre_bins, Cs, yerr=err_tot, fmt='o', color='#00b7eb', ecolor='#00b7eb',markersize=7, capsize=5, capthick=2, label = "1.30 < Q2 < 2.71 GeV2 ")

# Ajout des intervalles sous forme de ligne horizontale avec crochets plus grands
for i in range(len(borne) - 1):
    x_start, x_end = -borne[i], -borne[i+1]
    y_pos = 0  # Position sous l'axe des x pour afficher les crochets
    plt.hlines(y=y_pos, xmin=x_start, xmax=x_end, color='black', linewidth=2)  # Ligne horizontale
    plt.vlines(x=x_start, ymin=y_pos+0.0005*Cs[0], ymax=y_pos, color='black', linewidth=2)  # Crochet gauche plus grand
    plt.vlines(x=x_end, ymin=y_pos+0.0005*Cs[0], ymax=y_pos, color='black', linewidth=2)  # Crochet droit plus grand


# Personnalisation du graph
plt.xlabel("|t| [GeV^2]", fontsize=18)
plt.ylabel("Acceptance correction", size = 18)
plt.yscale('log')
plt.tick_params(axis='both', labelsize=14)
#plt.ylim(1e-8, 1e-1)
plt.ylim(1e-5, 1)  # Ajustement de l'échelle pour inclure les crochets
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()


plt.savefig("/Users/mr282803/Documents/analysisRGB/neutron/neutron_crosssection/acceptance_bin2.pdf", dpi = 500)



