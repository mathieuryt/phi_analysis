import math
import numpy as np
import matplotlib.pyplot as plt


L_int = 150843 #pb-1

Br = 2.954e-4 #channel e+ e-

wc = 1.0

nb_bins_gamma = 5

#BORNE 0 :


for j in range(nb_bins_gamma):


    with open(f"/local/home/mr282803/Documents/irfu/diff_cs_new/diff_bins/fichier_bornes_bins_{j}.txt", "r") as fichier:
        borne = [float(ligne.strip()) for ligne in fichier]

    with open(f"/local/home/mr282803/Documents/irfu/diff_cs_new/diff_bins/fichier_E_moy_{j}.txt", "r") as fichier:
        centre_bins = [float(ligne.strip()) for ligne in fichier]

    with open(f"/local/home/mr282803/Documents/irfu/diff_cs_new/fill_fit/nb_phi_bin_E_{j+1}.txt", "r") as fichier:
        N_phi = [float(ligne.strip()) for ligne in fichier]

    with open(f"/local/home/mr282803/Documents/irfu/diff_cs_new/fill_fit/error_phi_bin_E_{j+1}.txt", "r") as fichier:
        N_err = [float(ligne.strip()) for ligne in fichier]

    with open(f"/local/home/mr282803/Documents/irfu/diff_cs_new/flux/F.txt", "r") as fichier:
        F = [float(ligne.strip()) for ligne in fichier]

    with open(f"/local/home/mr282803/Documents/irfu/diff_cs_new/AE/AE_bin_E_{j+1}.txt", "r") as fichier:
        A = [float(ligne.strip()) for ligne in fichier]

    with open(f"/local/home/mr282803/Documents/irfu/diff_cs_new/diff_bins/fichier_bornes_bins.txt", "r") as fichier:
        borne_gamma = [float(ligne.strip()) for ligne in fichier]


    print("bornes des bins :", borne)
    print("centre bins :", centre_bins)
    print("nombre de phi :", N_phi)
    print("acceptance :", A)
    print("Flux de photon", F)
    print("bins en gamma", borne_gamma)


    Cs = []
    err_tot = []
    integral = 0


    for i in range(len(N_phi)):
        conv = 1e-3 #conversion en nb 
        diff_t = borne[i+1] - borne[i]
        Cs.append((N_phi[i]/(L_int*F[j]*A[i]*Br*wc*diff_t))*conv)  #nb
        integral += diff_t*Cs[i]



    for i in range(len(N_phi)):
        conv = 1e-3 #conversion en nb
        diff_t = borne[i+1] - borne[i]
        err_tot.append((N_err[i]/(L_int*F[j]*A[i]*Br*wc*diff_t))*conv)  #nb

    print("Cross section", Cs)
    print("erreur tot", err_tot)

    a = borne_gamma[j]
    b = borne_gamma[j+1]

    plt.figure(figsize=(8, 6))  # Taille de la figure
    plt.errorbar(centre_bins, Cs, yerr=err_tot, fmt='o', color='red', ecolor='red',markersize=7, capsize=5, capthick=2, label = f"E_gamma {a} - {b} GeV")

    # Ajout des intervalles sous forme de ligne horizontale avec crochets plus grands
    for i in range(len(borne) - 1):
        x_start, x_end = borne[i], borne[i+1]
        y_pos = 0  # Position sous l'axe des x pour afficher les crochets
        plt.hlines(y=y_pos, xmin=x_start, xmax=x_end, color='black', linewidth=2)  # Ligne horizontale
        plt.vlines(x=x_start, ymin=y_pos+10, ymax=y_pos, color='black', linewidth=2)  # Crochet gauche plus grand
        plt.vlines(x=x_end, ymin=y_pos+10, ymax=y_pos, color='black', linewidth=2)  # Crochet droit plus grand



    plt.xlabel("|t| [GeV²]", size = 18)
    plt.ylabel("Differential Cross section [nb/GeV²]", size = 18)
    #plt.yscale('log')
    plt.tick_params(axis='both', labelsize=14)
    #plt.ylim(1e-8, 1e-1)
    plt.ylim(0, 1700)  # Ajustement de l'échelle pour inclure les crochets
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(fontsize=16)


    plt.savefig(f"/local/home/mr282803/Documents/irfu/diff_cs_new/plot_dcs/bin{j}.pdf", dpi = 500)

    print(f"differential cross section integrée sur [{a} - {b} GeV] donne : ", integral, " nb")


    with open(f"/local/home/mr282803/Documents/irfu/form_factors/cs_bin{j}.txt", "w", encoding="utf-8") as f:
        for valeur in Cs:
            f.write(f"{valeur}\n")

    with open(f"/local/home/mr282803/Documents/irfu/form_factors/cs_error_bin{j}.txt", "w", encoding="utf-8") as f:
        for valeur in err_tot:
            f.write(f"{valeur}\n")

    with open(f"/local/home/mr282803/Documents/irfu/form_factors/bornes_bin{j}.txt", "w", encoding="utf-8") as f:
        for valeur in borne:
            f.write(f"{valeur}\n")

    with open(f"/local/home/mr282803/Documents/irfu/form_factors/centre_bin{j}.txt", "w", encoding="utf-8") as f:
        for valeur in centre_bins:
            f.write(f"{valeur}\n")

