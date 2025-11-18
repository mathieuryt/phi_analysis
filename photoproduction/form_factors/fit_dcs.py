import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Données

nb_bins_gamma = 5

ms = []
err_ms = []

with open("/local/home/mr282803/Documents/irfu/diff_cs_new/diff_bins/fichier_bornes_bins.txt", "r") as fichier:
    borne_gamma = [float(ligne.strip()) for ligne in fichier]

for j in range(nb_bins_gamma):

    with open(f"/local/home/mr282803/Documents/irfu/form_factors/bornes_bin{j}.txt", "r") as fichier:
        borne = [float(ligne.strip()) for ligne in fichier]

    with open(f"/local/home/mr282803/Documents/irfu/form_factors/centre_bin{j}.txt", "r") as fichier:
        centre_bins = [float(ligne.strip()) for ligne in fichier]

    with open(f"/local/home/mr282803/Documents/irfu/form_factors/cs_bin{j}.txt", "r") as fichier:
        Cs = [float(ligne.strip()) for ligne in fichier]

    with open(f"/local/home/mr282803/Documents/irfu/form_factors/cs_error_bin{j}.txt", "r") as fichier:
        err_tot = [float(ligne.strip()) for ligne in fichier]

    
    E_a = borne_gamma[j]
    E_b = borne_gamma[j+1]

    # Définition de la fonction à fitter
    def fit_func(t, ps, ms):
        return ps * (1 + t / (ms**2))**(-4)

    # Fit avec incertitudes (err_tot)
    popt, pcov = curve_fit(fit_func, centre_bins, Cs, sigma=err_tot, absolute_sigma=True, p0=[1000, 0.7])
    perr = np.sqrt(np.diag(pcov))  # Erreurs sur les paramètres


    ps_fit, ms_fit = popt


    ps_err, ms_err = perr

    print(f"Résultats du fit : ps = {ps_fit:.2f}, ms = {ms_fit:.2f} +- {ms_err:.4f} GeV")


    # Calcul de la courbe de fit
    t_fit = np.linspace(0.1, 1.2, 300)
    Cs_fit = fit_func(t_fit, *popt)


    # Plot avec les données et le fit
    plt.figure(figsize=(8, 6))
    plt.errorbar(centre_bins, Cs, yerr=err_tot, fmt='o', color='red', ecolor='red', markersize=7, capsize=5, capthick=2, label=f"data Eγ in [{E_a} - {E_b} GeV] ")
    plt.plot(t_fit, Cs_fit, label=f"Fit: $ps (1 - t/ms^2)^{{-4}}$\nps = {ps_fit:.1f}, ms = {ms_fit:.2f} GeV", color='blue', linewidth=2)

    # Ajout des crochets pour les bornes
    for i in range(len(borne) - 1):
        x_start, x_end = borne[i], borne[i+1]
        y_pos = 0
        plt.hlines(y=y_pos, xmin=x_start, xmax=x_end, color='black', linewidth=2)
        plt.vlines(x=x_start, ymin=y_pos+10, ymax=y_pos, color='black', linewidth=2)
        plt.vlines(x=x_end, ymin=y_pos+10, ymax=y_pos, color='black', linewidth=2)

    plt.xlabel("|t| [GeV²]", size=18)
    plt.ylabel("Differential Cross section [nb/GeV²]", size=18)
    plt.tick_params(axis='both', labelsize=14)
    plt.ylim(0, 1500)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(fontsize=14)
    plt.show()

    ms.append(ms_fit)
    err_ms.append(ms_err)


with open("/local/home/mr282803/Documents/irfu/form_factors/MS.txt", "w", encoding="utf-8") as f:
    for valeur in ms:
        f.write(f"{valeur}\n")

with open("/local/home/mr282803/Documents/irfu/form_factors/MS_err.txt", "w", encoding="utf-8") as f:
    for valeur in err_ms:
        f.write(f"{valeur}\n")








