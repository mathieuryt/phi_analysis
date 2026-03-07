import math
import numpy as np
import matplotlib.pyplot as plt


Q = 1.2157*1e7 
L_int = 1.51434*Q  # nb-1
Br = 0.492

# Largeurs des bins Q2
deltaQ2 = [0.5, 5.0]   # bin 0 et bin 1
deltaxb = [0.81*0.29, 0.49*0.72]

# Labels pour la légende
labels = [
    "1.0 < Q² < 1.5 GeV² with <Q²> = 1.22",
    "1.5 < Q² < 7.0 GeV² with <Q²> = 2.72"
]

# Chargement du flux (commun)
with open("/Users/mr282803/Documents/analysisRGB/neutron/neutron_crosssection/F/FluxQ2.txt", "r") as fichier:
    FluxQ2 = [float(x) for ligne in fichier for x in ligne.split()]

plt.figure(figsize=(8, 6))

# -----------------------
# BOUCLE SUR LES BINS Q2
# -----------------------

for i in range(2):   # ici 2 bins Q2

    # Chargement des fichiers dépendant de i
    with open(f"/Users/mr282803/Documents/analysisRGB/neutron/neutron_analysis/bornes_bins/bornes_t_{i}.txt", "r") as fichier:
        borne = [float(x) for ligne in fichier for x in ligne.split()]

    with open(f"/Users/mr282803/Documents/analysisRGB/neutron/neutron_analysis/bornes_bins/bornes_t_moy_{i}.txt", "r") as fichier:
        centre_bins = [float(x) for ligne in fichier for x in ligne.split()]

    with open(f"/Users/mr282803/Documents/analysisRGB/neutron/neutron_crosssection/N/integralN_t_{i}.txt", "r") as fichier:
        N_phi = [float(x) for ligne in fichier for x in ligne.split()]

    with open(f"/Users/mr282803/Documents/analysisRGB/neutron/neutron_crosssection/N/errorN_t_{i}.txt", "r") as fichier:
        err_Nphi = [float(x) for ligne in fichier for x in ligne.split()]

    with open(f"/Users/mr282803/Documents/analysisRGB/neutron/neutron_crosssection/A/integralNrec_mc_t_{i}.txt", "r") as fichier:
        A_rec = [float(x) for ligne in fichier for x in ligne.split()]

    with open(f"/Users/mr282803/Documents/analysisRGB/neutron/neutron_crosssection/A/integralNgen_mc_t_{i}.txt", "r") as fichier:
        A_gen = [float(x) for ligne in fichier for x in ligne.split()]

    Cs = []
    err_tot = []

    for j in range(len(N_phi)):

        denom = (
            L_int
            * (A_rec[j] / A_gen[j])
            * Br
            * abs(borne[j+1] - borne[j])
            * FluxQ2[i]
            * deltaQ2[i]
            * deltaxb[i]
        )

        Cs.append(N_phi[j] / denom)
        err_tot.append(err_Nphi[j] / denom)

    # Plot pour ce bin Q2
    plt.errorbar(
        centre_bins,
        Cs,
        yerr=err_tot,
        fmt='o',
        markersize=7,
        capsize=5,
        capthick=2,
        label=labels[i]
    )

# -----------------------
# Personnalisation
# -----------------------


t_model1, cs_model1 = np.loadtxt(
    "/Users/mr282803/Documents/analysisRGB/neutron/neutron_crosssection/modelCs/cs_model1.txt",
    unpack=True
)

t_model2, cs_model2 = np.loadtxt(
    "/Users/mr282803/Documents/analysisRGB/neutron/neutron_crosssection/ModelCs/cs_model2.txt",
    unpack=True
)

# -----------------------
# Ajout au plot
# -----------------------

plt.plot(
    t_model1,
    cs_model1,
    linestyle='--',
    linewidth=2,
    color='blue',
    label="Model with <Q²> = 1.22 GeV² and <W> = 3.0 GeV²"
)

plt.plot(
    t_model2,
    cs_model2,
    linestyle='--',
    linewidth=2,
    color='orange',
    label="Model with <Q²> = 2.72 GeV² and <W> = 3.0 GeV²"
)

for i in range(len(borne) - 1):
    x_start, x_end = -borne[i], -borne[i+1] 
    y_pos = 0 # Position sous l'axe des x pour afficher les crochets
    plt.hlines(y=y_pos, xmin=x_start, xmax=x_end, color='black', linewidth=2) # Ligne horizontale
    plt.vlines(x=x_start, ymin=y_pos+2e-3, ymax=y_pos, color='black', linewidth=2) # Crochet gauche plus grand 
    plt.vlines(x=x_end, ymin=y_pos+2e-3, ymax=y_pos, color='black', linewidth=2) # Crochet droit plus grand

plt.xlabel("|t| [GeV²]", fontsize=18)
plt.ylabel(r"$\frac{d\sigma}{dt}$ reduced [nb/GeV$^4$]", size=18)
plt.yscale('log')
plt.tick_params(axis='both', labelsize=12)
plt.ylim(1e-3, 2.1*10e1)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()

plt.savefig("/Users/mr282803/Documents/analysisRGB/neutron/neutron_crosssection/test3.pdf", dpi=500)
plt.show()