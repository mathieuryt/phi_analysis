import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad


me = 0.00051
Eb = 10.6  
Q2_max = 0.02
d = 5.0
X0 = 929.0
Mp = 0.9383

def N_EPA(Eb, Eg, Q2_max):
    alpha = 1.0 / 137.0
    PI = 3.14
    x = Eg / Eb
    me = 0.00051
    Q2_min = (me**2) * (x**2) / (1 - x)
    
    return (1 / Eb) * alpha / (PI * x) * ((1 - x + x**2 / 2) * math.log(Q2_max / Q2_min) - (1 - x))

def N_Brem(Eg, Eb, d, X0):
    return (0.5 * d / X0) * (1 / Eg) * ((4.0 / 3.0) - (4.0 / 3.0) * (Eg / Eb) + (Eg**2) / (Eb**2))

def Sum_dist(Eg):
    return (N_EPA(Eb, Eg, Q2_max) + N_Brem(Eg, Eb, d, X0))


E_gamma = np.linspace(1.51, 10.599999999, 100000)

Dist_EPA = []
Dist_Brem = []
Dist_sum = []

for i in range(len(E_gamma)):
    Dist_Brem.append(N_Brem(E_gamma[i], Eb, d, X0))

for i in range(len(E_gamma)):
    Dist_EPA.append(N_EPA(Eb, E_gamma[i], Q2_max))

for i in range(len(E_gamma)):
    Dist_sum.append(Sum_dist(E_gamma[i]))

filename = '/local/home/mr282803/Documents/irfu/flux_bin/fichier_bornes_bins.txt'

# Ouvre le fichier et lit les données dans une liste
with open(filename, 'r') as file:
    # Utilise list comprehension pour lire chaque ligne, convertir en float et les stocker dans une liste
    bornes = [float(line.strip()) for line in file]

# Affiche la liste des données
print("les bornes des bins sont : ")
print(bornes)

filename2 = '/local/home/mr282803/Documents/irfu/flux_bin/fichier_E_moy.txt'

# Ouvre le fichier et lit les données dans une liste
with open(filename2, 'r') as file:
    # Utilise list comprehension pour lire chaque ligne, convertir en float et les stocker dans une liste
    centre_bins = [float(line.strip()) for line in file]

# Affiche la liste des données
print("les centre des bins en energie sont :")
print(centre_bins)

Flux = []

for i in range(0, len(bornes) - 1):

    print(i)
    a = bornes[i]
    b = bornes[i+1]
    result, error = quad(Sum_dist, a, b)
    Flux.append(result)

print("l'integral sur chaque bins est :")
print(Flux)

plt.plot(E_gamma, Dist_EPA, label='EPA', color='blue')
plt.plot(E_gamma, Dist_Brem, label='Brem', color='green')
plt.plot(E_gamma, Dist_sum, label='Sum', color='red')
plt.xlabel('E gamma [Gev]')
plt.ylabel('Density')
plt.grid(True)
plt.legend()
plt.show()


# Création du plot
plt.figure(figsize=(8, 6))  # Taille de la figure
plt.plot(centre_bins, Flux, 'o')

# Ajout des intervalles sous forme de ligne horizontale avec crochets plus grands
for i in range(len(bornes) - 1):
    x_start, x_end = bornes[i], bornes[i+1]
    y_pos = 0  # Position sous l'axe des x pour afficher les crochets
    plt.hlines(y=y_pos, xmin=x_start, xmax=x_end, color='black', linewidth=2)  # Ligne horizontale
    plt.vlines(x=x_start, ymin=y_pos-0.0005, ymax=y_pos, color='black', linewidth=2)  # Crochet gauche plus grand
    plt.vlines(x=x_end, ymin=y_pos-0.0005, ymax=y_pos, color='black', linewidth=2)  # Crochet droit plus grand

# Personnalisation du graph
plt.xlabel("E_gamma [GeV]", size = 18)
plt.ylabel("Flux de photon", size = 18)
plt.ylim(-0.001, 0.02)

plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()

# Affichage
plt.show()


with open("/local/home/mr282803/Documents/irfu/flux_bin/F.txt", "w") as f:
    for valeur in Flux:
        f.write(f"{valeur}\n")





