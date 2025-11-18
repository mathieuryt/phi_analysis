import math
import numpy as np
import matplotlib.pyplot as plt

# Fonction Lambda en Python
def Lambda(x, y, z):
    return (x - y - z) * (x - y - z) - 4 * y * z

# Fonction T_min en Python
def T_min(E):
    Mp = 0.938
    M_phi = 1.019
    ma_2 = 0
    mb_2 = Mp * Mp
    m1_2 = 1.019 * 1.019
    m2_2 = Mp * Mp
    s = Mp * Mp + 2 * Mp * E
    return abs(ma_2 + m1_2 - (1 / (2 * s)) * ((s + ma_2 - mb_2) * (s + m1_2 - m2_2) - math.sqrt(Lambda(s, ma_2, mb_2) * Lambda(s, m1_2, m2_2))))

# Fonction T_max en Python
def T_max(E):
    Mp = 0.938
    M_phi = 1.019
    ma_2 = 0
    mb_2 = Mp * Mp
    m1_2 = 1.019 * 1.019
    m2_2 = Mp * Mp
    s = Mp * Mp + 2 * Mp * E
    return abs(ma_2 + m1_2 - (1 / (2 * s)) * ((s + ma_2 - mb_2) * (s + m1_2 - m2_2) + math.sqrt(Lambda(s, ma_2, mb_2) * Lambda(s, m1_2, m2_2))))

# Définir l'intervalle d'E entre 0 et 11
E_values = np.linspace(1.8, 11, 500)

# Calcul des valeurs de T_min et T_max pour chaque E
T_min_values = [T_min(E) for E in E_values]
T_max_values = [T_max(E) for E in E_values]

# Création du plot
plt.figure(figsize=(10, 6))

# Tracer T_min et T_max
plt.plot(E_values, T_min_values, label='$T_{min}(E)$', color='red')
plt.plot(E_values, T_max_values, label='$T_{max}(E)$', color='blue')

# Ajouter des titres et des labels
plt.title('$T_{min}(E)$ et $T_{max}(E)$ en fonction de $E$', fontsize=14)
plt.xlabel('E [GeV]', fontsize=12)
plt.ylabel('T [GeV²]', fontsize=12)

# Ajouter une légende
plt.legend()

# Afficher le plot
plt.grid(True)
plt.show()
