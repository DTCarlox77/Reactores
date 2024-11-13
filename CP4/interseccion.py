import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Constantes y valores calculados previamente
Cp_total = 858  # en J/mol K
Cp_extra = 100  # en J/h mol
H_rxn = 84600  # en J/mol
Delta_Cp = 16  # en J/mol K
k0 = 296105  # constante pre-exponencial en h⁻¹
Ea = 32500  # Energía de activación en J/mol
R = 8.314  # Constante de gas en J/mol K

# Definir las ecuaciones para BM y BE en función de T
def BM(T):
    T_ref1 = 398.15  # Temperatura de referencia en K para BM
    T_ref2 = 293.15  # Otra referencia de temperatura en K para BE
    return (Cp_total * (T - T_ref1) + Cp_extra * (T - T_ref2)) / (H_rxn + Delta_Cp * (T - 298.15))

def BE(T):
    u_A = 0.1266
    return (u_A * k0 * np.exp(-Ea / (R * T))) / (1 + u_A * k0 * np.exp(-Ea / (R * T)))

# Crear un rango de temperatura desde 398 K hasta que BE alcance 1
T_values = np.linspace(398, 550, 200)  # 200 puntos entre 398 y 550 K
BM_values = BM(T_values)
BE_values = BE(T_values)

# Encontrar el punto de intersección
def find_intersection(T_guess):
    return BM(T_guess) - BE(T_guess)

# Estimación inicial de la temperatura de intersección
T_intersect_guess = 450  # Aproximadamente en el centro del rango

# Resolver para el punto de intersección
T_intersect = fsolve(find_intersection, T_intersect_guess)[0]
BM_intersect = BM(T_intersect)
BE_intersect = BE(T_intersect)

# Graficar las curvas y el punto de intersección
plt.plot(T_values, BM_values, label='BM vs T', color='blue')
plt.plot(T_values, BE_values, label='BE vs T', color='red')
plt.plot(T_intersect, BM_intersect, 'go', label=f'Intersección: T={T_intersect:.2f}, BM=BE={BM_intersect:.4f}')

plt.xlabel('Temperatura (K)')
plt.ylabel('Valores de BM y BE')
plt.legend()
plt.title('Intersección de T vs BM y T vs BE')
plt.grid(True)
plt.show()

# Imprimir el punto de intersección
print(f'El punto de intersección ocurre a T = {T_intersect:.2f} K, donde BM = BE = {BM_intersect:.4f}')

# Elaborado por: Carlos Adrián Espinosa Luna | 4T1 Ing. Química