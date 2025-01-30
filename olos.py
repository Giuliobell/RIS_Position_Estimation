import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from sympy import symbols, Eq, solve, tan
from scipy.optimize import least_squares

# Costante: Velocità della luce
c = 299792458

def calculate_p_and_distances(r, tau_k, theta_F_k, phi_F_k, alpha_trial):
    d_k_1 = c * tau_k  # Distanza totale
    cos_theta = np.cos(theta_F_k)
    sin_theta = np.sin(theta_F_k)
    cos_phi = np.cos(phi_F_k + alpha_trial)
    sin_phi = np.sin(phi_F_k + alpha_trial)

    # Calcolo di p
    p_x = r[0] + d_k_1 * cos_theta + (c * tau_k - d_k_1) * cos_phi
    p_y = r[1] + d_k_1 * sin_theta - (c * tau_k - d_k_1) * sin_phi
    p = np.array([p_x, p_y])

    return p, d_k_1, c * tau_k - d_k_1


def calculate_scatter_positions(r, p, theta_F_k, phi_F_k, alpha):
    sx, sy = symbols('sx sy')

    # Sistemi lineari derivati da (4.11)
    eq1 = Eq(np.tan(np.pi - (phi_F_k + alpha)), (p[1] - sy) / (p[0] - sx))
    eq2 = Eq(np.tan(theta_F_k), (sy - r[1]) / (sx - r[0]))
    
    solution = solve([eq1, eq2], (sx, sy), dict=True)
    
    if not solution:
        raise ValueError("No solution found for the intersection points.")
    
    print(solution)
    return np.array([float(solution[sx]), float(solution[sy])])


def residual_function(variables, q, scatter_positions, measurements):
    """
    Funzione di residuo da minimizzare con LMA.
    
    :param variables: Variabili da ottimizzare [px, py, alpha]
    :param q: Posizione della Base Station (BS)
    :param scatter_positions: Posizioni degli scatter (array di coordinate)
    :param measurements: Misurazioni osservate (array di [tau, theta_F, phi_F])
    :return: Residuo
    """
    px, py, alpha = variables
    residui_totali = []
    p = np.array([px, py])

    for (sx, sy), (tau_k, theta_F_k, phi_F_k) in zip(scatter_positions, measurements):
        d_bs = np.linalg.norm([sx - q[0], sy - q[1]])  # Distanza BS -> scatter
        d_ue = np.linalg.norm([px - sx, py - sy])      # Distanza scatter -> UE

        # Residui
        res_toa = d_bs + d_ue - c * tau_k
        res_aod = np.arctan2(sy - q[1], sx - q[0]) - phi_F_k
        res_aoa = np.arctan2(py - sy, px - sx) - (theta_F_k + alpha)

        residui_totali.extend([res_toa, res_aod, res_aoa])

    return residui_totali


def olos_estimation(q, r, alpha_range, measurements, resolution):

    alpha_trials = np.arange(alpha_range[0], alpha_range[1], resolution)
    best_cost = float('inf')
    best_solution = None

    for alpha_trial in alpha_trials:
        p, d1_1, d2_1 = calculate_p_and_distances(r, *measurements[0], alpha_trial)
        scatter_positions = [calculate_scatter_positions(r, p, theta_F, phi_F, alpha_trial)
        for _, theta_F, phi_F in measurements]

        # Ottimizzazione con LMA
        initial_guess = [*p, alpha_trial]
        result = least_squares(
            residual_function,
            initial_guess,
            args=(q, scatter_positions, measurements),
            method='lm'
        )

        if result.cost < best_cost:
            best_cost = result.cost
            best_solution = (result.x[:2], result.x[2], scatter_positions)

    return best_solution


# Parametri di esempio
#q = np.array([0, 0])  # Posizione della Base Station
#r = np.array([10, 10])  # Posizione RIS
#alpha_range = [-np.radians(30), np.radians(30)]  # Intervallo per alpha trial
#resolution = np.radians(1)  # Risoluzione di alpha trial
#measurements = [
#    (1.33e-7, np.radians(45), np.radians(135)),
#    (1.67e-7, np.radians(60), np.radians(120)),
#    (2.00e-7, np.radians(30), np.radians(150))
#]

# Calcolo
best_p, best_alpha, best_scatter_positions = olos_estimation(q, r, alpha_range, measurements, resolution)
print("Posizione UE stimata:", best_p)
print("Angolo alpha stimato (in gradi):", np.degrees(best_alpha))
print("Posizioni degli scatter:", best_scatter_positions)
