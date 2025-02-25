from Algorithms.angle_estimation_algorithms import MBCEAlgorithm
from Algorithms.localization_algorithms import localization_algorithms
import auxiliary_functions.error_calculation as err_calc
from tabulate import tabulate
import numpy as np

#*****************************************STIMA DEGLI ANGOLI*****************************************
L = 1     #Numero di percorsi
N = 16    #Numero di antenne RIS
M = 4     #Numero di antenne BS
K = 4     #Numero di antenne UE
D = 15    #Numero di configurazioni di RIS

Theta_F = 1.1071
Theta_G = 0.9273
phi_F = 0.4036
phi_G = 0.6435
power = 0

#Stima angoli
mbce = MBCEAlgorithm(L, N, M, K, D)
Theta_G_hat, phi_F_hat, phi_G_hat, Theta_F_hat = mbce.run(Theta_F, Theta_G, phi_F, phi_G, power)


# Stampa i risultati
print("Valori Ottenuti")
print(tabulate([["Theta_G", Theta_G, Theta_G_hat],
                ["phi_F", phi_F, phi_F_hat],
                ["phi_G", phi_G, phi_G_hat],
                ["Theta_F", Theta_F, Theta_F_hat]],
                headers=[" ", "Valore Reale", "Valore Stimato"], tablefmt="pretty"))

print("Errori Percentuali")
print(tabulate([["Theta_G", err_calc.percentage_error(Theta_G, Theta_G_hat)],
                ["phi_F", err_calc.percentage_error(phi_F, phi_F_hat)],
                ["phi_G", err_calc.percentage_error(phi_G, phi_G_hat)],
                ["Theta_F", err_calc.percentage_error(Theta_F, Theta_F_hat)]],
                headers=["Angolo", "Errore Percentuale"], tablefmt="pretty"))



#*****************************************STIMA DELLA POSIZIONE*****************************************
c = 299792458           # Velocità della luce
q = np.array([0, 0])    # Posizione BS
r = np.array([30, 40])  # Posizione RIS
p = np.array([70,20])   # Posizione UE

tau_0 = 1.49174e-7      #Distanza RIS-UE nel caso LOS

RIS_Rotation = 3/2*np.pi



#Stima della posizione
loc = localization_algorithms(c)
#p_hat, alpha = loc.los_path_estimation(r, Theta_F_hat, phi_F_hat, tau_0)
p_hat, alpha = loc.los_path_estimation(r, RIS_Rotation + Theta_F, phi_F, tau_0)

print("Valori Ottenuti")
print(tabulate([["Posizione UE (X)", p[0], p_hat[0]],
                ["Posizione UE (Y)", p[1], p_hat[1]],
                ["Angolo alpha", "", alpha]],
                headers=["Valore Reale", "Valore Stimato"], tablefmt="pretty"))
print("Errori Percentuali")
print(tabulate([["Posizione UE (X)", err_calc.percentage_error(p[0], p_hat[0])],
                ["Posizione UE (Y)", err_calc.percentage_error(p[1], p_hat[1])]],
                headers=["Componente", "Errore Percentuale"], tablefmt="pretty"))












































def rotate_rectangle(center, width, height, angle):
    """Rotate a rectangle around its center by a given angle (in radians)."""
    # Rectangle vertices relative to its center
    rect = np.array([
        [-width / 2, -height / 2],  # Bottom-left
        [width / 2, -height / 2],  # Bottom-right
        [width / 2, height / 2],   # Top-right
        [-width / 2, height / 2]   # Top-left
    ])
    
    # Rotation matrix
    rotation_matrix = np.array([
        [np.cos(angle), -np.sin(angle)],
        [np.sin(angle), np.cos(angle)]
    ])
    
    # Rotate and translate back to the original center
    rotated_rect = np.dot(rect, rotation_matrix.T) + center
    return rotated_rect

#*****************************************************COPISCI COME FUNZIONANO GLI ANGOLI**********************************

def plot_angles(ax, point, angles, length=10, color='black', label='Angle'):
    """Plot AoD or AoA as arrows."""
    for angle in angles:
        end_point = point + length * np.array([np.cos(angle), np.sin(angle)])
        ax.annotate(
            '', xy=end_point, xytext=point,
            arrowprops=dict(facecolor=color, arrowstyle='->', lw=1.5)
        )
        ax.text(
            *end_point, f'{np.degrees(angle):.1f}°', color=color, fontsize=8,
            ha='center', va='center'
        )

def plot_graph(p, alpha, s_1x = 0, s_1y = 0):
    print(f"Posizione UE: {p}")
    print(f"Rotazione UE: {np.degrees(alpha)}")

    # --- Plot Setup ---
    fig, ax = plt.subplots(figsize=(8, 8))

    # Add rectangles for BS, RIS, and UE
    bs_rect = rotate_rectangle(q, 4, 2, 0)  # No rotation for BS
    ris_rect = rotate_rectangle(r, 4, 2, 0)  # No rotation for RIS
    ue_rect = rotate_rectangle(p, 4, 2, alpha)  # Rotate UE by alpha

    # Plot rectangles
    bs_polygon = Polygon(bs_rect, closed=True, color='blue', alpha=0.5, label="BS")
    ris_polygon = Polygon(ris_rect, closed=True, color='green', alpha=0.5, label="RIS")
    ue_polygon = Polygon(ue_rect, closed=True, color='red', alpha=0.5, label="UE")
    ax.add_patch(bs_polygon)
    ax.add_patch(ris_polygon)
    ax.add_patch(ue_polygon)

    # Plot BS, RIS, UE Points
    ax.scatter(*q, color='blue')
    ax.scatter(*r, color='green')
    ax.scatter(*p, color='red')

    if s_1x != 0:
        s = np.array([s_1x, s_1y])
        print(f"Posizione scatter: [{s_1x, s_1y}]")
        scatter_rect = rotate_rectangle(s, 2, 2, 0)
        scatter_polygon = Polygon(scatter_rect, closed=True, color='yellow', alpha=0.5, label='Scatter')
        ax.add_patch(scatter_polygon)
        ax.scatter(*s, color='yellow')
        ax.plot([r[0], s[0]], [r[1], s[1]], linestyle='--', color='green', label="RIS to Scatter")
        ax.plot([s[0], p[0]], [s[1], p[1]], linestyle='--', color='purple', label="Scatter to UE")
        ax.plot([r[0], p[0]], [r[1], p[1]], linestyle='--', color='orange', label="RIS to UE")
        
        # Plot AoD and AoA
        plot_angles(ax, s, [angles[1][0]], color='blue', label='AoD (Scatter)')
        plot_angles(ax, s, [angles[1][1]], color='green', label='AoA (Scatter)')
        
        ax.set_title("Canale di Comunicazione con RIS - Caso NLOS")
    else:    
        ax.set_title("Canale di Comunicazione con RIS - Caso LOS")
        # Connect points with lines
        ax.plot([r[0], p[0]], [r[1], p[1]], linestyle='--', color='green', label="RIS to UE")
        
        # Plot AoD and AoA
        plot_angles(ax, r, [angles[0][0], angles[0][1]], color='blue', label='AoD/AoA (RIS)')

    # Connect points with lines
    ax.plot([q[0], r[0]], [q[1], r[1]], linestyle='--', color='blue', label="BS to RIS")
    
    # Graph details
    ax.set_xlabel("Coordinate X")
    ax.set_ylabel("Coordinate Y")
    ax.axhline(0, color='black', linewidth=0.5)
    ax.axvline(0, color='black', linewidth=0.5)
    ax.grid(color='gray', linestyle='--', linewidth=0.5)
    ax.legend()
    plt.show()