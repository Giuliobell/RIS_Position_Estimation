'''
NOTE TESI
- usatolo solo 1 scatter per NLOS ma estendibile (numero minimo)
- applicata parte algoritmo NLOS per calcolare la posizione degli scatter 
- sistema nome d k-esimi nella formula 4.13
- aggiunti offset per le posizioni degli angoli finali a causa della geometria del problema
- verifica il sistema di riferimento angoli negli angoli d'arrivo

- Bisogna usare matlab for python per applicare l'algoritmo rootmusic
- Assumo D = lambda_c/2 per avere una semplificazione e per riprodurre fedelmente una giusta scelta
- Gli angoli li do io in inout all'algoritmo perchè in realtà sarebbero quelli che l'antenna misura
- Uso canale senza rumore gaussiano bianco
- Caso particolare in cui la matrice Theta ovvero la DFT non ha dimensione D*D con gli zeri ma ha proprio dimensione N*N come la matrice di antenne RIS
'''
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from sympy import symbols, Eq, solve, tan

# --- Constants ---
c = 299792458  # Speed of light
q = np.array([0, 0])  # BS position
r = np.array([30, 40])  # RIS position
angles = {
    0 : [(3/2*np.pi) + 1.1071, 0.4036],  #theta, phi LOS path
    1 : [(3/2*np.pi) + 1.373400767, 0.7264380874],  #theta, phi path 1
}
TOA = (1.49174e-7, 1.558021216e-7)  #tau_o, tau_1, tau_2, tau_3


# --- Functions ---
def los_param(tau_0: float, angles):
    """Calculate LOS parameters."""
    p = np.add(r, c * tau_0 * np.array((np.cos(angles[0][0]), np.sin(angles[0][0]))).T)
    alpha = np.pi + angles[0][0] - angles[0][1]
    return p, alpha

def nlos_param(tau_0: float, angles):
    """Calculate Non-Line-of-Sight (NLOS) parameters."""
    p, alpha = los_param(tau_0, angles)

    alpha = np.pi - alpha

    # Define symbolic variables for intersection points
    s_1x, s_1y = symbols('s_1x s_1y')

    # Define equations using SymPy's tan function
    eq_1 = Eq(tan(np.pi - (angles[1][1] + alpha)), (p[1] - s_1y) / (p[0] - s_1x))
    eq_2 = Eq(tan(angles[1][0]), (s_1y - r[1]) / (s_1x - r[0]))

    # Solve equations for s_1x and s_1y
    solutions = solve([eq_1, eq_2], (s_1x, s_1y), dict=True)
    
    if not solutions:
        raise ValueError("No solution found for the intersection points.")
    
    solution = solutions[0]  # Select the first solution
    return p, alpha, float(solution[s_1x]), float(solution[s_1y])

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



p, alpha = los_param(TOA[0], angles)
plot_graph(p, np.pi - alpha)

p, alpha, s_1x, s_1y = nlos_param(TOA[0], angles)
plot_graph(p, alpha, s_1x, s_1y)