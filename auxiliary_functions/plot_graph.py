import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np



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

def plot_graph(q, r, p, alpha, Theta_F_LOS, phi_F_LOS, graph_title, s = None):
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

    if s is not None:
        scatter_rect = rotate_rectangle(s, 2, 2, 0)
        scatter_polygon = Polygon(scatter_rect, closed=True, color='yellow', alpha=0.5, label='Scatter')
        ax.add_patch(scatter_polygon)
        ax.scatter(*s, color='yellow')
        ax.plot([r[0], s[0]], [r[1], s[1]], linestyle='--', color='green', label="RIS to Scatter")
        ax.plot([s[0], p[0]], [s[1], p[1]], linestyle='--', color='purple', label="Scatter to UE")
        ax.plot([r[0], p[0]], [r[1], p[1]], linestyle='--', color='orange', label="RIS to UE")
        
        ax.set_title("Canale di Comunicazione con RIS - Caso NLOS "+graph_title)
    else:    
        ax.set_title("Canale di Comunicazione con RIS - Caso LOS "+graph_title)
        # Connect points with lines
        ax.plot([r[0], p[0]], [r[1], p[1]], linestyle='--', color='green', label="RIS to UE")

    # Connect points with lines
    ax.plot([q[0], r[0]], [q[1], r[1]], linestyle='--', color='blue', label="BS to RIS")
    
    # Graph details
    ax.set_xlabel("Coordinate X")
    ax.set_ylabel("Coordinate Y")
    ax.axhline(0, color='black', linewidth=0.5)
    ax.axvline(0, color='black', linewidth=0.5)
    ax.grid(color='gray', linestyle='--', linewidth=0.5)
    ax.legend()
    plt.show(block=False)