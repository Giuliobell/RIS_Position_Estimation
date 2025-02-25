import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from mpl_toolkits.mplot3d import Axes3D

# Creazione di una griglia di valori per x e y tra 0 e 2π prendendo 200 campioni per intervallo
x_vals = np.linspace(0, 2 * np.pi, 200)
y_vals = np.linspace(0, 2 * np.pi, 200)
X, Y = np.meshgrid(x_vals, y_vals)

# Calcolo di c = cos(x) - cos(y) per ogni coppia (x, y) nella griglia
Z = np.cos(X) - np.cos(Y)

# Creazione del grafico 3D
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

#Su x e y plotto multipli di pi greco
ax.plot_surface(X/np.pi, Y/np.pi, Z, cmap='viridis', edgecolor='none')

# Aggiunta di etichette negli assi
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.set_zlabel(r'$c = \cos(x) - \cos(y)$')

#Imposto le etichette degli assi x e y con scala multipla di pi greco 
ax.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
ax.yaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))

plt.show()
