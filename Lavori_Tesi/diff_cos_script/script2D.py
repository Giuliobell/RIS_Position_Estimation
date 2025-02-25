import numpy as np
import matplotlib.pyplot as plt

# Creazione intervallo di valori per x tra 0 e 2π, prendendone 500
x_vals = np.linspace(0, 2 * np.pi, 500)

# Calcolo dei valori dei coseni e della differenza tra essi
y_fixed = np.pi / 4  # per esempio
cos_x = np.cos(x_vals)
cos_y_fixed = np.cos(y_fixed)
diff_cos = cos_x - cos_y_fixed

# Plot dei grafici di cos(x), cos(y) e della differenza cos(x) - cos(y)
plt.figure(figsize=(10, 6))
plt.plot(x_vals/np.pi, cos_x/np.pi, label=r'$\cos(x)$', color='blue')
plt.axhline(y=cos_y_fixed, color='green', linestyle='--', label=r'$\cos(y)$ con $y = \frac{\pi}{4}$')
plt.plot(x_vals/np.pi, diff_cos/np.pi, label=r'$\cos(x) - \cos(y)$', color='red')

# Aggiunta di legenda ed etichette per i grafici plottati
plt.xlabel(r'$x$')
plt.ylabel('Valori delle funzioni')
plt.title(r'Grafico di $\cos(x)$, $\cos(y)$ e $\cos(x) - \cos(y)$ con $y = \frac{\pi}{4}$')
plt.legend()
plt.grid(True)
plt.show()
