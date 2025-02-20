import matlab.engine
import numpy as np

# Avvia il motore MATLAB
eng = matlab.engine.start_matlab()

# Parametri fissi per il test
N = 8  # Numero di antenne riceventi
M = 1000  # Numero di snapshot
num_sources = 2  # Numero di sorgenti
angles_true = np.array([-20, 35])  # Angoli di arrivo fissi (in gradi)

# Costanti per la simulazione
lambda_wave = 1  # Lunghezza d'onda normalizzata
d = lambda_wave / 2  # Spaziatura tra le antenne
theta_rad = np.radians(angles_true)  # Conversione in radianti

# Matrice di steering (array response vectors)
A = np.exp(-1j * 2 * np.pi * d * np.arange(N).reshape(-1, 1) @ np.sin(theta_rad).reshape(1, -1))

# Segnale trasmesso fisso per garantire risultati ripetibili
S_fixed = np.ones((num_sources, M)) + 1j * np.ones((num_sources, M))

# Segnale ricevuto senza casualità
X = A @ S_fixed

# Matrice di correlazione corretta
R = (X @ X.conj().T) / M

# **Passa la matrice COMPLESSA correttamente a MATLAB**
R_real = np.real(R).tolist()  # Parte reale
R_imag = np.imag(R).tolist()  # Parte immaginaria
R_matlab = eng.complex(matlab.double(R_real), matlab.double(R_imag))  # CORRETTA costruzione in MATLAB

# **Chiamata a rootmusic di MATLAB**
estimated_angles = eng.rootmusic(R_matlab, num_sources, 'corr')  # Output in gradi

# Converti il risultato in un array NumPy
angles_python = np.array(np.rad2deg(estimated_angles))


# Output
print("Angoli reali:", angles_true)
print("Angoli stimati con MATLAB:", angles_python)

# Chiudi MATLAB Engine
eng.quit()
