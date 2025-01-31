#NOTE: Metti verifica che i parametri forniti in input rispettino tutti i vincoli richiesti senno lancia eccezzioni custom

import numpy as np
import matlab.engine as mtlab

class InvalidArgumentLength(Exception):
    def __init__(self, message):
        super().__init__(message)

#COSTANTI
K = 10    #Numero di antenne UE
M = 10    #Numero di antenne BS
N = 10    #Numero di antenne RIS
L= 1      #Numero di percorsi usati (L tra UE -> RIS = L tra RIS -> BS)   

#********USA PIU' PERCORSI PER IL CASO NLOS********
#*********CAPISCI COME FUNZIONANO GLI ANGOLI************
Theta_G = [0.9273]
Phi_G = [0.6435]
Theta_F = [1.1071]
Phi_F = [0.4636]

distance_BS_RIS = 50
distance_RIS_UE = 44.72135955

#***FUNZIONI***
def initialize_Gamma(dim, distance):
    matr = np.zeros((dim, dim), dtype=np.complex128)
    
    for i in range(dim):
        matr[i, i] = np.exp(1j*2*np.pi*np.random.uniform(0,1)) / distance
    
    return matr

def initialize_omega(n_ris_elems):
    matr = np.zeros((n_ris_elems, n_ris_elems), dtype=np.complex128)
    
    for i in range(n_ris_elems):
        matr[i, i] = np.exp(1j*2*np.pi*np.random.uniform(0,1))
    
    return matr

def initialize_A_matrix(row, column, angles):
    matr = np.zeros((row, column), dtype=np.complex128)
    
    #row = idx // cols  # Calculate the row index
    #col = idx % cols   # Calculate the column index
    
    for i in range(row*column):
        matr[i // column, i % column] = np.exp(-1j * (i // column) * np.pi * np.cos(angles[i % L]))
    
    return matr

def calculate_R_SS(matr,J, U, S):
    result = np.zeros((S, S), dtype=complex)  

    for u in range(U-1):
        submatrix = matr[u:u + S, u:u + S]
        term = submatrix + J @ submatrix.conjugate() @ J
        result += term

    result /= (2 * U)

    return result

def initialize_Theta(dim):
    matr = np.zeros((dim, dim), dtype=np.complex128)
    
    for i in range((dim*dim) - 1):
        matr[i // dim, i % dim] = np.exp(-1j * 2 * np.pi * ((i // dim) * (i % dim) / (dim*dim))) * (1/np.sqrt(dim*dim))
        
    return matr





# Funzione per generare configurazioni del RIS
def generate_ris_configurations(n_ris, num_configs):
    configurations = []
    for _ in range(num_configs):
        theta = np.exp(1j * 2 * np.pi * np.random.uniform(0, 1, size=(n_ris, n_ris)))
        configurations.append(theta)
    return configurations

# Funzione per calcolare la matrice delle osservazioni
def compute_observations(h_0, ris_configs, x_t, x_r):
    observations = []
    for theta in ris_configs:
        y = x_r.conj().T @ h_0 @ theta @ x_t
        observations.append(y)
    return np.array(observations)

# Funzione per stimare gli angoli con CEARC
def estimate_angles_cearc(observations, ris_configs, num_iterations=10):
    angles = []
    for _ in range(num_iterations):
        for theta in ris_configs:
            correlation = np.sum(np.abs(observations @ theta.conj().T), axis=0)
            estimated_angle = np.argmax(correlation)
            angles.append(estimated_angle)
    return np.array(angles)

# Funzione per affinare le configurazioni del RIS
def refine_ris_configurations(ris_configs, estimated_angles):
    refined_configs = []
    for theta, angle in zip(ris_configs, estimated_angles):
        refined_theta = theta * np.exp(1j * angle)
        refined_configs.append(refined_theta)
    return refined_configs






#Calculate inputs
print("***COMPUTE INPUTS***")
X_T = np.identity(K, dtype=np.complex128)         #Dimensione K*K
X_R = np.identity(M, dtype=np.complex128)         #Dimensione M*M
Gamma_G = initialize_Gamma(L, distance_BS_RIS)     #Guadagno di percorso RIS_BS dimensione L*L
Gamma_F = initialize_Gamma(L, distance_RIS_UE)     #Guadagno di percorso UE_RIS dimensione L*L

#Perform step 2
print("***PERFORM STEP 2***")
Omega_0 = initialize_omega(N)

A_R = initialize_A_matrix(M, L, Theta_G) 
A_I_H = initialize_A_matrix(N,L, Phi_G)
A_I_H = A_I_H.conj().T
A_I = initialize_A_matrix(N, L , Theta_F)
A_T_H = initialize_A_matrix(K, L, Phi_F)
A_T_H = A_T_H.conj().T

H_0 = A_R @ Gamma_G @ A_I_H @ Omega_0 @ A_I @ Gamma_F @ A_T_H

C_0 = X_R.conj().T @ H_0 @ X_T

H_0 = np.linalg.inv(X_R.conj().T) @ C_0 @ np.linalg.inv(X_T)


#Perform step 3
print("***PERFORM STEP 3***")
R_Theta_G = (1/M)*H_0*H_0.conj().T           #Dubbio su 1/K
J = np.fliplr(np.identity(M-L))
R_Theta_G_SS = calculate_R_SS(R_Theta_G, J, L+1, M-L)


#Perform step 4
print("***PERFORM STEP 4***")
eng = mtlab.start_matlab()
Theta_G_Ext, pow_Theta_G_Ext = eng.rootmusic(R_Theta_G_SS, 1, nargout=2)


#Perform step 5
print("***PERFORM STEP 5***")
R_Phi_F = (1/K)*H_0*H_0.conj().T
J = np.fliplr(np.identity(K-L))
R_Phi_F_SS = calculate_R_SS(R_Theta_G, J, L+1, K-L)


#Perform step 6
print("***PERFORM STEP 6***")
Phi_F_Ext, pow_Phi_F_Ext = eng.rootmusic(R_Phi_F_SS, 1, nargout=2)
eng.quit()

#Stima angoli interni implementata con algoritmo CEARC
print("***PERFORM STEP 7***")






# Parametri iterativi
max_iterations = 10
tolerance = 1e-3
num_configs = 1


# Inizializza configurazioni RIS
ris_configs = generate_ris_configurations(N, num_configs)
print(ris_configs)

# Processo iterativo
previous_angles = None
for iteration in range(max_iterations):
    #print(f"*** ITERAZIONE {iteration + 1} ***")

    # Step 1: Calcolo delle osservazioni
    observations = compute_observations(H_0, ris_configs, X_T, X_R)

    # Step 2: Stima angoli con CEARC
    estimated_angles = estimate_angles_cearc(observations, ris_configs)

    # Convergenza: verifica cambiamenti rispetto all'iterazione precedente
    if previous_angles is not None:
        angle_diff = np.linalg.norm(np.array(estimated_angles) - np.array(previous_angles))
        #print(f"Differenza tra angoli stimati: {angle_diff}")
        if angle_diff < tolerance:
            #print("Convergenza raggiunta.")
            break

    # Step 3: Affinamento configurazioni RIS
    ris_configs = refine_ris_configurations(ris_configs, estimated_angles)

    # Aggiorna gli angoli precedenti
    previous_angles = estimated_angles


reshaped_angles = np.array(estimated_angles).reshape(-1, num_configs)
final_internal_angles = reshaped_angles[-1]

# Risultati finali
print("*** RISULTATI FINALI ***")

print("Angoli di Arrivo (AoA) stimati:", Theta_G_Ext)
print("Angoli di Partenza (AoD) stimati:", Phi_F_Ext)


print("Angoli stimati finali:", np.radians(final_internal_angles))
#print("Configurazioni RIS finali:", ris_configs)
