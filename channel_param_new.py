import numpy as np

#****** DICHIARAZIONE DELLE COSTANTI ******
L = 1     #Numero di percorsi
N = 16    #Numero di antenne RIS
M = 4     #Numero di antenne BS
K = 4     #Numero di antenne UE

Theta_F = 1.1071
Theta_G = 0.9273
phi_F = 0.4036
phi_G = 0.6435
power = 1

#****** DICHIARAZIONE DELLE FUNZIONI ******
def khatri_rao_prod(A, B):
    if A.shape[1] != B.shape[1]:
        raise ValueError("Le matrici devono avere lo stesso numero di colonne")
    
    # Applica il prodotto di Kronecker colonna per colonna
    return np.hstack([np.kron(A[:, i], B[:, i]).reshape(-1, 1) for i in range(A.shape[1])])


def initialize_A_N_matrix(angle, row, column):
    matr = np.zeros((row, column), dtype=np.complex128)
    
    for i in range(row*column):
        matr[i // column, i % column] = np.exp(-1j * (i // column) * np.pi * np.cos(angle)) 
    
    return matr

def initialize_Gamma_Matrix(gain, row, column):
    matr = np.zeros((row, column), dtype=np.complex128)
    
    for i in range(row):
        matr[i, i] =  np.random.normal(0,np.sqrt(1), 1)[0]
    
    return matr

def initialize_omega(dim):
    arr = np.zeros(dim, dtype=np.complex128)
    
    for i in range(dim):
        arr[i] = np.exp(1j*2*np.pi*np.random.uniform(0,2*np.pi))
        
    return arr.T
    
def initialize_noise_Matrix(row, column):
    matr = np.zeros((row, column), dtype=np.complex128)
    
    for i in range(row*column):
        matr[i // column, i % column] = np.random.normal(0,np.sqrt(power), 1)[0]
    
    return matr

def initialize_pilotVector_matrix(row, column):
    matr = np.zeros((row, column), dtype=np.complex128)
    for i in range(row*column):
        matr[i // column, i % column] = 1
    
    return matr



#****** PARTE PRINCIPALE DEL CODICE ******
#++++ Modello di Sistema ++++
#Inizializzazione della Matrice F
A_N_Theta_F = initialize_A_N_matrix(Theta_F, N, L)
A_K_phi_F = initialize_A_N_matrix(phi_F, K, L)
Gamma_F = initialize_Gamma_Matrix(1, L, L)
F = A_N_Theta_F @ Gamma_F @ A_K_phi_F.T.conjugate()

#Inizializzazione della Matrice G
A_M_Theta_G = initialize_A_N_matrix(Theta_F, M, L)
A_N_phi_G = initialize_A_N_matrix(phi_G, N, L)
Gamma_G = initialize_Gamma_Matrix(1, L, L)
G = A_M_Theta_G @ Gamma_G @ A_N_phi_G.T.conjugate()

#Costruione della Matrice Omega
omega = initialize_omega(N)
Omega = np.diag(omega)

#Costruzione della matrice H
H = G @ Omega @ F

#Calcolo della vettorizzazione della matrice H
vec_H = H.ravel(order = 'F')
vec_H = vec_H.reshape(-1, 1)

#++++ Stima dei Parametri del Canale ++++
#Calcolo Angoli Esterni
X = initialize_pilotVector_matrix(K, K)
N_noise = initialize_noise_Matrix(M,K)
V = H @ X + N_noise

R_VV = V @ V.conjugate().T






np.set_printoptions(linewidth=np.inf,  precision=1)
print(V)