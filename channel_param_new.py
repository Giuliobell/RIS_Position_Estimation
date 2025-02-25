import numpy as np
import matlab.engine

#****** DICHIARAZIONE DELLE COSTANTI ******
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

def dft_matrix(N):
    """Crea una matrice DFT di dimensione NxN."""
    n = np.arange(N)
    k = n.reshape((N, 1))
    omega = (1/np.sqrt(N)) * np.exp(-2j * np.pi * k * n / N)
    return omega

def initialize_K_Matrix(D,N, Theta):
    matr = np.zeros((N,D), dtype=np.complex128)
    
    matr[0:D, 0:D] = Theta
    
    return matr

def estract_angles(cos_estim):
    # Assumendo che gli angoli di RIS siano compresi nell'intevallo  [30°,150°]
    start_angle = np.pi / 6
    end_angle = 5 * np.pi / 6
    
    best_phi_G = None
    best_Theta_F = None
    min_diff = float('inf')
    
    for phi_G in np.linspace(start_angle, end_angle, 5000):
        for Theta_F in np.linspace(start_angle, end_angle, 5000):
            diff = abs(np.cos(phi_G) - np.cos(Theta_F) - cos_estim)
            if diff < min_diff:
                min_diff = diff
                best_phi_G = phi_G
                best_Theta_F = Theta_F
                
    return best_phi_G, best_Theta_F

#****** PARTE PRINCIPALE DEL CODICE ******
#++++ Modello di Sistema ++++
#Inizializzazione della Matrice F
A_N_Theta_F = initialize_A_N_matrix(Theta_F, N, L)
A_K_phi_F = initialize_A_N_matrix(phi_F, K, L)
Gamma_F = initialize_Gamma_Matrix(1, L, L)
F = A_N_Theta_F @ Gamma_F @ A_K_phi_F.T.conjugate()

#Inizializzazione della Matrice G
A_M_Theta_G = initialize_A_N_matrix(Theta_G, M, L)
A_N_phi_G = initialize_A_N_matrix(phi_G, N, L)
Gamma_G = initialize_Gamma_Matrix(1, L, L)
G = A_M_Theta_G @ Gamma_G @ A_N_phi_G.T.conjugate()

#Costruione della Matrice Omega
omega = initialize_omega(N)
Omega = np.diag(omega)

#Costruzione della matrice H
H = G @ Omega @ F

#++++ Stima dei Parametri del Canale ++++
#Calcolo Angoli Esterni
#Theta_G
N_noise = initialize_noise_Matrix(M,K)
V = H + N_noise 

eng = matlab.engine.start_matlab()
R_VV = V @ V.conjugate().T
estimated_angle, pow = eng.rootmusic(R_VV, 1, 'corr', nargout=2)
Theta_G_hat = np.arccos(-estimated_angle/np.pi)

#phi_F
R_HH = V.conj().T @ V
estim_2, pow = eng.rootmusic(R_HH, 1, 'corr', nargout=2)
phi_F_hat = np.arccos(-estim_2/np.pi)


#Precoder e Combiner
P = initialize_A_N_matrix(phi_F_hat, K, L)
W = initialize_A_N_matrix(Theta_G_hat, M, L)

#Calcolo Angoli Interni
Theta = dft_matrix(D)
K = initialize_K_Matrix(D,N, Theta)

Gamma = np.kron(Gamma_F,Gamma_G)
A_N_psi = khatri_rao_prod(A_N_Theta_F.conj().T, A_N_phi_G.T)
A_N_psi = A_N_psi.T

Z = np.kron(P.T@A_K_phi_F.conj(), (W.conj().T@A_M_Theta_G)@Gamma@A_N_psi.conj().T)

N_noise_primo = initialize_noise_Matrix(L*L,D)

Y = Z@K +N_noise_primo

Z_hat = (1/D) * Y @ Theta.conj().T
R_ZZ = Z_hat.conj().T@Z_hat

#Stima angoli interni
estim_3, pow = eng.rootmusic(R_ZZ, 1, 'corr', nargout=2)

phi_G, Theta_F = estract_angles(np.cos(abs(estim_3)))

if phi_G is not None and Theta_F is not None:
    print(f"phi_G: {phi_G}, Theta_F: {Theta_F}")
else:
    print("Non è stato possibile trovare una combinazione di angoli che soddisfi l'equazione.")

print("Theta_G_hat: ", Theta_G_hat)
print("phi_F_hat: ", phi_F_hat)
eng.quit()


#valore reale di psi = 1.210294875