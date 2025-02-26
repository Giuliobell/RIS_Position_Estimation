import numpy as np
import matlab.engine


__L = None
__N = None
__M = None
__K = None
__D = None

class MBCEAlgorithm:
    def __init__(self, L, N, M, K, D):
        self.L = L
        self.N = N
        self.M = M
        self.K = K
        self.D = D
        self.eng = matlab.engine.start_matlab()



    #****** DICHIARAZIONE DELLE FUNZIONI ******
    def __khatri_rao_prod(self, A, B):
        if A.shape[1] != B.shape[1]:
            raise ValueError("Le matrici devono avere lo stesso numero di colonne")
    
        # Applica il prodotto di Kronecker colonna per colonna
        return np.hstack([np.kron(A[:, i], B[:, i]).reshape(-1, 1) for i in range(A.shape[1])])

    def __initialize_A_N_matrix(self, angle, row, column):  
        matr = np.zeros((row, column), dtype=np.complex128)
    
        for i in range(row*column):
            matr[i // column, i % column] = np.exp(-1j * (i // column) * np.pi * np.cos(angle)) 
    
        return matr

    def __initialize_Gamma_Matrix(self, gain, row, column):
        matr = np.zeros((row, column), dtype=np.complex128)
    
        for i in range(row):
            matr[i, i] =  np.random.normal(0,np.sqrt(1), 1)[0]
    
        return matr

    def __initialize_omega(self, dim):
        arr = np.zeros(dim, dtype=np.complex128)
    
        for i in range(dim):
            arr[i] = np.exp(1j*2*np.pi*np.random.uniform(0,2*np.pi))
        
        return arr.T
    
    def __initialize_noise_Matrix(self, row, column, power):
        matr = np.zeros((row, column), dtype=np.complex128)
    
        for i in range(row*column):
            matr[i // column, i % column] = np.random.normal(0,np.sqrt(power), 1)[0]
    
        return matr

    def __dft_matrix(self, N):
        #Crea una matrice DFT di dimensione NxN.
        n = np.arange(N)
        k = n.reshape((N, 1))
        omega = (1/np.sqrt(N)) * np.exp(-2j * np.pi * k * n / N)
        return omega

    def __initialize_K_Matrix(self, D,N, Theta):
        matr = np.zeros((N,D), dtype=np.complex128)
    
        matr[0:D, 0:D] = Theta
    
        return matr

    def __estract_angles(self, cos_estim):
        # Assumendo che gli angoli di RIS siano compresi nell'intevallo  [30°,150°]
        start_angle = np.pi / 6
        end_angle = 5 * np.pi / 6
    
        best_phi_G = None
        best_Theta_F = None
        min_diff = float('inf')
    
        for phi_G in np.linspace(start_angle, end_angle, 500):
            for Theta_F in np.linspace(start_angle, end_angle, 500):
                diff = abs(np.cos(phi_G) - np.cos(Theta_F) - cos_estim)
                if diff < min_diff:
                    min_diff = diff
                    best_phi_G = phi_G
                    best_Theta_F = Theta_F
                
        return best_phi_G, best_Theta_F

#****** PARTE PRINCIPALE DEL CODICE ******
    def run(self, Theta_F, Theta_G, phi_F, phi_G, power):

        #++++ Modello di Sistema ++++
        #Inizializzazione della Matrice F
        A_N_Theta_F = self.__initialize_A_N_matrix(Theta_F, self.N, self.L)
        A_K_phi_F = self.__initialize_A_N_matrix(phi_F, self.K, self.L)
        Gamma_F = self.__initialize_Gamma_Matrix(1, self.L, self.L)
        F = A_N_Theta_F @ Gamma_F @ A_K_phi_F.T.conjugate()

        #Inizializzazione della Matrice G
        A_M_Theta_G = self.__initialize_A_N_matrix(Theta_G, self.M, self.L)
        A_N_phi_G = self.__initialize_A_N_matrix(phi_G, self.N, self.L)
        Gamma_G = self.__initialize_Gamma_Matrix(1, self.L, self.L)
        G = A_M_Theta_G @ Gamma_G @ A_N_phi_G.T.conjugate()

        #Costruione della Matrice Omega
        omega = self.__initialize_omega(self.N)
        Omega = np.diag(omega)

        #Costruzione della matrice H
        H = G @ Omega @ F

        #++++ Stima dei Parametri del Canale ++++
        #Calcolo Angoli Esterni
        #Theta_G
        N_noise = self.__initialize_noise_Matrix(self.M,self.K, power)
        V = H + N_noise 

        R_VV = V @ V.conjugate().T
        estimated_angle, pow = self.eng.rootmusic(R_VV, 1, 'corr', nargout=2)
        Theta_G_hat = np.arccos(-estimated_angle/np.pi)

        #phi_F
        R_HH = V.conj().T @ V
        estim_2, pow = self.eng.rootmusic(R_HH, 1, 'corr', nargout=2)
        phi_F_hat = np.arccos(-estim_2/np.pi)


        #Precoder e Combiner
        P = self.__initialize_A_N_matrix(phi_F_hat, self.K, self.L)
        W = self.__initialize_A_N_matrix(Theta_G_hat, self.M, self.L)

        #Calcolo Angoli Interni
        Theta = self.__dft_matrix(self.D)
        K = self.__initialize_K_Matrix(self.D,self.N, Theta)

        Gamma = np.kron(Gamma_F,Gamma_G)
        A_N_psi = self.__khatri_rao_prod(A_N_Theta_F.conj().T, A_N_phi_G.T)
        A_N_psi = A_N_psi.T

        Z = np.kron(P.T@A_K_phi_F.conj(), (W.conj().T@A_M_Theta_G)@Gamma@A_N_psi.conj().T)

        N_noise_primo = self.__initialize_noise_Matrix(self.L*self.L,self.D, power)

        Y = Z@K +N_noise_primo

        Z_hat = (1/self.D) * Y @ Theta.conj().T
        R_ZZ = Z_hat.conj().T@Z_hat

        #B_phi_F_Theta_G = np.kron(A_K_phi_F.conj(),A_M_Theta_G)
        
        #A_N_psi = self.__khatri_rao_prod(A_N_Theta_F.conj().T, A_N_phi_G.T).T

        #Gamma = np.kron(Gamma_F,Gamma_G)
        
        #N_noise_primo = self.__initialize_noise_Matrix(self.M*self.K,self.N, power)
        
        #E = B_phi_F_Theta_G @ Gamma @ A_N_psi.conj().T + N_noise_primo
        #E = self.__khatri_rao_prod(F.T, G)
        
        #R_ZZ = E.conj().T @ E

        #Stima angoli interni
        estim_3, pow = self.eng.rootmusic(R_ZZ, 1, 'corr', nargout=2)

        phi_G_hat, Theta_F_hat = self.__estract_angles(np.cos(abs(estim_3)))

        self.eng.quit()
        
        return Theta_G_hat, phi_F_hat, phi_G_hat, Theta_F_hat