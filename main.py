from Algorithms.angle_estimation_algorithms import MBCEAlgorithm

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

mbce = MBCEAlgorithm(L, N, M, K, D)
Theta_G_hat, phi_F_hat, phi_G, Theta_F = mbce.run(Theta_F, Theta_G, phi_F, phi_G, power)

# Stampa i risultati
print(f"Theta_G_hat: {Theta_G_hat}, phi_F_hat: {phi_F_hat}, phi_G: {phi_G}, Theta_F: {Theta_F}")