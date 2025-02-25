from Algorithms.angle_estimation_algorithms import MBCEAlgorithm
from tabulate import tabulate

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
Theta_G_hat, phi_F_hat, phi_G_hat, Theta_F_hat = mbce.run(Theta_F, Theta_G, phi_F, phi_G, power)


# Stampa i risultati
print("Valori Ottenuti")
print(tabulate([["Theta_G", Theta_G, Theta_G_hat],
                ["phi_F", phi_F, phi_F_hat],
                ["phi_G", phi_G, phi_G_hat],
                ["Theta_F", Theta_F, Theta_F_hat]],
                headers=[" ", "Valore Reale", "Valore Stimato"], tablefmt="pretty"))

print("Errori Percentuali")
print(tabulate([["Theta_G", abs(Theta_G - Theta_G_hat) / Theta_G * 100],
                ["phi_F", abs(phi_F - phi_F_hat) / phi_F * 100],
                ["phi_G", abs(phi_G - phi_G_hat) / phi_G * 100],
                ["Theta_F", abs(Theta_F - Theta_F_hat) / Theta_F * 100]],
                headers=["Angolo", "Errore Percentuale"], tablefmt="pretty"))