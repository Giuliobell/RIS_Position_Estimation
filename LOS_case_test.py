from Algorithms.angle_estimation_algorithms import MBCEAlgorithm
from Algorithms.localization_algorithms import localization_algorithms
import auxiliary_functions.error_calculation as err_calc
import auxiliary_functions.plot_graph as plot
from tabulate import tabulate
import numpy as np

#*****************************************STIMA DEGLI ANGOLI*****************************************
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

#Stima angoli
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
print(tabulate([["Theta_G", err_calc.percentage_error(Theta_G, Theta_G_hat)],
                ["phi_F", err_calc.percentage_error(phi_F, phi_F_hat)],
                ["phi_G", err_calc.percentage_error(phi_G, phi_G_hat)],
                ["Theta_F", err_calc.percentage_error(Theta_F, Theta_F_hat)]],
                headers=["Angolo", "Errore Percentuale"], tablefmt="pretty"))



#*****************************************STIMA DELLA POSIZIONE*****************************************
c = 299792458           # Velocità della luce
q = np.array([0, 0])    # Posizione BS
r = np.array([30, 40])  # Posizione RIS
p = np.array([70,20])   # Posizione UE

tau_0 = 1.49174e-7      #Distanza RIS-UE nel caso LOS

RIS_Rotation = 3/2*np.pi



#Stima della posizione
loc = localization_algorithms(c)
p_hat_estim_angle, alpha_estim_angle = loc.los_path_estimation(r, RIS_Rotation + Theta_F_hat, phi_F_hat, tau_0)
p_hat, alpha = loc.los_path_estimation(r, RIS_Rotation + Theta_F, phi_F, tau_0)

print("Valori Ottenuti")
print(tabulate([["Posizione UE (X) (angoli reali)", p[0], p_hat[0]],
                ["Posizione UE (Y) (angoli reali)", p[1], p_hat[1]],
                ["Posizione UE (X) (angoli stimati)", p[0], p_hat_estim_angle[0]],
                ["Posizione UE (Y) (angoli stimati)", p[1], p_hat_estim_angle[1]],
                ["Angolo alpha angoli reali", "", RIS_Rotation + alpha],
                ["Angolo alpha angoli stimati", "", RIS_Rotation + alpha_estim_angle]],
                headers=["Valore Reale", "Valore Stimato"], tablefmt="pretty"))
print("Errori Percentuali")
print(tabulate([["Posizione UE (X) (angoli reali)", err_calc.percentage_error(p[0], p_hat[0])],
                ["Posizione UE (Y) (angoli reali)", err_calc.percentage_error(p[1], p_hat[1])],
                ["Posizione UE (X) (angoli stimati)", err_calc.percentage_error(p[0], p_hat_estim_angle[0])],
                ["Posizione UE (Y) (angoli stimati)", err_calc.percentage_error(p[1], p_hat_estim_angle[1])]],
                headers=["Componente", "Errore Percentuale"], tablefmt="pretty"))


if input("Stampare i grafici [Y/N]: ") == "Y":
    plot.plot_graph(q, r, p_hat_estim_angle, RIS_Rotation + alpha_estim_angle, Theta_F_hat, phi_F_hat, "Posizione Stimata Angoli Stimati")

