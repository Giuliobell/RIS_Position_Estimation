from Algorithms.angle_estimation_algorithms import MBCEAlgorithm
from Algorithms.localization_algorithms import localization_algorithms
import auxiliary_functions.error_calculation as err_calc
import auxiliary_functions.plot_graph as plot
from tabulate import tabulate
import numpy as np
import matplotlib.pyplot as plt

#*****************************************STIMA DEGLI ANGOLI*****************************************
L = 1     #Numero di percorsi
N = 16    #Numero di antenne RIS
M = 4     #Numero di antenne BS
K = 4     #Numero di antenne UE
D = 15    #Numero di configurazioni di RIS

Theta_F_0 = 1.1071
Theta_F_1 = 1.3734
Theta_G = 0.9273
phi_F_0 = 0.4036
phi_F_1 = 0.7853
phi_G = 0.6435
power = 0

c = 299792458           # Velocità della luce
q = np.array([0, 0])    # Posizione BS
r = np.array([30, 40])  # Posizione RIS
p = np.array([70,20])   # Posizione UE
s = np.array([55,35])

tau_0 = 1.49174e-7      #Distanza RIS-UE nel caso LOS
tau_1 = 1.5580e-7

RIS_Rotation = 3/2*np.pi



num_iterations = 5  # Numero di ripetizioni

# Genera 10 valori di SNR tra -30dB e 0dB
samples = 10
dB_values = np.linspace(-30, 0, samples)

# Inizializza array vuoti per i risultati
iteration_params_extimated = []
iteration_params_extimated_error = []

for i in range(num_iterations):
    for val in dB_values:
        # Stima angoli
        print(f"Iterazione {i} SNR {val} power: {10**(val/10)}")
        mbce = MBCEAlgorithm(L, N, M, K, D)
        Theta_G_hat, phi_F_0_hat, phi_G_hat, Theta_F_0_hat = mbce.run(Theta_F_0, Theta_G, phi_F_0, phi_G, 10**(val/10))
        _, phi_F_1_hat, _, Theta_F_1_hat = mbce.run(Theta_F_1, Theta_G, phi_F_1, phi_G, 10**(val/10))

        loc = localization_algorithms(c)
        #p_hat, alpha_hat, s_hat = loc.nlos_path_estimation(r, Theta_F_0_hat, Theta_F_1_hat, phi_F_0_hat, phi_F_1_hat, tau_0)
        p_hat, alpha_hat, s_hat = loc.nlos_path_estimation(r, RIS_Rotation + Theta_F_0_hat, RIS_Rotation + Theta_F_1_hat, phi_F_0_hat, phi_F_1_hat, tau_0)

        # Salva i valori stimati
        iteration_params_extimated.append([val, Theta_G_hat, phi_F_0_hat, phi_F_1_hat, phi_G_hat, Theta_F_0_hat, Theta_F_1_hat, p_hat[0], p_hat[1], s_hat[0], s_hat[1], alpha_hat])

        # Calcola gli errori percentuali
        iteration_params_extimated_error.append([
            val,
            err_calc.percentage_error(Theta_G, Theta_G_hat),
            err_calc.percentage_error(phi_F_0, phi_F_0_hat),
            err_calc.percentage_error(phi_F_1, phi_F_1_hat),
            err_calc.percentage_error(phi_G, phi_G_hat),
            err_calc.percentage_error(Theta_F_0, Theta_F_0_hat),
            err_calc.percentage_error(Theta_F_1, Theta_F_1_hat),
            err_calc.percentage_error(p[0], p_hat[0]),
            err_calc.percentage_error(p[1], p_hat[1]),
            err_calc.percentage_error(s[0], s_hat[0]),
            err_calc.percentage_error(s[1], s_hat[1])
        ])

# Calcola le medie degli angoli stimati e degli errori per ogni livello di rumore
mean_params_extimated = []
mean_params_extimated_error = []

for val in dB_values:
    filtered_estimations = [x for x in iteration_params_extimated if x[0] == val]
    filtered_errors = [x for x in iteration_params_extimated_error if x[0] == val]

    if filtered_estimations:
        mean_estimations = np.mean(filtered_estimations, axis=0)
    else:
        mean_estimations = [val] + [np.nan] * (len(iteration_params_extimated[0]) - 1)

    if filtered_errors:
        mean_errors = np.mean(filtered_errors, axis=0)
    else:
        mean_errors = [val] + [np.nan] * (len(iteration_params_extimated_error[0]) - 1)

    mean_params_extimated.append(mean_estimations)
    mean_params_extimated_error.append(mean_errors)

# Stampa i risultati in tabella
print(tabulate(np.round(mean_params_extimated, 4), headers=["SNR", "Theta_G_hat", "phi_F_0_hat" , "phi_F_1_hat", "phi_G_hat", "Theta_F_0_hat", "Theta_F_1_hat", "p_hat_x", "p_hat_y", "s_hat_x", "s_hat_y" , "alpha_hat"], tablefmt="pretty", floatfmt=".4f"))
print(tabulate(np.round(mean_params_extimated_error, 4), headers=["SNR", "Theta_G", "phi_F_0" , "phi_F_1", "phi_G", "Theta_F_0", "Theta_F_1", "p_x", "p_y", "s_x", "s_y"], tablefmt="pretty", floatfmt=".4f"))


# Plot dei risultati
mean_params_extimated = np.array(mean_params_extimated)
mean_params_extimated_error = np.array(mean_params_extimated_error)

plt.figure(figsize=(12, 8))

# Plot degli errori percentuali degli angoli
plt.subplot(2, 1, 1)
plt.plot(dB_values, mean_params_extimated_error[:, 1], label="Theta_G_error")
plt.plot(dB_values, mean_params_extimated_error[:, 2], label="phi_F_0_error")
plt.plot(dB_values, mean_params_extimated_error[:, 3], label="phi_F_1_error")
plt.plot(dB_values, mean_params_extimated_error[:, 4], label="phi_G_error")
plt.plot(dB_values, mean_params_extimated_error[:, 5], label="Theta_F_0_error")
plt.plot(dB_values, mean_params_extimated_error[:, 6], label="Theta_F_1_error")
plt.xlabel("SNR (dB)")
plt.ylabel("Errori Percentuali (%)")
plt.title("Medie degli Errori Percentuali degli Angoli")
plt.legend()
plt.grid(True)

# Plot degli errori percentuali delle posizioni
plt.subplot(2, 1, 2)
plt.plot(dB_values, mean_params_extimated_error[:, 7], label="p_x_error")
plt.plot(dB_values, mean_params_extimated_error[:, 8], label="p_y_error")
plt.plot(dB_values, mean_params_extimated_error[:, 9], label="s_x_error")
plt.plot(dB_values, mean_params_extimated_error[:, 10], label="s_y_error")
plt.xlabel("SNR (dB)")
plt.ylabel("Errori Percentuali (%)")
plt.title("Medie degli Errori Percentuali delle Posizioni")
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()



# Stampa i risultati
#print("Valori Ottenuti")
#print(tabulate([["Theta_G", Theta_G, Theta_G_hat],
#                ["phi_F", phi_F, phi_F_hat],
#                ["phi_G", phi_G, phi_G_hat],
#                ["Theta_F", Theta_F, Theta_F_hat]],
#                headers=[" ", "Valore Reale", "Valore Stimato"], tablefmt="pretty"))
#
#print("Errori Percentuali")
#print(tabulate([["Theta_G", err_calc.percentage_error(Theta_G, Theta_G_hat)],
#                ["phi_F", err_calc.percentage_error(phi_F, phi_F_hat)],
#                ["phi_G", err_calc.percentage_error(phi_G, phi_G_hat)],
#                ["Theta_F", err_calc.percentage_error(Theta_F, Theta_F_hat)]],
#                headers=["Angolo", "Errore Percentuale"], tablefmt="pretty"))
#
#print("Valori Ottenuti")
#print(tabulate([["Posizione UE (X) (angoli reali)", p[0], p_hat[0]],
#                ["Posizione UE (Y) (angoli reali)", p[1], p_hat[1]],
#                ["Posizione UE (X) (angoli stimati)", p[0], p_hat_estim_angle[0]],
#                ["Posizione UE (Y) (angoli stimati)", p[1], p_hat_estim_angle[1]],
#                ["Angolo alpha angoli reali", "", RIS_Rotation + alpha],
#                ["Angolo alpha angoli stimati", "", RIS_Rotation + alpha_estim_angle]],
#                headers=["Valore Reale", "Valore Stimato"], tablefmt="pretty"))
#print("Errori Percentuali")
#print(tabulate([["Posizione UE (X) (angoli reali)", err_calc.percentage_error(p[0], p_hat[0])],
#                ["Posizione UE (Y) (angoli reali)", err_calc.percentage_error(p[1], p_hat[1])],
#                ["Posizione UE (X) (angoli stimati)", err_calc.percentage_error(p[0], p_hat_estim_angle[0])],
#                ["Posizione UE (Y) (angoli stimati)", err_calc.percentage_error(p[1], p_hat_estim_angle[1])]],
#                headers=["Componente", "Errore Percentuale"], tablefmt="pretty"))
#

#if input("Stampare i grafici [Y/N]: ") == "Y":
#    plot.plot_graph(q, r, p_hat_estim_angle, RIS_Rotation + alpha_estim_angle, Theta_F_hat, phi_F_hat, "Posizione Stimata Angoli Stimati")