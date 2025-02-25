import numpy as np
from sympy import symbols, Eq, solve, tan

__c = None

class localization_algorithms:
    
    def __init__(self,c):
        self.c = c
    
    #****** DICHIARAZIONE DELLE FUNZIONI ******
    def los_path_estimation(self, r, Theta_F_hat, phi_F_hat, tau_0):
        """Calculate LOS parameters."""
        p = np.add(r, self.c * tau_0 * np.array((np.cos(Theta_F_hat), np.sin(Theta_F_hat))).T)
        alpha = np.pi + Theta_F_hat - phi_F_hat
        return p, alpha

    def nlos_path_Estimation(self, tau_0: float, angles):
        """Calculate Non-Line-of-Sight (NLOS) parameters."""
        p, alpha = self.los_path_estimation(tau_0, angles)

        alpha = np.pi - alpha

        # Define symbolic variables for intersection points
        s_1x, s_1y = symbols('s_1x s_1y')

        # Define equations using SymPy's tan function
        eq_1 = Eq(tan(np.pi - (angles[1][1] + alpha)), (p[1] - s_1y) / (p[0] - s_1x))
        eq_2 = Eq(tan(angles[1][0]), (s_1y - r[1]) / (s_1x - r[0]))

        # Solve equations for s_1x and s_1y
        solutions = solve([eq_1, eq_2], (s_1x, s_1y), dict=True)

        if not solutions:
            raise ValueError("No solution found for the intersection points.")

        solution = solutions[0]  # Select the first solution
        return p, alpha, float(solution[s_1x]), float(solution[s_1y])