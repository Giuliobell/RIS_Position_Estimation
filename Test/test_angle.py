import numpy as np
import matlab.engine 




angle = 0.4567

matr = np.zeros((4, 4), dtype=np.complex128)
for i in range(4):
    for j in range(4):
        matr[i, j] = np.exp(-1j * i * angle)

np.set_printoptions(precision=4, linewidth=np.inf)

# matr @ matr.conj().T risulta in una matrice 4x4
print("matr @ matr.conj().T:")
print(matr @ matr.conj().T)

print("****************")

# matr.conj().T @ matr risulta in una matrice 4x4
print("matr.conj().T @ matr:")
print(matr.conj().T @ matr)


matlab_engine = matlab.engine.start_matlab()
estimated_angle, pow = matlab_engine.rootmusic(matr@matr.conj().T, 1, 'corr', nargout=2)


print(abs(estimated_angle))
print(pow)

estimated_angle2, pow = matlab_engine.rootmusic(matr.conj().T@matr, 1, 'corr', nargout=2)
matlab_engine.quit()

print(abs(estimated_angle2))