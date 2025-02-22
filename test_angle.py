import numpy as np
import matlab.engine 




angle = 0.4567

matr = np.zeros((5, 1), dtype=np.complex128)
for i in range(5):
    matr[i, 0] = np.exp(-1j * i *angle)

print(matr)


matlab_engine = matlab.engine.start_matlab()
estimated_angle, pow = matlab_engine.rootmusic(matr@matr.conj().T, 1, 'corr', nargout=2)
matlab_engine.quit()


print(abs(estimated_angle))
print(pow)