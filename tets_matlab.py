import matlab.engine as mtlab

eng = mtlab.start_matlab()
[w, pow] = eng.rootmusic()
eng.quit()

print(w)
print(pow)