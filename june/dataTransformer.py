import numpy as np



### REWRITE TRAINING DATA INTO MORE PYTORCH FRIENDLY FORMAT
### FOR TRAINING

file_path = 'june_data.txt'


with open(file_path, 'r') as f:
    lines = f.readlines()

cutoff = int(len(lines)/2)

# kick rate 1000
data_1000 = [line.strip().split(' ') for line in lines[3:cutoff+3]]
# kick rate 3000
data_3000 = [line.strip().split(' ') for line in lines[cutoff+3:2*cutoff+3]]
# kick rate 5000
data_5000 = [line.strip().split(' ') for line in lines[2*cutoff+3:]]

dta = [coupled_dta for coupled_dta in zip(data_1000, data_3000, data_5000)]


new_dta = []
for index, e in enumerate(dta):
  E_1000, I_1000 = float(e[0][8]), float(e[0][9])
  E_3000, I_3000 = float(e[1][8]), float(e[1][9])
  E_5000, I_5000 = float(e[2][8]), float(e[2][9])

  SEE, SIE, SEI, SII = float(e[0][0]), float(e[0][1]), float(e[0][2]), \
  float(e[0][3])

  new_dta.append([E_1000, I_1000, E_3000, I_3000,
                         E_5000, I_5000,
                         SEE, SIE, SEI, SII])

print('new_dta', new_dta)

np.savetxt('data_june.txt', new_dta)
