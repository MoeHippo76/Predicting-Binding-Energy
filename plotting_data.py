from matplotlib import pyplot as plt

lines = []
with open('/home/hammy/VSLAT/1_dft/data/ene_beta2.agr') as f:
    lines = f.readlines()

Beta = []
Energy = []

for line in lines:
    data = line.split()
    Beta.append(float(data[0]))
    Energy.append(float(data[1]))


plt.plot(Beta,Energy)
plt.xlabel("Beta")
plt.ylabel("Energy (MeV)")
plt.show()