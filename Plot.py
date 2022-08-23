import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter

lri_files = ["ice", "AC", "silicate"]

for lri_file_name in lri_files:
    lri_file = np.genfromtxt(f"{lri_file_name}.RI")
    fig = plt.figure()
    ax = fig.add_subplot(1, 3, 1)
    ax.set_xscale('log')
    l = lri_file[:,0]
    r = lri_file[:,1]
    i = lri_file[:,2]
    plt.plot(l,r, color='g')
    plt.plot(l,i, color='m')
    plt.title(f'Refractive index, {lri_file_name}')
    plt.xlabel(r'${\lambda},{\mu m}$')

    Temps = [3500, 6000, 10000, 18000, 24000, 30000, 50000]

    a = np.linspace(0, 1, 20)
    q_file = np.genfromtxt(f"results_Q_{lri_file_name}", delimiter=',')
    ax = fig.add_subplot(1, 3, 2)
    plt.title(r'${\overline{Q}_{pr}}$' + f', {lri_file_name}')
    plt.xlabel(r'$a,{\mu m}$')
    for i in range(len(Temps)):
        plt.plot(a, q_file[i], label = f"T={Temps[i]}")
        plt.legend()

    lablels = ["A5 I", "A5 V", "G0 I","G0 III", "G0 V", "M0 I", "M0 III", "M0 V"]
        
    b_file = np.genfromtxt(f"results_b_{lri_file_name}", delimiter=',')
    ax = fig.add_subplot(1, 3, 3)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    plt.title(r'${\beta}$' + f', {lri_file_name}')
    plt.xlabel(r'$a,{\mu m}$')
    for i in range(len(lablels)):
        plt.plot(a, b_file[i], label = f"T={lablels[i]}")
        plt.legend()
    plt.show()