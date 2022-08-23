import numpy as np
import math
from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt

# Данная программа считает планковские средние факторы лучевого давления
# и эффективности выметания пыллинок для звезд различных типов, после
# выполнения расчетов -- строит графики.

# Функция Планка
def B(l, T):
    h = 6.626e-34
    c = 3e+8
    try:
        e = math.exp( (h*c)/(l*T*1.38e-23) )
    except:
        print(f"l = {l}, T={T}")
        return 0
    return (2*h*c*c)/(l**5) * 1/(e - 1)

# Тут мы считаем Q_pr-ы
def func(l, q_pr, T):
    return q_pr*math.pi*B(l, T)
# Подсчет двух интегралов, тот, что в числителе, и тот
# что в знаменателе
def integrate_upper(l, q_pr, T):
    summa = 0
    for (l_i, q_i) in zip(l, q_pr):
        summa += func(l_i, q_i, T)
    return (l[-1] - l[0])/len(l) * summa

def integrate_lower(l, T):
    summa = 0
    for l_i in l:
        summa += math.pi*B(l_i, T)
    return (l[-1] - l[0])/len(l) * summa

Temps = [8510, 8200, 5550, 5850, 6030, 3650, 3800, 3850]
ri_files = ["ice", "AC", "silicate"]

# вносим результаты в файл
for ri_file_name in ri_files:
   ri_file = np.genfromtxt(ri_file_name + ".RI", skip_header=1, skip_footer=0)
   q_file = np.genfromtxt(ri_file_name + ".Q")
   l = ri_file[:,0]*1e-6
   with open("results_Q_" + ri_file_name, 'w') as outfile:
       for T in Temps:
           Q = [0.0]
           lower = integrate_lower(l, T)
           for a in range(q_file.shape[1]):
               upper = integrate_upper(l, q_file[:,a], T)
               Q.append(upper/lower)
           outfile.write(", ".join(map(str, Q)) + "\n")

lablel = ["A5 I", "A5 V", "G0 I","G0 III", "G0 V", "M0 I", "M0 III", "M0 V"]
Ms = [13, 2, 10, 1, 1.05, 13, 1.2, 0.51]
Rs = [60, 1.7, 120, 6, 1.1, 500, 40, 0.6] 
ros = [0.92, 2.25, 2.285]
a = np.linspace(0.05, 1, 19)

# Тут мы считаем beta
def betta(R, T, M, Q_pr, a, ro):
    M *= 2e30
    R *= 7e8
    a *= 1e-6
    return 2.12e-8 * (R**2 * T**4)/M * Q_pr/(a*ro)

# Вносим результат по beta в файл
for i in range(len(ri_files)):
    q_file = np.genfromtxt("results_Q_{0}".format(ri_files[i]), delimiter = ',')

    with open("results_b_" + ri_files[i], 'w') as outfile:
        for k in range(len(Temps)):
            b = [0.0]
            for j in range(len(a)):
                b.append(betta(Rs[k], Temps[k], Ms[k], q_file[k, j+1], a[j], ros[i]))
            outfile.write(", ".join(map(str, b)) + "\n")

# Строим графики
for ri_file_name in ri_files:
    ri_file = np.genfromtxt(f"{ri_file_name}.RI")
    fig = plt.figure()
    ax = fig.add_subplot(1, 3, 1)
    ax.set_xscale('log')
    l = ri_file[:,0]
    r = ri_file[:,1]
    i = ri_file[:,2]
    plt.plot(l,r, color='g')
    plt.plot(l,i, color='m')
    plt.title(f'Refractive index, {ri_file_name}')
    plt.xlabel(r'${\lambda},{\mu m}$')

    Temps = [3500, 6000, 10000, 18000, 24000, 30000, 50000]

    a = np.linspace(0, 1, 20)
    q_file = np.genfromtxt(f"results_Q_{ri_file_name}", delimiter=',')
    ax = fig.add_subplot(1, 3, 2)
    plt.title(r'${\overline{Q}_{pr}}$' + f', {ri_file_name}')
    plt.xlabel(r'$a,{\mu m}$')
    for i in range(len(Temps)):
        plt.plot(a, q_file[i], label = f"T={Temps[i]}")
        plt.legend()

    labels = ["A5 I", "A5 V", "G0 I","G0 III", "G0 V", "M0 I", "M0 III", "M0 V"]
        
    b_file = np.genfromtxt(f"results_b_{ri_file_name}", delimiter=',')
    ax = fig.add_subplot(1, 3, 3)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    plt.title(r'${\beta}$' + f', {ri_file_name}')
    plt.xlabel(r'$a,{\mu m}$')
    for i in range(len(labels)):
        plt.plot(a, b_file[i], label = f"T={labels[i]}")
        plt.legend()
    plt.show()