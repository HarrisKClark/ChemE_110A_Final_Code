import math
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['text.usetex'] = False
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['figure.figsize'] = [9, 5]
plt.rcParams.update({'font.size': 20})


def simulate(r, f2a, speed):
    f2a = f2a  # kgfuel/kgair
    t1 = -65  # f
    p1_out = 18  # kpa
    u1 = speed  # mph
    r = r

    fuel_type = "Eth"

    R = 8.314  # j/mol K
    R_kpa = 0.0083145
    CP = (7 / 2) * R

    t1 = (5 / 9) * (t1 - 32) + 273  # conv to k
    p1 = 0.00986923 * p1_out  # conv to atm
    u1 = u1 * 0.44704  # conv to m/s
    v1 = u1  # m^3/s
    # n1 = (p1*V1*1000)/(t1*R_L_atm)
    n1 = (p1_out * v1) / (R_kpa * t1)

    m1 = n1 * 28.96 / 1000  # kg/s
    mm_fuel = 0
    H_r = 0

    if fuel_type == "dodec":
        H_r = -15151000  # j/mole
        mm_fuel = 170.34


    elif fuel_type == "Eth":
        H_r = -1273000  # j/mole
        mm_fuel = 46.068

    m_fuel = f2a * m1
    n_fuel = m_fuel * (1000 / mm_fuel)

    Q = -H_r * n_fuel  # j/mole

    t2 = ((u1 ** 2 * m1) / (2 * CP * n1)) + t1
    p2 = p1 * ((t2 / t1) ** (CP / R))

    p3 = r * p1
    t3_ideal = t2 * ((p3 / p2) ** (R / CP))
    W_comp_ideal = n1 * CP * (t3_ideal - t2)
    n_comp = 0.95 - 0.003 * r
    W_comp_real = W_comp_ideal / n_comp
    t3 = (W_comp_real / (CP * n1)) + t2

    t4 = Q / (CP * n1) + t3
    p4 = p3  # isobaric

    W_turb_real = -W_comp_real
    t5 = W_turb_real / (CP * n1) + t4
    p5 = p4 * ((t5 / t4) ** (CP / R))

    p6 = p1  # same as atmospheric p

    t6 = t5 * ((p6 / p5) ** (R / CP))

    u6 = math.sqrt(((-2 * (n1 * CP) * (t6 - t5)) / (m1)))

    d_EK = (0.5) * (u6 ** 2 - u1 ** 2) * m1
    d_H61 = n1 * CP * (t6 - t1)

    f = m1 * (u6 - u1)

    n_th = d_EK / Q

    n_p = (f * u1) / (d_EK)

    n_final = n_p * n_th

    return n_final


num_val = 50

speed = np.linspace(0, 650, num_val)

n = []

for i in speed:
    n.append(simulate(30, 0.002, i))

plt.plot(speed, n)
plt.xlabel(r"Plane speed $\frac{miles}{hour}$")
plt.ylabel(r"Efficiency")
plt.title("Plane Speed vs Efficiency (Ethanol Fuel)")

plt.show()

