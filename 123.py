import matplotlib.pyplot as plt
import math

import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve

# Элементы выборки
T_LIST = [
        151.3466, 242.2403, 191.665, 185.5005, 271.1247,
        139.1302, 224.5614, 223.8082, 156.4592, 188.6266,
        85.7753, 76.2705, 182.2532, 229.4682, 279.1218,
        103.6403, 188.8004, 171.79, 45.7548, 148.7198,
        112.2839, 233.0046, 143.9457, 181.8659, 113.154,
        194.6064, 134.6452, 204.0216, 210.6647, 222.6096,
        168.5761, 88.8053, 127.6561, 269.4633, 109.7671,
        240.9062, 183.5455, 354.223, 86.6773, 167.232,
        96.629, 296.7801, 33.3702, 228.4972, 238.6932,
        253.2475, 89.0315, 159.8398, 134.0248, 234.3971
    ]

# N- число элементов в нашей выборке
N = 50


def program_1():
    # Расчет точечных оценок четырех выборочных начальных моментов.
    # расчет первых четырех начальных моментов
    m1 = sum(T_LIST) / N
    # T-мат ожидание и оно равно м1
    T = m1
    print("T is mat ozhidanie ", round(T, 2))
    print("расчет первых четырех начальных моментов")
    # создаем новый список со всеми элем-ми возведенными в ^2
    m2_list = [x ** 2 for x in T_LIST]
    m2 = sum(m2_list) / N
    # создаем новый список со всеми элем-ми возведенными в ^3
    m3_list = [x ** 3 for x in T_LIST]
    m3 = sum(m3_list) / N
    # создаем новый список со всеми элем-ми возведенными в ^4
    m4_list = [x ** 4 for x in T_LIST]
    m4 = sum(m4_list) / N
    print("m1 is ", round(m1, 2))
    print("m2 is ", round(m2, 2))
    print("m3 is ", round(m3, 2))
    print("m4 is ", round(m4, 2))
    # Строим график эмпирической функции распределения
    n = sum(T_LIST)
    fet = []
    start = 0
    for t in T_LIST:
        fet.append(round((start + t) / n, 2))
        start = start + t
    fet.sort()
    T_LIST.sort()
    # T_LIST.insert(0, 0)
    print(f"fet = {fet}")
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot()
    ax.set_ylabel("Fe")
    ax.set_xlabel("t")
    ax.step(T_LIST, fet, where="post")
    ax.grid()
    plt.title("График эмпирической функции распределения случайной величины")
    plt.show()

    return m1, m2, m3, m4, T, fet


def program_2(m1, m2, m3, m4):
    # 1.2.	Расчет точечных оценок четырех центральных моментов
    # с учетом поправочных коэффициентов для обеспечения
    # несмещенности вычисляемых моментов.

    # вычислены оценки первых четырех центральных моментов
    # по известным связующим соотношениям между ними
    # оценка первых 4 центтральных моментов
    print("оценка первых 4 центральных моментов")
    mu1 = T = m1
    print("mu1 is ", round(mu1, 2))
    mu2 = m2 - m1 ** 2
    print("mu2 is ", round(mu2, 2))
    mu3 = m3 - 3 * m1 * m2 + 2 * m1 ** 3
    print("mu3 is ", round(mu3, 2))
    mu4 = m4 - 4 * m1 * m3 + 6 * (m1 ** 2) * m2 - 3 * m1 ** 4
    print("mu4 is ", round(mu4, 2))

    # Для более удобного восприятия характеристик выборки
    # вводятся нормированные коэффициенты
    # асимметрии и островершинности:
    # Skewness и Excess:
    # Sk - коэф ассиметрии
    # Ex - коэф островершинности

    # ro - дисперсия
    ro = math.sqrt(mu2)
    print("ro is dispersiya ", round(ro, 2))
    # Sk - коэф ассиметрии
    Sk = mu3 / ro ** 3
    print("Sk is kf assimetrii ", round(Sk, 7))
    # Ex - коэф островершинности
    Ex = (mu4 / ro ** 4) - 3
    print("Ex is kf ostroveshinnosti ", round(Ex, 7))
    # Расчет несмещенных центральных моментов
    D_1 = [(x - T) ** 2 for x in T_LIST]
    D = sum(D_1) / N
    print("D is ", round(D, 2))
    print("Расчет несмещенных центральных моментов")
    mu1 = T = m1
    print("mu1 is ", round(mu1, 2))
    mu2H = D
    print("mu2H is ", round(mu2H, 2))
    roH = math.sqrt(D)
    print("roH is ", round(roH, 2))
    SkH = (math.sqrt(N * (N - 1)) / (N - 2)) * Sk
    print("SkH is ", round(SkH, 7))
    ExN = 0.00237103
    print("ExN is kf ostroveshinnosti ", round(ExN, 7))
    ExH = ((N - 1) / ((N - 2) * (N - 3))) * ((N + 1) * Ex + 6)
    print("ExH is kf ostroveshinnosti ", round(ExH, 7))
    return Sk, SkH, Ex, ExH, ro, D


def program_3(m1, m2, m3, m4, T, ro, fet):
    # Критерий стареющего распределения
    print("Критерий стареющего распределения ")
    M1 = m1 / 1
    M2 = m2 / 2
    M3 = m3 / 6
    M4 = m4 / 24
    print("M1 is ", round(M1, 7))
    print("M2 is ", round(M2, 7))
    print("M3 is ", round(M3, 7))
    print("M4 is ", round(M4, 7))

    # Для проверки неравенств (9) рассчитываем значения выражений
    # M3*M1 - M2**2 =
    # M4*M2 - M3**2 =
    temp_1 = M3 * M1 - M2 ** 2
    temp_2 = M4 * M2 - M3 ** 2
    print("temp_1 is ", round(temp_1, 7))
    print("temp_2 is ", round(temp_2, 7))
    print("Расчет параметров предполагаемых")
    print("теоретических распределений методом моментов")
    # Расчет параметров предполагаемых
    # теоретических распределений методом моментов
    # для Гамма-Распределения
    # alfa - параметр мастшаба alfa > 0
    # a - параметр формы a > 0
    v = ro / T
    a = 1 / v ** 2
    print("a is ", round(a, 7))
    alfa = a / T
    print("alfa is ", round(alfa, 7))
    # аналитическое выражение для плотности вероятности:
    # f(t) = ((alfa**a)/G(a)) * (t**a-1) * e**-alfa*t
    D_2 = a / alfa ** 2
    print("D_2 is ", round(D_2, 7))
    # G(a) - гамма функция или эйлеров интеграл 2 рода
    G = ((math.exp(-a)) * (a ** (a - 1 / 2)) *
         (math.sqrt(2 * math.pi)) *
         (1 + (1 / (12 * a)) +
          (1 / (288 * a ** 2)) -
          (139 / (51840 * a ** 3)) -
          (571 / (2488320 * a ** 4))))
    print("G is ", round(G, 7))
    # Плотность вероятности наработки до отказа
    f_ot_t = [((alfa ** a) / G) * (t ** (a - 1)) * math.exp(-alfa * t) for t in T_LIST]
    print("f is ", f_ot_t)
    # Расчет нумерического выражения для плотности вероятности:

    def f(t):
        return ((alfa ** a) / G) * (t ** (a - 1)) * math.exp(-alfa * t)
    F = [quad(f, 0, t)[0] for t in T_LIST]
    print(f"F is {F}")
    # Строим график Аппроксимации эмпирической функции
    F.sort()
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot()
    ax.set_ylabel("Fe, F")
    ax.set_xlabel("t")
    ax.plot(T_LIST, F, label="Эмпирическая функция распределения")
    ax.step(T_LIST, fet, where="post", label="Аппроксимация эмпирической функции")
    ax.grid("Аппроксимация эмпирической функции распределения "
            "гамма распределением с параметрами, найденными методом моментов")
    plt.title("График эмпирической функции распределения случайной величины")
    plt.show()

    Sk_g = 2 / math.sqrt(a)
    print("Sk_g is ", Sk_g)
    Ex_g = 6 / a
    print("Ex_g is ", Ex_g)

    print(f"Sk is {Sk} and Sk_g is {Sk_g}")
    print(f"Ex is {Ex} and Ex_g is {Ex_g}")

    print(f"Sk is {SkH} and Sk_g is {Sk_g}")
    print(f"Ex is {ExH} and Ex_g is {Ex_g}")
    return F


def kolmogorov(Function, function):
    import numpy
    y_t = [abs(F - f) for F, f in zip(Function, function)]
    print(f"y(t) is {y_t}")
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot()
    ax.set_ylabel("y(t)")
    ax.set_xlabel("t")
    plt.title("Изменение модуля «рассогласования» "
              "\nэмпирического и теоретического "
              "\nраспределений во времени")
    ax.plot(T_LIST, y_t)
    plt.show()

    D = numpy.amax(y_t, 0)
    print(f" D is { D}")

    alfa = D * math.sqrt(N)
    print(f"alfa is {alfa}")

def ggg(a, data): #a / data - неизвестные переменные . наши а и б
    a1 = 1 + 1 / a[0]  # 3
    a2 = 1 + 2 / a[0]  # 5

    G = [((math.exp(-a)) * (a ** (a - 1 / 2)) *
          (math.sqrt(2 * math.pi)) *
          (1 + (1 / (12 * a)) +
           (1 / (288 * 1 ** 2)) -
           (139 / (51840 * a ** 3)) -
           (571 / (2488320 * a ** 4)))) for a in [a1, a2]]

    return [data[0] / G[0] - a[1], data[0] / data[1] + 1 - G[1] / G[0]**2] # b , a


def veibula(D, T, fet): #
    x = fsolve(ggg, np.array([1, 1]), [T, D])
    print("___________________________________________")
    print("Veibula")
    a, b = [float(v) for v in x]
    print(f"b = {b}, a = {a}")
    f = [1 - math.exp(-(t/b) ** a) for t in T_LIST]
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot()
    ax.set_ylabel("Fe(t), Fw(t)")
    ax.set_xlabel("t")
    plt.title("Аппроксимация эмпирической функции распределения Вейбулла "
              "\nс параметрами, найденными методом моментов")
    ax.plot(T_LIST, f)
    ax.plot(T_LIST, fet)
    plt.show()
    Sk = 2 / math.sqrt(a)
    print(f"Sk = {Sk}")
    Ex = 6 / float(a)
    print(f"Ex  = {Ex }")
    return f

# ⦁	Расчет параметров предполагаемых теоретических распределений методом максимального правдоподобия


def aaa(a):
    return (-1 + (a - 0.5) / a +
          (-(1 / (12 * a**2)) - 2 / (288 * a ** 3) + (3 * 139) / (51840 * a**4) + (571 * 4) / (2488320 * a**5)) /
          (1 + (1 / (12 * a)) + (1 / (288 * a ** 2)) - (139 / (51840 * a ** 3)) - (571 / (2488320 * a ** 4))))


def gamma_max(T): #⦁	Для гамма-распределения

    a = fsolve(aaa, np.array([1]))
    # print(f"a  = {a}")
    def G(a):
        return (-1 + (a - 0.5) / a +
          (-(1 / (12 * a**2)) - 2 / (288 * a **3) + (3 * 139) / (51840 * a**4) + (571 * 4) / (2488320 * a**5)) /
          (1 + (1 / (12 * a)) + (1 / (288 * 1 ** 2)) - (139 / (51840 * a ** 3)) - (571 / (2488320 * a ** 4))))
    G = [quad(G, 0.5, t)[0] for t in T_LIST]
    print(f"G  = {G}")
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot()
    ax.set_ylabel("Fe(t), Fw(t)")
    ax.set_xlabel("t")
    plt.title("Аппроксимация эмпирической функции распределения Вейбулла "
              "\nс параметрами, найденными методом моментов")
    ax.plot(T_LIST, G)
    plt.show()

    a = 1 / N * math.fsum([math.log(t / T) for t in T_LIST])
    print(f"a  = {a}")
    alfa = a / T
    print(f"alfa  = {alfa}")

    Ga = (-1 + (a - 0.5) / a +
          (-1 / (12 * a**2) - 2 / (288 * a **3) + (3 * 139) / (51840 * a**4) + (571 * 4) / (2488320 * a**5)) /
          (1 + (1 / (12 * a)) + (1 / (288 * a ** 2)) - (139 / (51840 * a ** 3)) - (571 / (2488320 * a ** 4))))
    print(f"Ga  = {Ga}")


if __name__ == "__main__":
    # Методом моентов гамма
    m1, m2, m3, m4, T, fet = program_1()
    Sk, SkH, Ex, ExH, ro, D = program_2(m1, m2, m3, m4,)
    F = program_3(m1, m2, m3, m4, T, ro, fet)
    kolmogorov(fet, F)
    # методом моментв Вейбула
    FW = veibula(D, T, fet)
    kolmogorov(fet, FW)
    gamma_max(T)
