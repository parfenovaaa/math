# coding=utf-8
import math

# import matplotlib.pyplot as plt

# выборка из 50 элементов
import sympy as sy
from matplotlib import pyplot as plt

# t_list = [
#     151.3466, 242.2403, 191.665, 185.5005, 271.1247,
#     139.1302, 224.5614, 223.8082, 156.4592, 188.6266,
#     85.7753, 76.2705, 182.2532, 229.4682, 279.1218,
#     103.6403, 188.8004, 171.79, 45.7548, 148.7198,
#     112.2839, 233.0046, 143.9457, 181.8659, 113.154,
#     194.6064, 134.6452, 204.0216, 210.6647, 222.6096,
#     168.5761, 88.8053, 127.6561, 269.4633, 109.7671,
#     240.9062, 183.5455, 354.223, 86.6773, 167.232,
#     96.629, 296.7801, 33.3702, 228.4972, 238.6932,
#     253.2475, 89.0315, 159.8398, 134.0248, 234.3971
# ]

t_list = [140, 223, 58, 113, 222, 192, 168,
          225, 182, 239, 149, 53, 126, 66,
          242, 165, 145, 159, 205, 196, 130,
          122, 143, 244, 78, 244, 160, 198,
          175, 76, 162, 147, 211, 225, 203,
          153, 92, 117, 133, 109, 142, 128,
          132, 202, 156, 126, 192, 264, 197, 118]

# N- число элементов в нашей выборке
N = 50
# Расчет точечных оценок четырех выборочных начальных моментов.
# расчет первых четырех начальных моментов
m1 = sum(t_list) / N
# T-мат ожидание и оно равно м1
T = m1
print("T is mat ozhidanie ", round(T, 2))
print("расчет первых четырех начальных моментов")
print("m1 is ", round(m1, 2))
# создаем новый список со всеми элем-ми возведенными в ^2
m2_list = [x ** 2 for x in t_list]
m2 = sum(m2_list) / N
# создаем новый список со всеми элем-ми возведенными в ^3
m3_list = [x ** 3 for x in t_list]
m3 = sum(m3_list) / N
# создаем новый список со всеми элем-ми возведенными в ^4
m4_list = [x ** 4 for x in t_list]
m4 = sum(m4_list) / N
print("m2 is ", round(m2, 2))
print("m3 is ", round(m3, 2))
print("m4 is ", round(m4, 2))

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
# 49 / 2256 = 0.02171986
# 51*Ex = - 15.160515
# (51*Ex+6) = -9.160515
# 49 / 2256 *(51*Ex+6) = 0.00237103
# ExN = 0.00237103
ExN = -0.591176
print("Ex is kf ostroveshinnosti ", round(Ex, 7))
print("ExN is kf ostroveshinnosti ", round(ExN, 7))

# Расчет несмещенных центральных моментов
D_1 = [(x - T) ** 2 for x in t_list]
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
print("ExN is ", round(ExN, 7))

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
# G(a) - гамма функция или эйлеров интеграл 2 рода
# a = 6 / Ex
# print("a is ", round(a, 7))
# alfa = a / T
# print("alfa is ", round(alfa, 7))
v = ro / T
a = 1 / v ** 2
print("a is ", round(a, 7))
alfa = a / T
print("alfa is ", round(alfa, 7))

# аналитическое выражение для плотности вероятности:
# f(t) = ((alfa**a)/G(a)) * (t**a-1) * e**-alfa*t
D_2 = a / alfa ** 2
print("D_2 is ", round(D_2, 7))

G = ((math.exp(-a)) * (a ** (a - 1 / 2)) * (math.sqrt(2 * math.pi)) * (1 + (1 / (12 * a)) +
                                                                       (1 / (288 * a ** 2)) -
                                                                       (139 / (51840 * a ** 3)) -
                                                                       (571 / (2488320 * a ** 4))))
print("G is ", round(G, 7))
# f(t) = ((alfa**a)/G(a)) * (t**a-1) * e**-alfa*t
# f(t) = ((alfa**a)/G(a)) * (t**a-1) * e**-alfa*t
t_list.sort()
f = [((alfa**a) / G) * (t**(a - 1)) * math.exp(-alfa * t) for t in t_list]
print("f is ", f)

x = sy.Symbol("x")
l = [sy.integrate(x1, (x, 0, 2)) for x1 in f]

print(l)
n = sum(t_list)
fet = [0]
start = 0
for t in t_list:
    fet.append(round((start+t)/n, 2))
    start = start + t
t_list.insert(0, 0)
f.insert(0, 0)
t_list.sort()
l.sort()
f.sort()
fig = plt.figure(figsize=(8, 8))
# plt.plot(t_list, f)
ax = fig.add_subplot()
# ax.set_xlim(0, 370)
# ax.set_ylim(0, 1.05)
ax.step(t_list, fet, where="post")
# ax.step(t_list, l, where="post")
ax.grid()
plt.show()

















# a = 0.5
#
# a_1 = 1 + 1/a
# G_1 = 1 + (1 / 12 * a_1**(-1)) + (1 / 288 * a_1**(-2)) - (139 / 51840 * a_1**(-3)) - (571 / 2488320 * a_1**(-4))
# A_1 = math.exp(-a_1) * a_1**(a_1 - 0.5) * math.sqrt(2 * math.pi) * G
#
# a_2 = 1 + 2/a
# G_2 = 1 + (1 / 12 * a_2**(-1)) + (1 / 288 * a_2**(-2)) - (139 / 51840 * a_2**(-3)) - (571 / 2488320 * a_2**(-4))
# A_2 = math.exp(-a_2) * a_2**(a_2 - 0.5) * math.sqrt(2 * math.pi) * G
#
# V = 1.10685 - A_2 / A_1**2
# b = D/T**2
# print(f"b: {b}")