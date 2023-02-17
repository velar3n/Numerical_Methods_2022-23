import matplotlib.pyplot as plt
import numpy as np

X = np.array([-0.875, -0.625, -0.375, -0.125, 0.125, 0.375, 0.625, 0.875])
Y = np.array([0.20712, 0.338624, 0.587156, 0.927536, 0.927536, 0.587156, 0.338624, 0.20712])

x1 = np.arange(-0.875, -0.625, 0.001)
x2 = np.arange(-0.625, -0.375, 0.001)
x3 = np.arange(-0.375, -0.125, 0.001)
x4 = np.arange(-0.125, 0.125, 0.001)
x5 = np.arange(0.125, 0.375, 0.001)
x6 = np.arange(0.375, 0.625, 0.001)
x7 =  np.arange(0.625, 0.875, 0.001)

F1 = 0.20712 + 0.606237*(x1 + 0.875) - 0.96263*(x1 + 0.875)**2 + 2.56701*(x1 + 0.875)**3
F2 = 0.338624 + 0.686457*(x2 + 0.625) + 0.96263*(x2 + 0.625)**2 + 1.0722*(x2 + 0.625)**3
F3 = 0.587156 + 1.36881*(x3 + 0.375) + 1.76678*(x3 + 0.375)**2 + -7.18371*(x3 + 0.375)**3
F4 = 0.927536 + 0.905252*(x4 + 0.125) + -3.62101*(x4 + 0.125)**2
F5 = 0.927536 + -0.905252*(x5 - 0.125) + -3.62101*(x5 - 0.125)**2 + 7.18371*(x5 - 0.125)**3
F6 = 0.587156 + -1.36881*(x6 - 0.375) + 1.76678*(x6 - 0.375)**2 + -1.0722*(x6 - 0.375)**3
F7 = 0.338624 + -0.686457*(x7 - 0.625) + 0.96263*(x7 - 0.625)**2 + -1.28351*(x7 - 0.625)**3

plt.plot(X, Y, 'o', label='data')
plt.title('Natural Cubic Spline Interpolation')
plt.plot(x1, F1, x2, F2, x3, F3, x4, F4, x5, F5, x6, F6, x7, F7)
plt.show()