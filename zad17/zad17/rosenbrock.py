import numpy as np
import matplotlib.pyplot as plt


def gradient(x):
    grad = np.array([0, 0], dtype = float)
    grad[0] = 400.0 * x[0] ** 3 - 400.0 * x[0] * x[1] + 2.0 * x[0] - 2.0
    grad[1] = 200.0 * (x[1] - x[0] ** 2)
    return grad

def hessian(x, l = 0):
    hess = np.eye(2, dtype = float)
    hess[0, 0] = (1200.0 * x[0] ** 2 - 400.0 * x[1] + 2.0) * (1 + l)
    hess[0, 1] = -400.0 * x[0]
    hess[1, 0] = -400.0 * x[0]
    hess[1, 1] = 200.0 * (1 + l)
    return hess

def rosenbrock(x):
    return (1.0 - x[0])**2 + 100.0 * (x[1] - x[0]**2)**2


testRange = 12
totalSteps = 0
gridSize = 200


for loop in range(testRange):

    randomX = np.random.uniform(-2.0, 2.0)
    randomY = np.random.uniform(-1.0, 3.0) 

    print("Random starting point: (", randomX, ", ", randomY, ")")
    x = np.array([randomX, randomY], dtype = float)
    traceX = np.array([randomX], dtype = float)
    traceY = np.array([randomY], dtype = float)

    dx = 1
    l = 2 ** (-10)
    step = 0
    gradX = gradient(x)

    while dx > 2 ** (-10):
        someX = x - np.linalg.solve(hessian(x, l), gradX)
        if rosenbrock(someX) > rosenbrock(x):
            l = l * 8
        else:
            step += 1
            l = l / 8
            dx = np.linalg.norm(x - someX)
            x = someX
            traceX = np.append(traceX, x[0])
            traceY = np.append(traceY, x[1])
            gradX = gradient(x)


    totalSteps = totalSteps + step
    print("Estimated number of steps: ", round(totalSteps / (loop + 1), 0))


    x, y = np.meshgrid(np.linspace(-2.0, 2.0, gridSize), np.linspace(-1.0, 3.0, gridSize))
    plt.contourf(x, y, rosenbrock([x, y]))
    plt.plot(traceX, traceY, 'ro-')
    plt.title("Searching minimum for Rosenbrock function \nwith Levenberg-Marquardta algoritm (" + str(step) + " steps) from point " + str(round(randomX, 1)))
    plt.xlabel("X")
    plt.ylabel("Y")
    stepNumber = 0
    for x, y in zip(traceX, traceY):
        plt.text(x, y, str(stepNumber), color="black", fontsize=4)
        stepNumber += 1
    plt.savefig("(" + str(round(randomX, 1)) + ", " + str(round(randomY, 1)) + ") in " + str(step) + " steps.png", dpi=500)

    plt.clf()
    plt.cla()
    plt.close()
