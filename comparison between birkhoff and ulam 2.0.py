import numpy as np
import matplotlib.pyplot as plt
import time
from scipy import linalg
import multiprocessing as mp

start = time.time()


def markov(pullback, row0, boxwidth, Nbox, NHistory, tfinal, box, dt, xAll, return_dict):
    # calculating markov matrix

    p = np.zeros((Nbox, Nbox))

    for j in range(Nbox):
        total1 = 0
        for k in range(int((pullback + tfinal) / dt - 1)):

            if box[j] <= xAll[row0, NHistory + k] < box[j + 1]:
                total1 += 1
                for i in range(Nbox):
                    if box[i] <= xAll[row0, NHistory + k + 1] < box[i + 1]:
                        p[i, j] += 1
        if total1 != 0:
            p[:, j] = p[:, j] / total1

        if (j + 1) % 10 == 0:
            print('markov chain ' + str(100 * (j + 1) / Nbox) + "% completed")

    for i in range(Nbox):
        print(np.sum(p[:, i]))

    [evals, invariant] = linalg.eig(p)

    density = abs(invariant[:, evals.argmax()])
    density = density / np.sum(density)
    density = density / boxwidth

    print(invariant[:, evals.argmax()])

    plt.figure()
    plt.plot(range(len(evals)), evals, '.')  # plotting eigenvalues

    print(evals)

    print('markov chain finished ' + str(time.time() - start))

    plt.figure()
    plt.imshow(p)
    plt.colorbar()

    return_dict[0] = density


def counting(pullback, row1, boxwidth, Nbox, NHistory, tfinal, box, dt, xAll, N, L, return_dict):
    # calculating markov matrix

    pvector = np.zeros(Nbox)

    for j in range(Nbox):
        total1 = 0
        for k in range(int((pullback + tfinal) / dt)):

            if box[j] <= xAll[row1, NHistory + k] < box[j + 1]:
                total1 += 1
        pvector[j] = total1
        print(str(100 * (j + 1) / Nbox) + "% boxes counted")
    pvector = (pvector / (N + L)) / boxwidth  # transform probability to probability density

    return_dict[1] = pvector


'''''''main'''''''''
'''This file compares the invariant measures calculated by using two different methods, namely Birkhoff theorem and ulam's method.'''
if __name__ == '__main__':
    '''parameters'''
    tfinal = 10000
    R = 1  # factor of the unit time size ( integer)
    tau = 9
    sigma = 0
    alpha = 0.2
    beta = 0.5
    pullback = [0]
    seed = 0  # random seed
    tinitial = 0
    xinitial = 1  # initial condition of x
    tDelta = 0.01  # the unit time size (resolution of the graph)
    dt = R * tDelta
    tObs = tfinal  # the time you want to observe, MUST be smaller than tfinal and non-negative,
    Nbox = 100
    M = 2  # number of different initial condition
    lowerboundary = 0

    for ii in pullback:

        '''initialisation of history'''
        NHistory = int(tau / dt)
        tHistory = np.linspace((tinitial - tau - ii), (tinitial - ii), NHistory, endpoint=False)
        xHistory = np.zeros((M, NHistory))

        for i in range(M):
            np.random.seed(seed + 2 ** i)

            xHistory[i, :] = np.random.uniform(0, xinitial, NHistory)


        '''parameters'''
        N = int(tfinal / dt)  # number of nodes for non-negative t

        L = int(ii / dt)  # number of nodes for negative t

        '''initialisation of t, x and v for non-negative t'''
        t = np.linspace(tinitial, tfinal, N, endpoint=False)
        x = np.zeros((M, N))


        '''initialisation of t, x and v when pullback exists'''
        if pullback != 0:
            t1 = np.linspace(tinitial - ii, tinitial, L, endpoint=False)
            x1 = np.zeros((M, L))

            t = np.concatenate((t1, t))
            x = np.concatenate((x1, x), axis=-1)


        '''combining t, x and v'''
        tAll = np.concatenate((tHistory, t))
        xAll = np.concatenate((xHistory, x), axis=-1)


        '''initialisation of random variable'''

        dW = np.zeros((2, R * N))
        for i in range(M):
            np.random.seed(seed + i + 100)
            dW[i, :] = np.sqrt(tDelta) * np.random.randn(R * N)  # Brownian increments

        '''initialisation of random variable when pullback exists'''
        if pullback != 0:
            for i in range(M):
                np.random.seed(i + 2 ** 15 + seed)
                temp = np.sqrt(tDelta) * np.random.randn(R * L)
                dW1 = temp[::-1]
                dW[i, :] = np.concatenate((dW1, dW[i, :]))

        dW[:, 0] = 0  # condition of Brownian motion

        '''solving logistic equation'''
        for i in range(N + L):
            xAll[:, i + NHistory] = xAll[:, i + NHistory - 1] + xAll[:, i + NHistory - 1] * (
                    alpha - beta * xAll[:, i]) * dt + sigma * \
                                    xAll[:, i + NHistory - 1] * np.sum(
                dW[:, R * i:R * (i + 1)], axis=-1)

        """calculating invariant measure"""
        print("orbits generated" + str(time.time() - start))  # a flag

        max1 = xAll[0, NHistory + 200:].max()
        max2 = xAll[1, NHistory + 200:].max()

        upperboundary = max(max1, max2)
        box = np.linspace(lowerboundary, upperboundary, Nbox + 1)
        if max1 >= max2:
            row = [0, 1]
        else:
            row = [1, 0]

        boxwidth = (upperboundary - lowerboundary) / Nbox

        manager = mp.Manager()
        return_dict = manager.dict()

        process1 = mp.Process(target=markov,
                              args=(ii, row[0], boxwidth, Nbox, NHistory, tfinal, box, dt, xAll, return_dict))

        process1.start()

        '''applying birkhorff theorem'''

        process2 = mp.Process(target=counting,
                              args=(ii, row[1], boxwidth, Nbox, NHistory, tfinal, box, dt, xAll, N, L, return_dict))

        process2.start()

        process1.join()
        process2.join()

        '''plotting '''

        plt.figure()

        plt.plot(box[:-1], return_dict[0], '.')
        plt.title('invariant measure from markov chain with tfinal = ' + str(tfinal))

        plt.figure()
        plt.plot(box[:-1], return_dict[1], '.')
        plt.title('invariant measure by Birkhoff theorem with tfinal = ' + str(tfinal))

        plt.figure()
        plt.plot(box[:-1], return_dict[0] - return_dict[1], '.')
        plt.title('difference between markov chain and Birkhoff theorem with tfinal = ' + str(tfinal))

    plt.figure(0)
    for i in range(M):
        plt.plot(tAll, xAll[i, :])

    end = time.time()
    print(end - start)

    plt.show()
