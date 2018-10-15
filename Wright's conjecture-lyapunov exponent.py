import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
import multiprocessing as mp


def euler(xAll, v, alpha, sigma, dt, R, NHistory, dW, return_dict, N, L):
    '''solving logistic equation'''

    for i in range(N + L):
        temp = xAll[:, (i - 1) % NHistory] - alpha * xAll[:, i % NHistory] * (
                1 + xAll[:, (i - 1) % NHistory]) * dt + sigma * \
               xAll[:, (i - 1) % NHistory] * np.sum(
            dW[R * i:R * (i + 1)], axis=-1)

        '''solving tangent vector'''

        v[:, i % NHistory] = v[:, (i - 1) % NHistory] - alpha * (
                xAll[:, i % NHistory] * v[:, (i - 1) % NHistory] + v[:, i % NHistory] * (
                1 + temp)) * dt + sigma * v[:, (
                                                       i - 1) % NHistory] * np.sum(
            dW[R * i:R * (i + 1)])

        xAll[:, i % NHistory] = temp
    ordered = np.zeros((NHistory, NHistory))
    for j in range(NHistory):
        ordered[:, j] = v[:, (j + N + L ) % NHistory]

    return_dict[0] = ordered


def heun(xAll, v, alpha, sigma, dt, R, NHistory, dW, return_dict, N, L):
    '''solving logistic equation'''

    grad1 = -alpha * xAll[:, 0] * (1 + xAll[:, NHistory - 1])

    temp1 = xAll[:, NHistory - 1] + grad1 * dt + sigma * \
            xAll[:, NHistory - 1] * np.sum(
        dW[0:R], axis=-1)

    grad2 = -alpha * xAll[:, 1] * (1 + temp1)

    temp1 = xAll[:, NHistory - 1] + (grad1 + grad2) * dt / 2 + sigma * \
            xAll[:, NHistory - 1] * np.sum(
        dW[0:R], axis=-1)

    for i in range(N + L):
        xcurrent = temp1
        '''solving logistic equation'''
        grad1 = -alpha * xAll[:, (i + 1) % NHistory] * (1 + xcurrent)

        temp1 = xcurrent + grad1 * dt + sigma * \
                xcurrent * np.sum(
            dW[R * (i + 1):R * (i + 2)], axis=-1)

        grad2 = -alpha * xAll[:, (i + 2) % NHistory] * (1 + temp1)

        temp1 = xcurrent + (grad1 + grad2) * dt / 2 + sigma * \
                xcurrent * np.sum(
            dW[R * (i + 1):R * (i + 2)], axis=-1)

        '''solving tangent vector'''

        grad3 = -alpha * (xAll[:, i % NHistory] * v[:, (i - 1) % NHistory] + v[:, i % NHistory] * (1 + xcurrent))

        temp2 = v[:, (i - 1) % NHistory] + grad3 * dt + sigma * v[:, (i - 1) % NHistory] * np.sum(
            dW[R * i:R * (i + 1)])

        grad4 = -alpha * (xAll[:, (i + 1) % NHistory] * temp2 + v[:, (i + 1) % NHistory] * (1 + temp1))

        v[:, i % NHistory] = v[:, (i - 1) % NHistory] + (grad3 + grad4) * dt / 2 + sigma * v[:,
                                                                                           (i - 1) % NHistory] * np.sum(
            dW[R * i:R * (i + 1)])

        xAll[:, i % NHistory] = xcurrent
        xcurrent = temp1
        if (i + 1) % 10000 == 0:
            print('Heun' + str(100 * (i + 1) / (N + L)) + '% completed')

    ordered = np.zeros((NHistory, NHistory))
    for j in range(NHistory):
        ordered[:, j] = v[:, (j + N + L) % NHistory]
    return_dict[1] = ordered


if __name__ == '__main__':
    '''parameters'''
    tfinal = 1000
    R = 10  # factor of the unit time size (preferably integer)
    tau = 1
    sigma = 0.2
    alpha = 2

    M = 1  # number of paths
    pullback = [0]
    seed = 0  # random seed
    tinitial = 0
    xinitial = np.linspace(1, 5.1, M, endpoint=False)  # initial condition of x
    tDelta = 0.001  # the unit time size (resolution of the graph)
    dt = R * tDelta
    NeedHeun = 1 # 1 for yes, 0 for no

    for ii in range(len(pullback)):

        '''initialisation of history'''
        NHistory = int(tau / dt)
        # tHistory = np.linspace((tinitial - tau - pullback[ii]), (tinitial - pullback[ii]), NHistory, endpoint=False)
        xHistory = np.zeros((M, NHistory))
        # vHistory = np.full((M,NHistory, NHistory), 1)
        for k in range(M):
            xHistory[k, :] = xinitial[k]

        '''parameters'''

        N = int(tfinal / dt)  # number of nodes for non-negative t

        L = int(pullback[ii] / dt)  # number of nodes for negative t

        '''initialisation of t, x and v for non-negative t'''
        # t = np.linspace(tinitial, tfinal, N, endpoint=False)

        # v = np.zeros((M, N))  # tangent vector of x

        '''initialisation of t, x and v when pullback exists'''
        # if pullback[ii] != 0:
        #   t1 = np.linspace(tinitial - pullback[ii], tinitial, L, endpoint=False)

        # v1 = np.zeros((M, L))
        #  t = np.concatenate((t1, t), axis=0)
        # del t1
        # del x1,t1
        # v = np.concatenate((v1, v), axis=-1)

        # v[:, 0] = -2

        '''combining t, x and v'''
        # tAll = np.concatenate((tHistory, t), axis=0)
        xAll = xHistory
        # vAll = np.concatenate((vHistory, v), axis=-1)

        # del t
        del xHistory  # , tHistory

        '''initialisation of random variable'''

        np.random.seed(seed)
        dW = np.sqrt(tDelta) * np.random.randn(R * N)  # Brownian increments

        # all rows of the random variable matrix should be the same

        '''initialisation of random variable when pullback exists'''
        if pullback[ii] != 0:
            np.random.seed(2 ** 15 + seed)
            temp = np.sqrt(tDelta) * np.random.randn(R * L)
            dW1 = temp[::-1]
            dW = np.concatenate((dW1, dW), axis=-1)
            del dW1, temp

        dW[0] = 0  # condition of Brownian motion

        '''initialisation of tangent vector'''

        v = np.eye(NHistory)

        manager = mp.Manager()
        return_dict = manager.dict()

        process1 = mp.Process(target=euler, args=(xAll, v, alpha, sigma, dt, R, NHistory, dW, return_dict, N, L))

        if NeedHeun == 1:
            process2 = mp.Process(target=heun, args=(xAll, v, alpha, sigma, dt, R, NHistory, dW, return_dict, N, L))

        process1.start()

        if NeedHeun == 1:
            process2.start()

        process1.join()

        if NeedHeun == 1:
            process2.join()

        tspan = np.linspace(tfinal - tau, tfinal, NHistory, endpoint=False)
        ''''plotting lyapunov spectrum for euler method'''

        [eigenvalue, eigenvector] = linalg.eig(return_dict[0])
        plt.figure(1)
        lyaE = [np.log(eigenvalue[kk]) / tspan[kk] for kk in range(NHistory)]
        plt.scatter(np.linspace(1, NHistory, NHistory), lyaE, marker='.', label='Euler')
        plt.xlabel('the i th lyapunov exponent')
        plt.ylabel('value of Lyapunov exponent')

        '''plotting lyapunov spectrum for heun method'''
        [eigenvalue, eigenvector] = linalg.eig(return_dict[1])
        plt.figure(1)
        lyaH = [np.log(eigenvalue[kk]) / tspan[kk] for kk in range(NHistory)]
        plt.scatter(np.linspace(1, NHistory, NHistory), lyaH, marker='.', label='Heun')
        plt.legend()

        '''plotting errors'''

        plt.figure()
        plt.scatter(np.linspace(1, NHistory, NHistory), np.array(lyaH) - np.array(lyaE), marker='.')

    plt.show()
