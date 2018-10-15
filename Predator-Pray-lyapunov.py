import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
import multiprocessing as mp


def euler(xAll, v, r_1, r_2, a_11, a_12, a_21, a_22, sigma, dt, R, NHistory, dW, return_dict, N, L):
    '''solving logistic equation'''

    for i in range(N + L):
        temp = np.zeros((2, 1))
        temp[0] = xAll[0, (i - 1) % NHistory] + xAll[0, (i - 1) % NHistory] * (
                r_1 - a_11 * xAll[0, i % NHistory] - a_12 * xAll[1, (i - 1) % NHistory]) * dt + sigma * \
                  xAll[0, (i - 1) % NHistory] * np.sum(
            dW[R * i:R * (i + 1)])
        temp[1] = xAll[1, (i - 1) % NHistory] + xAll[1, (i - 1) % NHistory] * (
                -r_2 + a_21 * xAll[0, (i - 1) % NHistory] - a_22 * xAll[1, i % NHistory]) * dt + sigma * \
                  xAll[1, (i - 1) % NHistory] * np.sum(
            dW[R * i:R * (i + 1)])

        '''solving tangent vector'''

        v[0, :, i % NHistory] = v[0, :, (i - 1) % NHistory] + (
                    -a_11 * temp[0] * v[0, :, i % NHistory] + v[0, :, (i - 1) % NHistory] * (
                        r_1 - a_11 * xAll[0, i % NHistory] - a_12 * xAll[1, (i - 1) % NHistory])) * dt + sigma * v[0,:, (
                                                                                                                              i - 1) % NHistory] * np.sum(
            dW[R * i:R * (i + 1)])

        v[1, :, i % NHistory] = v[1, :, (i - 1) % NHistory] + (
                -a_12 * temp[1] * v[1, :, i % NHistory] + v[1, :, (i - 1) % NHistory] * (
                -r_2 + a_21 * xAll[0, i % NHistory] - a_22 * xAll[1, (i - 1) % NHistory])) * dt + sigma * v[1,:, (
                                                                                                                      i - 1) % NHistory] * np.sum(
            dW[R * i:R * (i + 1)])



        xAll[0, i % NHistory] = temp[0]
        xAll[1, i % NHistory] = temp[1]

        if (i+1)%10000==0:
            print('Euler'+str(100*(i+1)/(N+L))+'% completed')
    ordered = np.zeros((2,NHistory, NHistory))
    for j in range(NHistory):
        ordered[:,:, j] = v[:, (j + N + L ) % NHistory]

    return_dict[0] = ordered


if __name__ == '__main__':
    '''parameters'''
    tfinal = 50
    R = 100  # factor of the unit time size (preferably integer)
    tau = 1.2
    sigma = 0
    r_1 = 1
    a_11 = 1
    a_12 = 1
    r_2 = 1
    a_21 = 3
    a_22 = 1

    pullback = [0]
    seed = 0  # random seed
    tinitial = 0
    xinitial = 0.6  # initial condition of x
    tDelta = 0.0001  # the unit time size (resolution of the graph)
    dt = R * tDelta
    NeedHeun = 1  # 1 for yes, 0 for no

    for ii in range(len(pullback)):

        '''initialisation of history'''
        NHistory = int(tau / dt)
        # tHistory = np.linspace((tinitial - tau - pullback[ii]), (tinitial - pullback[ii]), NHistory, endpoint=False)
        xHistory = np.zeros((2, NHistory))
        # vHistory = np.full((M,NHistory, NHistory), 1)
        xHistory[0, :] = xinitial
        xHistory[1, :] = xinitial

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

        v = np.zeros((2, NHistory, NHistory))
        v[0, :, :] = np.eye(NHistory)
        v[1, :, :] = np.eye(NHistory)

        manager = mp.Manager()
        return_dict = manager.dict()

        process1 = mp.Process(target=euler, args=(
        xAll, v, r_1, r_2, a_11, a_12, a_21, a_22, sigma, dt, R, NHistory, dW, return_dict, N, L))

        process1.start()

        process1.join()

        tspan = np.linspace(tfinal - tau, tfinal, NHistory, endpoint=False)
        ''''plotting lyapunov spectrum for euler method'''

        [eigenvaluex, eigenvectorx] = linalg.eig(return_dict[0][0,:,:])
        [eigenvaluey,eigenvectory] = linalg.eig(return_dict[0][1,:,:])
        plt.figure(1)
        lyaX = [np.log(eigenvaluex[kk]) / tspan[kk] for kk in range(NHistory)]
        plt.scatter(np.linspace(1, NHistory, NHistory), lyaX, marker='.', label='x')
        plt.xlabel('the i th lyapunov exponent')
        plt.ylabel('value of Lyapunov exponent')
        plt.title('lyapunov spectrum of x')
        plt.figure(2)
        lyaY = [np.log(eigenvaluey[kk]) / tspan[kk] for kk in range(NHistory)]
        plt.scatter(np.linspace(1, NHistory, NHistory), lyaY, marker='.', label='y')
        plt.xlabel('the i th lyapunov exponent')
        plt.ylabel('value of Lyapunov exponent')
        plt.title('lyapunov spectrum of y')

    plt.show()
