import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg


def pullbacklogistic(pullback, xinitial):
    """"""

    '''parameters'''
    N = int(tfinal / dt)  # number of nodes for non-negative t

    L = int(pullback / dt)  # number of nodes for negative t

    '''initialisation of t, x and v for non-negative t'''
    t = np.linspace(tinitial, tfinal, N, endpoint=False)
    x = np.zeros(N)
    v = np.zeros((NHistory, N))  # tangent vector of x

    '''initialisation of t, x and v when pullback exists'''
    if pullback != 0:
        t1 = np.linspace(tinitial - pullback, tinitial, L, endpoint=False)
        x1 = np.zeros(L)
        v1 = np.zeros((NHistory, L))
        t = np.concatenate((t1, t))
        x = np.concatenate((x1, x))
        v = np.concatenate((v1, v), axis=1)

        del t1,x1,v1

    '''combining t, x and v'''
    tAll = np.concatenate((tHistory, t))
    xAll = np.concatenate((xHistory, x))
    vAll = np.concatenate((vHistory, v), axis=1)

    del tHistory,t
    del xHistory,x
    del vHistory,v

    '''initialisation of random variable'''

    np.random.seed(seed)
    dW = np.sqrt(tDelta) * np.random.randn(R * N)  # Brownian increments

    '''initialisation of random variable when pullback exists'''
    if pullback != 0:
        np.random.seed(2 ** 15 + seed)
        temp = np.sqrt(tDelta) * np.random.randn(R * L)
        dW1 = temp[::-1]
        dW = np.concatenate((dW1, dW))

    dW[0] = 0  # condition of Brownian motion

    '''solving logistic equation'''
    for i in range(N + L):
        xAll[i + NHistory] = xAll[i + NHistory - 1] + xAll[i + NHistory - 1] * (alpha - beta * xAll[i]) * dt + sigma * \
                             xAll[i + NHistory - 1] * np.sum(
            dW[R * i:R * (i + 1)])
        '''solving tangent vector'''

        vAll[:, i + NHistory] = vAll[:, i + NHistory - 1] + vAll[:, i + NHistory - 1] * (
                    alpha - beta * xAll[i]) * dt - vAll[:, i] * beta * xAll[i + NHistory] * dt + sigma * vAll[:,
                                                                                                         i + NHistory - 1] * np.sum(
            dW[R * i:R * (i + 1)])

    '''plotting lyapunov spectrum '''
    [eigenvalue, eigenvector] = linalg.eig(vAll[:, L + int(tObs/dt):L + int(tObs/dt) + NHistory])
    plt.figure()
    lya = [np.log(eigenvalue[kk]) / tAll[kk + L + int(tObs/dt)] for kk in range(NHistory)]
    plt.plot(np.linspace(1, NHistory, NHistory), lya,'.')
    plt.xlabel('the i th lyapunov exponent')
    plt.ylabel('value of Lyapunov exponent')
    plt.title('lyapunov spectrum at t = ' +str(tObs)+'with xinitial = ' + str(xinitial) + ', pullback = ' + str(pullback))



    return [tAll, xAll]

'''This file is used for calculating and plotting trace and lyapunov spectrum of a stochastic logistic differential equation. 
You should not set the time step too large or the integration step too small,as it may lead to a memory insufficiency. Use the file,
comparison between euler and heun-lyapunov spectrum, instead, if you does want to explore lyapunov spectrum. In that case, the trace is not stored.'''
if __name__ == '__main__':

        '''parameters'''
        tfinal = 1000
        R = 1  # factor of the unit time size ( integer)
        tau = 7
        sigma =0
        alpha = 1
        beta =10
        pullback = [0]
        seed = 0  # random seed
        tinitial = 0
        xinitial = [2]  # initial condition of x
        tDelta = 0.01  # the unit time size (resolution of the graph)
        dt = R * tDelta
        tObs = tfinal-1     #the time you want to observe, MUST be smaller than tfinal and non-negative,

        plt.figure(0)  # constructing figure

        for ii in pullback:
            for jj in xinitial:
                '''initialisation of history'''
                NHistory = int(tau / dt)
                tHistory = np.linspace((tinitial - tau - ii), (tinitial - ii), NHistory, endpoint=False)
                xHistory = np.full(NHistory, jj)
                vHistory = np.eye(NHistory)

                [tAll, xAll] = pullbacklogistic(ii, jj)

                '''plotting '''

                plt.figure(0)
                plt.plot(tAll, xAll, label='pullback = ' + str(ii) + ', xinitial=' + str(jj))

        plt.figure(0)
        plt.xlabel('t')
        plt.ylabel('x')
        plt.legend()
        plt.show()
