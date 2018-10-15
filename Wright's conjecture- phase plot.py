import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

'''''''main'''''''''

'''parameters'''
tfinal = 1000
R = 10  # factor of the unit time size (preferably integer)
tau = 1
sigma = 0
alpha = 2

M = 1000  # number of paths
pullback = [0]
seed = 0  # random seed
tinitial = 0
xinitial = np.linspace(1, 5.1, M, endpoint=False)  # initial condition of x
tDelta = 0.001  # the unit time size (resolution of the graph)
dt = R * tDelta
threeD = 1  # 1 for yes, 0 for no

for ii in range(len(pullback)):

    '''initialisation of history'''
    NHistory = int(tau / dt)
    tHistory = np.linspace((tinitial - tau - pullback[ii]), (tinitial - pullback[ii]), NHistory, endpoint=False)
    xHistory = np.zeros((M, NHistory))
    # vHistory = np.full((M, NHistory), 1)
    for k in range(M):
        xHistory[k, :] = xinitial[k]

    '''parameters'''

    N = int(tfinal / dt)  # number of nodes for non-negative t

    L = int(pullback[ii] / dt)  # number of nodes for negative t

    '''initialisation of t, x and v for non-negative t'''
    t = np.linspace(tinitial, tfinal, N, endpoint=False)

    # v = np.zeros((M, N))  # tangent vector of x

    '''initialisation of t, x and v when pullback exists'''
    if pullback[ii] != 0:
        t1 = np.linspace(tinitial - pullback[ii], tinitial, L, endpoint=False)

        # v1 = np.zeros((M, L))
        t = np.concatenate((t1, t), axis=0)

        # del x1,t1
    # v = np.concatenate((v1, v), axis=-1)

    # v[:, 0] = -2

    '''combining t, x and v'''
    tAll = np.concatenate((tHistory, t), axis=0)
    xAll = xHistory
    # vAll = np.concatenate((vHistory, v), axis=-1)

    del t
    del xHistory, tHistory

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

    '''solving logistic equation'''

    for current in range(N + L):
        xAll[:, current % NHistory] = xAll[:, (current - 1) % NHistory] -alpha*xAll[:,current%NHistory]*(1+xAll[:,
                                                                                                           (current-1)%NHistory]) * dt + sigma * \
                                      xAll[:, (current - 1) % NHistory] * np.sum(
            dW[R * current:R * (current + 1)], axis=-1)

        '''solving tangent vector'''

        # v[j, i + 1] = v[j, i] + v[j, i] * (
        #           alpha - beta * xAll[j,i]) * dt - vAll[j,i]) * beta * \
        #           x[j, i] * dt + sigma * v[j, i] * np.sum(dW[j, R * i:R * (i + 1)])
        # vAll = np.concatenate((vHistory, v), axis=-1)

        '''plotting lyapunov exponent '''

        # plt.figure()
        # lya = [np.log(abs(v[j,kk]))/t[kk] for kk in range(len(v[j,:]))]
        # plt.plot(t, lya,'.')
        # plt.xlabel('time')
        # plt.ylabel('Lyapunov exponent')
        # plt.title('xinitial = ' + str(xinitial) + ', pullback = ' + str(pullback)

    '''plotting '''
    if threeD == 0:
        for i in range(M):
            plt.figure(ii + 1)
            plt.plot(xAll[i, current % NHistory], xAll[i, (current+1) % NHistory], '.')
    else:
        fig = plt.figure(ii + 1)
        ax = fig.add_subplot(111, projection='3d')
        for i in range(M):
            ax.scatter(xAll[i, current % NHistory], xAll[i, (current+1) % NHistory],
                       xAll[i, (current + NHistory // 2) % NHistory], marker='.')

    plt.figure(ii + 1)
    plt.xlabel('x')
    plt.ylabel('x-tau')
    plt.title('phase plot when pullback is ' + str(pullback[ii]) + 'with tfinal = ' + str(tfinal))

plt.show()
