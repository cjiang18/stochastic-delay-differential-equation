import numpy as np
import matplotlib.pyplot as plt

'''''''main'''''''''

'''parameters'''
tfinal = 1000
R = 100  # factor of the unit time size (preferably integer)
tau = 1
sigma = 0
alpha = 2

M = 2  # number of paths
pullback = [0]
seed = 0  # random seed
tinitial = 0
xinitial = np.linspace(1, 5.1, M, endpoint=False)  # initial condition of x
tDelta = 0.0001  # the unit time size (resolution of the graph)
dt = R * tDelta

plt.figure(0)  # constructing figure

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
    x = np.zeros((M, N))
    # v = np.zeros((M, N))  # tangent vector of x

    '''initialisation of t, x and v when pullback exists'''
    if pullback[ii] != 0:
        t1 = np.linspace(tinitial - pullback[ii], tinitial, L, endpoint=False)
        x1 = np.zeros((M, L))
        # v1 = np.zeros((M, L))
        t = np.concatenate((t1, t), axis=0)
        x = np.concatenate((x1, x), axis=-1)
        # del x1,t1
    # v = np.concatenate((v1, v), axis=-1)

    # x[:, 0] = xinitial
    # v[:, 0] = -2

    '''combining t, x and v'''
    tAll = np.concatenate((tHistory, t), axis=0)
    xAll = np.concatenate((xHistory, x), axis=-1)
    # vAll = np.concatenate((vHistory, v), axis=-1)

    del x, t
    del xHistory, tHistory

    '''initialisation of random variable'''

    np.random.seed(seed)
    dW = np.sqrt(tDelta) * np.random.randn(R * N)  # Brownian increments

    '''initialisation of random variable when pullback exists'''
    if pullback[ii] != 0:
        np.random.seed(2 ** 15 + seed)
        temp = np.sqrt(tDelta) * np.random.randn(R * L)
        dW1 = temp[::-1]
        dW = np.concatenate((dW1, dW), axis=-1)
        del dW1, temp

    dW[0] = 0  # condition of Brownian motion

    '''solving logistic equation'''


    for i in range(N + L):
        xAll[:, i + NHistory] = xAll[:, i - 1 + NHistory] - alpha * xAll[:, i] * (
                1 + xAll[:,
                    i - 1 + NHistory]) * dt + sigma * \
                                xAll[:, i - 1 + NHistory] * np.sum(
            dW[R * i:R * (i + 1)])

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
    for i in range(M):
        plt.figure(0)
        plt.plot(tAll, xAll[i, :])

    # plt.figure(ii+1)
    # plt.xlabel('x')
    # plt.ylabel('x-tau')
    # plt.title('phase plot when pullback is ' + str(pullback[ii])+'with tfinal = '+str(tfinal))

plt.figure(0)
plt.xlabel('t')
plt.ylabel('x')
plt.title('trace')
plt.show()
