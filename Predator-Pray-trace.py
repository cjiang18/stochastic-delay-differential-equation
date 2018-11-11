import numpy as np
import matplotlib.pyplot as plt

'''''''main'''''''''

'''parameters'''
tfinal = 100
R = 1  # factor of the unit time size (preferably integer)
tau = 1.2
sigma = 0
r_1 = 1
a_11 = 1
a_12 = 1
r_2 = 1
a_21 = 3
a_22 = 1

M = 1  # number of different initial conditions, not too large, recommend below 5.
pullback = [0]
seed = 0  # random seed
tinitial = 0
np.random.seed(seed+2)
#xinitial = np.random.uniform(0.1,0.6,M)  # initial condition of x
xinitial=np.full(M,0.6)
tDelta = 0.01  # the unit time size (resolution of the graph)
dt = R * tDelta

plt.figure(0)  # constructing figure

for ii in range(len(pullback)):

    '''initialisation of history'''
    NHistory = int(tau / dt)
    tHistory = np.linspace((tinitial - tau - pullback[ii]), (tinitial - pullback[ii]), NHistory, endpoint=False)
    xHistory = np.zeros((M, 2, NHistory))
    # vHistory = np.full((M, NHistory), 1)
    for k in range(M):
        xHistory[k, 0, 0] = xinitial[k]
        xHistory[k,1,0]=xinitial[k]
    for k in range(1,NHistory):
        xHistory[:,:,k]=xHistory[:,:,k-1]

    '''parameters'''

    N = int(tfinal / dt)  # number of nodes for non-negative t

    L = int(pullback[ii] / dt)  # number of nodes for negative t

    '''initialisation of t, x and v for non-negative t'''
    t = np.linspace(tinitial, tfinal, N, endpoint=False)
    x = np.zeros((M, 2, N))
    # v = np.zeros((M, N))  # tangent vector of x

    '''initialisation of t, x and v when pullback exists'''
    if pullback[ii] != 0:
        t1 = np.linspace(tinitial - pullback[ii], tinitial, L, endpoint=False)
        x1 = np.zeros((M, 2, L))
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
        xAll[:, 0, i + NHistory] = xAll[:, 0, i - 1 + NHistory] + xAll[:, 0, i - 1 + NHistory] * (
                    r_1 - a_11 * xAll[:, 0, i] - a_12 * xAll[:, 1, i - 1 + NHistory]) * dt + sigma * \
                                   xAll[:, 0, i - 1 + NHistory] * np.sum(
            dW[R * i:R * (i + 1)])
        xAll[:, 1, i + NHistory] = xAll[:, 1, i - 1 + NHistory] + xAll[:, 1, i - 1 + NHistory] * (
                -r_2 + a_21 * xAll[:, 0, i - 1 + NHistory] - a_22 * xAll[:, 1, i])*dt + sigma * \
                                   xAll[:, 1, i - 1 + NHistory] * np.sum(
            dW[R * i:R * (i + 1)])

        

    '''plotting '''
    for i in range(M):
        plt.figure(0)
        plt.plot(tAll, xAll[i,0, :])
        plt.figure(1)
        plt.plot(tAll,xAll[i,1,:])
        plt.figure(2)
        plt.plot(xAll[i,0,:],xAll[i,1,:])

    # plt.figure(ii+1)
    # plt.xlabel('x')
    # plt.ylabel('x-tau')
    # plt.title('phase plot when pullback is ' + str(pullback[ii])+'with tfinal = '+str(tfinal))

plt.figure(0)
plt.xlabel('t')
plt.ylabel('x')
plt.title('trace of x')
plt.figure(1)
plt.xlabel('t')
plt.ylabel('y')
plt.title('trace of y')
plt.figure(2)
plt.xlabel('x')
plt.ylabel('y')
plt.title('phase plot')
plt.show()
