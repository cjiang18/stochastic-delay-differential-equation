import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import time
from scipy import linalg
import multiprocessing as mp
from mpl_toolkits.mplot3d import Axes3D

start = time.time()


def counting(pullback,rank, xlowerboundary, ylowerboundary, xboxwidth, yboxwidth, Nbox, NHistory, tfinal, dt, xAll, N, L,
             return_dict):
    # calculating markov matrix

    pvector = np.zeros((Nbox+1, Nbox+1))

    for k in range(int((pullback + tfinal) / dt)):
        xco = int(np.floor((xAll[0, NHistory + k] - xlowerboundary) / xboxwidth))
        yco = int(np.floor((xAll[1, NHistory + k] - ylowerboundary) / yboxwidth))
        pvector[xco, yco] += 1

    pvector = (pvector[:-1,:-1] / (N + L)) / (xboxwidth * yboxwidth)  # transform probability to probability density

    return_dict[rank] = pvector


'''''''main'''''''''

if __name__ == '__main__':
    '''parameters'''
    tfinal = 100
    R = 1  # factor of the unit time size (preferably integer)
    tau = 1.2
    sigma = 0.3
    r_1 = 1
    a_11 = 1
    a_12 = 1
    r_2 = 1
    a_21 = 3
    a_22 = 1

    M = 1  # do not change
    pullback = 0
    seed = 0  # random seed
    tinitial = 0
    np.random.seed(seed + 2)
    # xinitial = np.random.uniform(0.1,0.6,M)  # initial condition of x
    xinitial = 0.6
    tDelta = 0.01  # the unit time size (resolution of the graph)
    dt = R * tDelta
    Nbox = 100

    '''initialisation of history'''
    NHistory = int(tau / dt)
    tHistory = np.linspace((tinitial - tau - pullback), (tinitial - pullback), NHistory, endpoint=False)
    xHistory = np.zeros((M, 2, NHistory))

    for i in range(M):
        np.random.seed(seed + 2 ** i)
        # xHistory[i,:] = np.random.uniform(1,3,NHistory)
        xHistory[i, :, :]  = xinitial

    # vHistory = np.eye(NHistory)

    '''parameters'''
    N = int(tfinal / dt)  # number of nodes for non-negative t

    L = int(pullback / dt)  # number of nodes for negative t

    '''initialisation of t, x and v for non-negative t'''
    t = np.linspace(tinitial, tfinal, N, endpoint=False)
    x = np.zeros((M, 2,N))
    # v = np.zeros((NHistory, N))  # tangent vector of x

    '''initialisation of t, x and v when pullback exists'''
    if pullback != 0:
        t1 = np.linspace(tinitial - pullback, tinitial, L, endpoint=False)
        x1 = np.zeros((M, 2,L))
        # v1 = np.zeros((NHistory, L))
        t = np.concatenate((t1, t))
        x = np.concatenate((x1, x), axis=-1)
        # v = np.concatenate((v1, v), axis=1)
        del t1, x1
    '''combining t, x and v'''
    tAll = np.concatenate((tHistory, t))
    xAll = np.concatenate((xHistory, x), axis=-1)
    # vAll = np.concatenate((vHistory, v), axis=1)
    del t, x
    del tHistory, xHistory
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
        xAll[:, 0, i + NHistory] = xAll[:, 0, i - 1 + NHistory] + xAll[:, 0, i - 1 + NHistory] * (
                r_1 - a_11 * xAll[:, 0, i] - a_12 * xAll[:, 1, i - 1 + NHistory]) * dt + sigma * \
                                   xAll[:, 0, i - 1 + NHistory] * np.sum(
            dW[R * i:R * (i + 1)])
        xAll[:, 1, i + NHistory] = xAll[:, 1, i - 1 + NHistory] + xAll[:, 1, i - 1 + NHistory] * (
                -r_2 + a_21 * xAll[:, 0, i - 1 + NHistory] - a_22 * xAll[:, 1, i]) * dt + sigma * \
                                   xAll[:, 1, i - 1 + NHistory] * np.sum(
            dW[R * i:R * (i + 1)])

    """calculating invariant measure"""
    print("orbits generated" + str(time.time() - start))  # a flag

    xupperboundary = np.max(np.max(xAll[:, 0, NHistory :],axis=-1))
    xlowerboundary = np.min(np.min(xAll[:, 0, NHistory :],axis=-1))

    yupperboundary = np.max(np.max(xAll[:, 1, NHistory :],axis=-1))
    ylowerboundary = np.min(np.min(xAll[:, 1, NHistory :],axis=-1))

    xboxwidth = (xupperboundary - xlowerboundary) / Nbox
    yboxwidth = (yupperboundary - ylowerboundary) / Nbox

    manager = mp.Manager()
    return_dict = manager.dict()

    process = []
    '''applying birkhorff theorem'''
    for i in range(M):
        process.append(mp.Process(target=counting,
                                args=(
                                pullback,i, xlowerboundary, ylowerboundary, xboxwidth, yboxwidth, Nbox, NHistory, tfinal, dt,
                                xAll[i, :, :], N, L, return_dict)))

    for i in range(M):
        process[i].start()

    for i in range(M):
        process[i].join()
    prob=return_dict[0]
    for i in range(1,M):
        prob+=return_dict[i]
    prob = prob/M
    '''plotting '''

    fig=plt.figure()
    ax=fig.gca(projection='3d')
    X=np.linspace(xlowerboundary,xupperboundary,Nbox,endpoint=False)
    Y=np.linspace(ylowerboundary,yupperboundary,Nbox,endpoint=False)
    X,Y=np.meshgrid(X,Y)
    surf=ax.plot_surface(X,Y,prob,cmap=cm.coolwarm)


    fig.colorbar(surf,shrink=0.5, aspect=5)





    end = time.time()
    print(end - start)

    plt.show()
