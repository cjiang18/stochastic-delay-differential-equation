

'''this script is used to plot trace of a particle with different pullback and initial values'''

import numpy as np
import matplotlib.pyplot as plt






'''parameters'''
tfinal = 3000
R = 1   # factor of the unit time size (preferably integer)
tau = 7  #delay
sigma = 0 
alpha = 0.2
beta = 0.5
M = 2  # number of paths with same noise realisation but different initial values
pullback = [0] # a list of numbers
seed = 0  # random seed
tinitial = 0 #do not change
xinitial = np.linspace(1,5.1,M,endpoint=False)  # a list of initial conditions of x
tDelta = 0.01  # the unit time size (resolution of the graph)
dt = R * tDelta



for ii in range(len(pullback)):

        '''initialisation of history'''
        NHistory = int(tau / dt)
        tHistory = np.linspace((tinitial - tau - pullback[ii]), (tinitial - pullback[ii]), NHistory, endpoint=False)
        xHistory = np.zeros((M, NHistory))

        for k in range(M):
            xHistory[k,:]=xinitial[k]

        '''parameters'''

        N = int(tfinal / dt)  # number of nodes for non-negative t

        L = int(pullback[ii] / dt)  # number of nodes for negative t

        '''initialisation of t, x and v for non-negative t'''
        t = np.linspace(tinitial, tfinal, N, endpoint=False)
        x = np.zeros((M, N))


        '''initialisation of t, x and v when pullback exists'''
        if pullback[ii] != 0:
            t1 = np.linspace(tinitial - pullback[ii], tinitial, L, endpoint=False)
            x1 = np.zeros((M, L))

            t = np.concatenate((t1, t), axis=0)
            x = np.concatenate((x1, x), axis=-1)



        '''combining t, x and v'''
        tAll = np.concatenate((tHistory, t), axis=0)
        xAll = np.concatenate((xHistory, x), axis=-1)


        del x,t
        del xHistory,tHistory

        '''initialisation of random variable'''
        dW = np.zeros((M, R * N))
        np.random.seed(seed)
        dW[0,:]=np.sqrt(tDelta) * np.random.randn(R * N)# Brownian increments

        for i in range(M-1):
            dW[i+1, :] = dW[i,:]   #all rows of the random variable matrix should be the same

        '''initialisation of random variable when pullback exists'''
        if pullback[ii] != 0:
            dW1 = np.zeros((M, R * L))
            np.random.seed( 2 ** 15 + seed)
            temp = np.sqrt(tDelta) * np.random.randn(R * L)
            dW1[0, :] = temp[::-1]
            for i in range(M-1):
                dW1[i+1,:]=dW1[i,:]
            dW = np.concatenate((dW1, dW), axis=-1)
            del dW1,temp

        dW[:, 0] = 0  # condition of Brownian motion


        '''solving logistic equation'''

        for i in range(N + L ):
            xAll[:, i + NHistory] = xAll[:, i + NHistory - 1] + xAll[:, i + NHistory - 1] * (
                    alpha - beta * xAll[:, i]) * dt + sigma * \
                                    xAll[:, i + NHistory - 1] * np.sum(
                dW[:, R * i:R * (i + 1)], axis=-1)




        '''plotting '''
        for i in range(M):
            plt.figure(0)
            plt.plot(tAll, xAll[i,:])
            plt.figure(ii+1)
            plt.plot(xAll[i,-1],xAll[i,-1-NHistory],'.')






plt.figure(0)
plt.xlabel('t')
plt.ylabel('x')
plt.title('trace')
plt.show()
