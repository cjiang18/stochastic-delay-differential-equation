
# Introduction  
This is a summer project I did in 2018 summer in Imperial College London. The main task is to simulate the following logistic stochastic delay differential equation. 

$dW=X( \alpha +\beta X_\tau )dt+\sigma X dW$

where $\alpha$ and $\beta$ are parameters of the deterministic delay logistic equation, and $X_\tau(t)=X(t-\tau)$.

We were also ask to simulate the system when a pullback in time is set, in order to find pullback attractors.

I am going to explain some of the algorithms I used in the scripts.

# How Noise Realisation Is Implemented When Pullback Exists

A pullback means we have negative time. For example, if pullback is 200, then the system starts at t = -200.  We want to make sure that the noise realisation should be the same at each time point, so we can adjust the value of pullback and explore pullback attractors.  

Thus, in order to allow pullback to vary without altering the noise realisation, we should not generate noise forward from the very initial time point (e.g. t = -200). Instead, we should generate **two** streams of noise, starting from t = 0. 

Stream A  runs forward from t = 0 to $+\infty$.

Stream B runs backward from t = 0 to $-\infty$.

**NOTICE** : we need two different seeds to generate those streams, otherwise there is a symmetry in noise.

This makes sure that once the random seed is fixed, we can change adjust pullback and the final time as we wish without altering the noise.

# Stochastic Integration Schemes
The philosophy of numeric integration is to discretise time, and use summation to replace integration.
 
Two integration schemes are used for integrations. In most of the scripts, Euler-Maruyama method are used to save computing time. Heun's Method  is only used when stated in the title of the scripts. More sophisticated integration schemes like Runge-Kutta, requires fractional time step, which I found infeasible for stochastic delay differential equation. 
## Euler-Maruyama
The Euler-Maruyama method is basically a stochastic version of the Euler's method for deterministic equation. Under Euler-Maruyama method, our equation becomes

$X(t+dt)-X(t)=X(t)* \left(\alpha +\beta X_\tau(t)\right)*dt+$

## Heun's Method




<!--stackedit_data:
eyJoaXN0b3J5IjpbNzYxMjc2MDAzXX0=
-->