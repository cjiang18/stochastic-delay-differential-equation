# Introduction  
This is a summer project I did in 2018 summer in Imperial College London. The main task is to simulate the following logistic stochastic delay differential equation. 

$dW=X( \alpha +\beta X_\tau )dt+\sigma X dW$

where $\alpha$ and $\beta$ are parameters of the deterministic delay logistic equation, and $X_\tau(t)=X(t-\tau)$.

We were also ask to simulate the system when a pullback in time is set, in order to find pullback attractors.

I am going to explain some of the algorithms I used in the scripts.

# How noise realisation is implemented when pullback exists

A pullback means we have negative time. For example, if pullback is 200, then the system starts at t=-200.  We want to make sure that the noise realisation should be the same at each time point, so we can adjust the value of pullback and explore pullback attractors.  

Thus, in order to allow pullback to vary without altering the noise realisation, we should not generate noise 


<!--stackedit_data:
eyJoaXN0b3J5IjpbMzkyNjc3MzIzLDEzOTkxNzk5ODAsMTAyOD
QyNTEyNiwtMzgzMTM4MzQ0LDYzNjU5MDYzNF19
-->