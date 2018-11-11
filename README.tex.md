
# 1. Introduction  
This is a summer project I did in 2018 summer in Imperial College London. The main task is to simulate the following logistic stochastic delay differential equation. 

$dX=X( \alpha -\beta X_\tau )dt+\sigma X dW$

where $\alpha$ and $\beta$ are parameters of the deterministic delay logistic equation, and $X_\tau(t)=X(t-\tau)$.

We were also ask to simulate the system when a pullback in time is set, in order to find pullback attractors.

I am going to explain some of the algorithms I used in the scripts.

# 2. How Noise Realisation Is Implemented When Pullback Exists

A pullback means we have negative time. For example, if pullback is 200, then the system starts at t = -200.  We want to make sure that the noise realisation should be the same at each time point, so we can adjust the value of pullback and explore pullback attractors.  

Thus, in order to allow pullback to vary without altering the noise realisation, we should not generate noise forward from the very initial time point (e.g. t = -200). Instead, we should generate **two** streams of noise, starting from t = 0. 

Stream A  runs forward from t = 0 to $+\infty$.

Stream B runs backward from t = 0 to $-\infty$.

**NOTICE** : we need two different seeds to generate those streams, otherwise there is a symmetry in noise.

This makes sure that once the random seed is fixed, we can change adjust pullback and the final time as we wish without altering the noise.

# 3. Stochastic Integration Schemes
The philosophy of numeric integration is to discretise time, and use summation to replace integration.
 
Two integration schemes are used for integrations. In most of the scripts, Euler-Maruyama method are used to save computing time. Heun's Method  is only used when stated in the title of the scripts. More sophisticated integration schemes like Runge-Kutta, requires fractional time step, which I found infeasible for stochastic delay differential equation. 
## 3.1 Euler-Maruyama
The Euler-Maruyama method is basically a stochastic version of the Euler's method for deterministic equation. Under Euler-Maruyama method, our equation becomes

$X(t+dt)-X(t)=X(t)\left[\alpha +\beta X_\tau(t)\right]dt+\sigma X(t)[W(t+dt)-W(t)]$

$W(t+dt)-W(t)$ follows a normal distribution with variance $dt$. Thus, $W(t+dt)-W(t)$ is realised by drawing a sample from the normal distribution with variance $dt$. In fact, in my implementation, a smaller times step called **tDelta** is set, and $dt=R*tDelta$, where $R$ is an integer. Now

$W(t+dt)-W(t)=\displaystyle\Sigma_{i=0}^{R-1}[W(t+(i+1)*tDelta)-W(t+i*tDelta)]$
## 3.2  Heun's Method
Heun's Method is supposed to be more accurate than Euler's Method for integrating the deterministic equation, but it is more time- consuming. The scheme for integrating the random variable is the same as Euler-Maruyama, as Heun's Method only improves evaluation of the deterministic gradient. 

For the sake of simplicity, suppose our equation is 

$dX=\phi(X(t)) dt+\theta(X(t))dW$

Then, under Heun's method, we first use the same technique as Euler's method to find the value of X(t+dt) by

$X(t+dt)=X(t)+\phi(X(t)) dt+\theta(X(t))dW$

However, this $X(t+dt)$ is only an intermediate value. The purpose is to use this to evaluate $\phi(X(t+dt))$ and then use the average gradient to evaluate $X(t+dt)$ again.

$X(t+dt)=X(t)+\frac{1}{2}[\phi(X(t))+\phi(X(t+dt))] dt+\theta(X(t))dW$

# 4. Lyapunov spectrum
The Lyapunov specturm is very useful in determining types of attractors. Negative Lyapunov spectra mean **stable** attractors, and positive Lyapunov spectra mean **unstable** attractors, and Lyapunov spectra with both positive and negative values mean **strange** attractors. 

Although the delay function $X:[-\tau,0]\rightarrow\mathbb{R}$ is infinite dimensional, we estimate lyapunov spectrum using finite points to approximate. 

For the purpose of more robust computation, we use $NHistory=\tau/dt$ points, which are evenly spaced in the time domain, to represent the function, where $dt$ is set to guarantee that $NHistory$ is an integer. Now the infinite dimensional functions is approximated by a $NHistory$-dimensional object.  Thus, this $NHistory$-dimensional object can be represented by a $NHistory$-dimensional vector, indicating the values the object takes at each time point. 

Since we are dealing with $NHistory$ -dimension, we choose the canonical basis $\{e_1,e_2,\dots, e_{Nhistory}\}$.  Aligning them together gives us the identity matrix. To estimate the Lyapunov spectrum,  we are going to see how this canonical basis develops when time span is very large, under the linearisation of the of stochastic delay differential equation. Since this system is autonomous, the dynamics of the tangent equation can be written in the form 

$V(t)=e^{tA}V_0,\quad A\in\mathbb{R}^{NHistory\times NHistory}$

Thus starting with $V_o=I$, we have $V(t)=e^{tA}$. The Lyapunov spectrum is a vector :

$L=\lim(t\rightarrow \infty)\{\frac{1}{t}\log ( \lambda_V(t))\}$

where $\lambda_X(t)$ is the vector of all real parts of eigenvalues of $V(t)$. 

Working from first principles, the tangent equation is found to be

$dV=[V(\alpha-\beta X_\tau)-\beta V_\tau X] dt+\sigma V dW$

Now, implementing the integration scheme we can find $V(t)$, and estimate the Lyapunov Spectrum by choosing a large positive t .

# 5. Invariant Measure
The invariant measure could be understood as a (probability) distribution of particles in the dynamical system when time tends to infinity. Two methods are used to cross-validate the calculation of the invariant measure.
## 5.1 Ulam's Method
The Ulam's method uses Markov Chain. 
## 5.2 Birkhoff's Theorem











 
<!--stackedit_data:
eyJoaXN0b3J5IjpbODYyNTI1MTE4LC0yNzQ5NzgwNjYsMTkxMT
YzMDk1OCwtMTAxMzgzNzk1MCwtNjA4ODM1MzQyLC04Njc5NTE2
NSwxMzQyNjcxODY0LDI2NTg3NDE0MCwxNDQ2MjAzNDUxLC02Mj
E3MDIwMzUsLTI1OTIwODQzMiwtMjEzMjE2MDM0NSwtMTQwNTA4
MzcxMSwtMTM2NzgxNzc3MSwtODAyNTg1MjcxLDQ3MzM3MDA4MV
19
-->