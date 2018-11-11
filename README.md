
# Introduction  
This is a summer project I did in 2018 summer in Imperial College London. The main task is to simulate the following logistic stochastic delay differential equation. 

<img src="/tex/9fb1ead4f46bfbe356eed30dec7b748f.svg?invert_in_darkmode&sanitize=true" align=middle width=2281445pt height=24.65753399999998pt/>

where <img src="/tex/c745b9b57c145ec5577b82542b2df546.svg?invert_in_darkmode&sanitize=true" align=middle width=10.57650494999999pt height=14.15524440000002pt/> and <img src="/tex/8217ed3c32a785f0b5aad4055f432ad8.svg?invert_in_darkmode&sanitize=true" align=middle width=10.16555099999999pt height=22814486pt/> are parameters of the deterministic delay logistic equation, and <img src="/tex/647ddedd0d2f600c40dbbe8108056d5d.svg?invert_in_darkmode&sanitize=true" align=middle width=125.24022225pt height=24.65753399999998pt/>.

We were also ask to simulate the system when a pullback in time is set, in order to find pullback attractors.

I am going to explain some of the algorithms I used in the scripts.

# How Noise Realisation Is Implemented When Pullback Exists

A pullback means we have negative time. For example, if pullback is 200, then the system starts at t = -200.  We want to make sure that the noise realisation should be the same at each time point, so we can adjust the value of pullback and explore pullback attractors.  

Thus, in order to allow pullback to vary without altering the noise realisation, we should not generate noise forward from the very initial time point (e.g. t = -200). Instead, we should generate **two** streams of noise, starting from t = 0. 

Stream A  runs forward from t = 0 to <img src="/tex/701fa44621fd283e3f2c5468958859d8.svg?invert_in_darkmode&sanitize=true" align=middle width=29.223836399999986pt height=19.1781018pt/>.

Stream B runs backward from t = 0 to <img src="/tex/1d5ba78bbbafd3226f371146bc348363.svg?invert_in_darkmode&sanitize=true" align=middle width=29.223836399999986pt height=19.1781018pt/>.

**NOTICE** : we need two different seeds to generate those streams, otherwise there is a symmetry in noise.

This makes sure that once the random seed is fixed, we can change adjust pullback and the final time as we wish without altering the noise.

# Stochastic Integration Schemes
The philosophy of numeric integration is to discretise time, and use summation to replace integration.
 
Two integration schemes are used for integrations. In most of the scripts, Euler-Maruyama method are used to save computing time. Heun's Method  is only used when stated in the title of the scripts. More sophisticated integration schemes like Runge-Kutta, requires fractional time step, which I found infeasible for stochastic delay differential equation. 
## Euler-Maruyama
The Euler-Maruyama method is basically a stochastic version of the Euler's method for deterministic equation. Under Euler-Maruyama method, our equation becomes

<img src="/tex/4c6fda2b84bf19d1ef8c7d29f1c50d75.svg?invert_in_darkmode&sanitize=true" align=middle width=488.5447545pt height=24.65753399999998pt/>

<img src="/tex/f11145648cff3ba9c4465b8461448c77.svg?invert_in_darkmode&sanitize=true" align=middle width=127.73402069999999pt height=24.65753399999998pt/> follows a normal distribution with variance <img src="/tex/5a8af6f173febd968ef4c52695efcf85.svg?invert_in_darkmode&sanitize=true" align=middle width=14.492060549999989pt height=22.831056599999986pt/>. Thus, <img src="/tex/f11145648cff3ba9c4465b8461448c77.svg?invert_in_darkmode&sanitize=true" align=middle width=127.73402069999999200b0a90bb002ce5f59368cd4a8fb0b73639592c8dfbbd6bef7fe53edf1872af439d5b12dc6c860995693aa45e4255d1.svg?invert_in_darkmode&sanitize=true" align=middle width=127.7340206378.04251044217.02084683.46465164999999958pt height=24.65753392.831056599999986pt/> is realised by drawing a sample from the

W follows a normal distribution with variance <img src="/tex/5a8af6f173febd968ef4c52695efcf85.svg?invert_in_darkmode&sanitize=true" align=middle width=14.492060549999989pt height=22.831056599999986pt/>. In fact, in my implementation, a smaller times step called **tDelta** is set, and <img src="/tex/e65c57c88cc602403a9760a73adca1ec.svg?invert_in_darkmode&sanitize=true" align=middle width=112.05319124999998pt height=22.831056599999986pt/>, where <img src="/tex/1e438235ef9ec72fc51ac5025516017c.svg?invert_in_darkmode&sanitize=true" align=middle width=12.60847334999999pt height=22.465723500000017pt/> is an integer. Now

<img src="/tex/c5b43f0282cd8e2163c09c50dbc558a43a7611d3205be2e25ef167b61ff5bcf8.svg?invert_in_darkmode&sanitize=true" align=middle width=510.09202365176.58099689999997pt height=29.78995477.8211450000000243pt/>
## Heun's Methodi
Heun's Method is supposed to be more accurate than Euler's Method for integrating the deterministic equation, but it is more time- consuming. The scheme for integrating the random variable is the same as Euler-Maruyama, as Heun's Method only improves evaluation of the deterministic gradient. 

For the sake of simplicity, suppose our equation is 

<img src="/tex/331b59c0123e3ea1aa2c0f4710955fc5276736f3fa54ff8824088ff97508ce5636998de7effa663e2e390430a2c4ffe9.svg?invert_in_darkmode&sanitize=true" align=middle width=217.129004.34347945161.7408703499999967pt height=24.65753399999998pt/>

Then, under Heun's method, we first use the same technique as Euler's method to find the value of X(t+dt) by

<img src="/tex/09cf12e61c237b88b551100a84c3023c3670dc0c8303a17ad177aadde8179653.svg?invert_in_darkmode&sanitize=true" align=middle width=315.599278402.813689199999974pt height=24.65753399999998pt/>

However, this <img src="/tex/2d5e9e9001f4057fdc75304f69d1b973e200e34da364d469d9368cf343f10844.svg?invert_in_darkmode&sanitize=true" align=middle width=68.21345024999998260.2110786pt height=24.65753399999998pt/> is only an intermediate value. The purpose is to use this to evaluate <img src="/tex/18da5344aa0718fc0560cf835cbbb5eab2367977fd59de49ecc94c7e87b2ab75.svg?invert_in_darkmode&sanitize=true" align=middle width=90.763.0993426349999998pt height=24.65753399999998pt/> and then use the average gradient to evaluate <img src="/tex/2d5e9e9001f4057fdc75304f69d1b973.svg?invert_in_darkmode&sanitize=true" align=middle width=68.21345024999998pt height=24.65753399999998pt/> again.

<img src="/tex/f0ab9226dc1b7ec0ca0e2a9cbdc2d5b3.svg?invert_in_darkmode&sanitize=true" align=middle width=446.11408215000006pt height=27.77565449999998pt/>

# Lyapunov spectrum
The Lyapunov specturm is very useful in determining types of attractors. Negative Lyapunov spectra mean **stable** attractors, and positive Lyapunov spectra mean **unstable** attractors, and Lyapunov spectra with both positive and negative values mean **strange** attractors. 

Although the delay funcequation <img src="/tex/057fde3677e10e4628746048e05a0584.svg?invert_in_darkmode&sanitize=true" align=middle width=111.62643029999998pt height=24.65753399999998pt/> is infinite dimensional, we estimate lyapunov spectrum using finite points to approximate. 

For the purpose of more robust computation, we use <img src="/tex/00e90768cd54cf7237b1c509f34fd44fe7370d092fe0ec1bc1e0b94742b0e461.svg?invert_in_darkmode&sanitize=true" align=middle width=127.470687254.183642149999998pt height=24.65753399999998pt/> 
points, which are evenly spaced in the time domain, to represent the function, where <img src="/tex/5a8af6f173febd968ef4c52695efcf85.svg?invert_in_darkmode&sanitize=true" align=middle width=14.492060549999989pt height=22.831056599999986pt/> is set to guarantee that <img src="/tex/c7f4bd27295dced6069c17f5d8b3a740.svg?invert_in_darkmode&sanitize=true" align=middle width=73.79497124999999pt height=22.465723500000017pt/> is an integer. Now the infinite dimensional functions is approximated by a <img src="/tex/c7f4bd27295dced6069c17f5d8b3a740.svg?invert_in_darkmode&sanitize=true" align=middle width=73.79497124999999pt height=22.465723500000017pt/>-dimensional object.  Thus, this <img src="/tex/c7f4bd27295dced6069c17f5d8b3a740.svg?invert_in_darkmode&sanitize=true" align=middle width=73.79497124999999pt height=22.465723500000017pt/>-dimensional object can be represented by a <img src="/tex/c7f4bd27295dced6069c17f5d8b3a740.svg?invert_in_darkmode&sanitize=true" align=middle width=73.79497124999999pt height=22.465723500000017pt/>-dimensional vector, indicating the values the object takes at each time point. 

Since we are dealing with <img src="/tex/c7f4bd27295dced6069c17f5d8b3a740.svg?invert_in_darkmode&sanitize=true" align=middle width=73.79497124999999pt height=22.465723500000017pt/> -dimension, we 
vector, indicating the value 
$dW=X( \alpha +\beta X_\tau )dt+\sigma X dWis 


 
<!--stackedit_data:
eyJoaXN0b3J5IjpbMTE0MTE0MDQ0MCwtMTU5MDY3OTEsLTI1OT
IwODQzMiwtMjEzMjE2MDM0NSwtMTQwNTA4MzcxMSwtMTM2Nzgx
Nzc3MSwtODAyNTg1MjcxLDQ3MzM3MDA4MV19
-->