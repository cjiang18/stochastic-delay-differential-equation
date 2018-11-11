
# 1. Introduction  
This is a summer project I did in 2018 summer in Imperial College London. The main task is to simulate the following logistic stochastic delay differential equation. 

<img src="/tex/d54e858e85ff6acf67687130083eb7969fb1ead4f46bfbe356eed30dec7b748f.svg?invert_in_darkmode&sanitize=true" align=middle width=2281445pt height=24.65753399999998pt/>

where <img src="/tex/c745b9b57c145ec5577b82542b2df546.svg?invert_in_darkmode&sanitize=true" align=middle width=10.57650494999999pt height=14.15524440000002pt/> and <img src="/tex/8217ed3c32a785f0b5aad4055f432ad8.svg?invert_in_darkmode&sanitize=true" align=middle width=10.16555099999999pt height=22814486pt/> are parameters of the deterministic delay logistic equation, and <img src="/tex/647ddedd0d2f600c40dbbe8108056d5d.svg?invert_in_darkmode&sanitize=true" align=middle width=125.24022225pt height=24.65753399999998pt/>.

We were also ask to simulate the system when a pullback in time is set, in order to find pullback attractors.

I am going to explain some of the algorithms I used in the scripts.

# 2. How Noise Realisation Is Implemented When Pullback Exists

A pullback means we have negative time. For example, if pullback is 200, then the system starts at t = -200.  We want to make sure that the noise realisation should be the same at each time point, so we can adjust the value of pullback and explore pullback attractors.  

Thus, in order to allow pullback to vary without altering the noise realisation, we should not generate noise forward from the very initial time point (e.g. t = -200). Instead, we should generate **two** streams of noise, starting from t = 0. 

Stream A  runs forward from t = 0 to <img src="/tex/701fa44621fd283e3f2c5468958859d8.svg?invert_in_darkmode&sanitize=true" align=middle width=29.223836399999986pt height=19.1781018pt/>.

Stream B runs backward from t = 0 to <img src="/tex/1d5ba78bbbafd3226f371146bc348363.svg?invert_in_darkmode&sanitize=true" align=middle width=29.223836399999986pt height=19.1781018pt/>.

**NOTICE** : we need two different seeds to generate those streams, otherwise there is a symmetry in noise.

This makes sure that once the random seed is fixed, we can change adjust pullback and the final time as we wish without altering the noise.

# 3. Stochastic Integration Schemes
The philosophy of numeric integration is to discretise time, and use summation to replace integration.
 
Two integration schemes are used for integrations. In most of the scripts, Euler-Maruyama method are used to save computing time. Heun's Method  is only used when stated in the title of the scripts. More sophisticated integration schemes like Runge-Kutta, requires fractional time step, which I found infeasible for stochastic delay differential equation. 

## 3.1 Euler-Maruyama

The Euler-Maruyama method is basically a stochastic version of the Euler's method for deterministic equation. Under Euler-Maruyama method, our equation becomes

<img src="/tex/4c6fda2b84bf19d1ef8c7d29f1c50d75.svg?invert_in_darkmode&sanitize=true" align=middle width=488.5447545pt height=24.65753399999998pt/>

<img src="/tex/f11145648cff3ba9c4465b8461448c77.svg?invert_in_darkmode&sanitize=true" align=middle width=127.73402069999999pt height=24.65753399999998pt/> follows a normal distribution with variance <img src="/tex/5a8af6f173febd968ef4c52695efcf85.svg?invert_in_darkmode&sanitize=true" align=middle width=14.492060549999989pt height=22.831056599999986pt/>. Thus, <img src="/tex/f11145648cff3ba9c4465b8461448c77.svg?invert_in_darkmode&sanitize=true" align=middle width=127.73402069999999200b0a90bb002ce5f59368cd4a8fb0b73639592c8dfbbd6bef7fe53edf1872af439d5b12dc6c860995693aa45e4255d1.svg?invert_in_darkmode&sanitize=true" align=middle width=127.7340206378.04251044217.02084683.46465164999999958pt height=24.65753392.831056599999986pt/> is realised by drawing a sample from the

W follows a normal distribution with variance <img src="/tex/5a8af6f173febd968ef4c52695efcf85.svg?invert_in_darkmode&sanitize=true" align=middle width=14.492060549999989pt height=22.831056599999986pt/>. In fact, in my implementation, a smaller times step called **tDelta** is set, and <img src="/tex/e65c57c88cc602403a9760a73adca1ec.svg?invert_in_darkmode&sanitize=true" align=middle width=112.05319124999998pt height=22.831056599999986pt/>, where <img src="/tex/1e438235ef9ec72fc51ac5025516017c.svg?invert_in_darkmode&sanitize=true" align=middle width=12.60847334999999pt height=22.465723500000017pt/> is an integer. Now

<img src="/tex/c5b43f0282cd8e2163c09c50dbc558a43a7611d3205be2e25ef167b61ff5bcf8.svg?invert_in_darkmode&sanitize=true" align=middle width=510.09202365176.58099689999997pt height=29.78995477.8211450000000243pt/>

## 3.2  Heun's Method
i
Heun's Method is supposed to be more accurate than Euler's Method for integrating the deterministic equation, but it is more time- consuming. The scheme for integrating the random variable is the same as Euler-Maruyama, as Heun's Method only improves evaluation of the deterministic gradient. 

For the sake of simplicity, suppose our equation is 

<img src="/tex/331b59c0123e3ea1aa2c0f4710955fc5276736f3fa54ff8824088ff97508ce5636998de7effa663e2e390430a2c4ffe9.svg?invert_in_darkmode&sanitize=true" align=middle width=217.129004.34347945161.7408703499999967pt height=24.65753399999998pt/>

Then, under Heun's method, we first use the same technique as Euler's method to find the value of X(t+dt) by

<img src="/tex/09cf12e61c237b88b551100a84c3023c3670dc0c8303a17ad177aadde8179653.svg?invert_in_darkmode&sanitize=true" align=middle width=315.5992784999999702.8136891999997974pt height=24.65753399999998pt/>

However, this <img src="/tex/2d5e9e9001f4057fdc75304f69d1b973e200e34da364d469d9368cf343f10844.svg?invert_in_darkmode&sanitize=true" align=middle width=68.21345024999998260.2110786pt height=24.65753399999998pt/> is only an intermediate value. The purpose is to use this to evaluate <img src="/tex/18da5344aa0718fc0560cf835cbbb5eab2367977fd59de49ecc94c7e87b2ab75.svg?invert_in_darkmode&sanitize=true" align=middle width=90.763.0993426349999998pt height=24.65753399999998pt/> and then use the average gradient to evaluate <img src="/tex/2d5e9e9001f4057fdc75304f69d1b973.svg?invert_in_darkmode&sanitize=true" align=middle width=68.21345024999998pt height=24.65753399999998pt/> again.

<img src="/tex/f0ab9226dc1b7ec0ca0e2a9cbdc2d5b3.svg?invert_in_darkmode&sanitize=true" align=middle width=446.11408215000006pt height=27.77565449999998pt/>

# 4. Lyapunov spectrum

The Lyapunov specturm is very useful in determining types of attractors. Negative Lyapunov spectra mean **stable** attractors, and positive Lyapunov spectra mean **unstable** attractors, and Lyapunov spectra with both positive and negative values mean **strange** attractors. 

Although the delay funcequation <img src="/tex/057fde3677e10e4628746048e05a0584.svg?invert_in_darkmode&sanitize=true" align=middle width=111.62643029999998pt height=24.65753399999998pt/> is infinite dimensional, we estimate lyapunov spectrum using finite points to approximate. 

For the purpose of more robust computation, we use <img src="/tex/00e90768cd54cf7237b1c509f34fd44fe7370d092fe0ec1bc1e0b94742b0e461.svg?invert_in_darkmode&sanitize=true" align=middle width=127.470687254.183642149999998pt height=24.65753399999998pt/> 
points, which are evenly spaced in the time domain, to represent the function, where <img src="/tex/5a8af6f173febd968ef4c52695efcf85.svg?invert_in_darkmode&sanitize=true" align=middle width=14.492060549999989pt height=22.831056599999986pt/> is set to guarantee that <img src="/tex/c7f4bd27295dced6069c17f5d8b3a740.svg?invert_in_darkmode&sanitize=true" align=middle width=73.79497124999999pt height=22.465723500000017pt/> is an integer. Now the infinite dimensional functions is approximated by a <img src="/tex/c7f4bd27295dced6069c17f5d8b3a740.svg?invert_in_darkmode&sanitize=true" align=middle width=73.79497124999999pt height=22.465723500000017pt/>-dimensional object.  Thus, this <img src="/tex/c7f4bd27295dced6069c17f5d8b3a740.svg?invert_in_darkmode&sanitize=true" align=middle width=73.79497124999999pt height=22.465723500000017pt/>-dimensional object can be represented by a <img src="/tex/c7f4bd27295dced6069c17f5d8b3a740.svg?invert_in_darkmode&sanitize=true" align=middle width=73.79497124999999pt height=22.465723500000017pt/>-dimensional vector, indicating the values the object takes at each time point. 

Since we are dealing with <img src="/tex/c7f4bd27295dced6069c17f5d8b3a740.svg?invert_in_darkmode&sanitize=true" align=middle width=73.79497124999999pt height=22.465723500000017pt/> -dimension, we choose the canonical basis <img src="/tex/94e46ff3acf48b47afcef8b7dc467385.svg?invert_in_darkmode&sanitize=true" align=middle width=153.99570945pt height=24.65753399999998pt/>.  Aligning them together gives us the identity matrix. To estimate the Lyapunov spectrum,  we are going to see how this canonical basis develops when time span is very large, under the linearisation of the of stochastic delay differential equation. Since this system is autonomous, the dynamics of the tangent equation can be written in the form 

<img src="/tex/0e797463df34a1d249fa77d9a0681d94a33a047cc77417ca26fa9fe855632c6a.svg?invert_in_darkmode&sanitize=true" align=middle width=290.760764249999976.4571005pt height=27.6567522pt/>

Thus starting with <img src="/tex/d6ca07c8f420b618e165faa2b0de3548450b5f7853bf1fd5a4036795c5456869b1cd4637e30911a2f47435730f2e24dblinearisation could be 

<img src="/tex/4b23095a3b08e297428f7aa0fb35d446.svg?invert_in_darkmode&sanitize=true" align=middle width=47.3331688451280.01871974106.36128681436299999869857pt height=22.465723500000017pt/>, we have <img src="/tex/9e973a531e7be3193ea1c99e5e69c0efbf09f350ab7974d7d9c8f33db56a9ede.svg?invert_in_darkmode&sanitize=true" align=middle width=76.387033799999988.0536691pt height=27.6567522pt/>. The Lyapunov spectrum is a vector :

<img src="/tex/661a63e5988afc00d5f43115086f5be38a2237fc39c43226f65e31c00b32716e.svg?invert_in_darkmode&sanitize=true" align=middle width=220.086871950000021.16792224999995pt height=27.77565449999998pt/>

where <img src="/tex/b7680b03af3ad50555b9981995ff1ae5.svg?invert_in_darkmode&sanitize=true" align=middle width=40.80723569999999pt height=24.65753399999998pt/> 



is the vector of all real parts of eigenvalues of 







<img src="/tex/02fa8b922302d21626546fa1121e02f1.svg?invert_in_darkmode&sanitize=true" align=middle width=31.96356239999999739fe69d785eb129eeb2746580e64377.svg?invert_in_darkmode&sanitize=true" align=middle width=33.63019769999999ecd73ea956ea7b048e7d8e1736ac6c7fa4b5a7ec9abf2ab526f120b1012df42b9a2d6f5a1121e02f1.svg?invert_in_darkmode&sanitize=true" align=middle width=31.96356239999999b5cefcead9197a1551befc0c6e564d8302ac2a1c53276b77edd633bbc4f8e880b4dbb9c7758e135379ae60.svg?invert_in_darkmode&sanitize=true" align=middle width=221.167922208.6713007499999726.776316763.8349900147.533462981.94564849999959378pt height=24.65753397.77565449999998pt/>. 

Working from first principles, the tangent equation is found to be









<img src="/tex/306dc770a69a0e1b8243fe9cf150bf8365a8c688bf4c94c5ff86a20c2a0f7a0.svg?invert_in_darkmode&sanitize=true" align=middle width=288.68251499999997pt height=2vAll[:, i + NHistory] = vAll[:, i + NHistory - 1] + vAll[:, i + NHistory - 1] * (





74.6567522pt/>

3399999998pt/>

Now, implementing the integration scheme we can find <img src="/tex/02fa8b922302d21626546fa1121e02f1 


vector, indicating the value 
$dW=X( \alpha - +\beta * xAll[i]) * dt - vAll[:, i] * beta * xAll[i + NHistory] * dt + sigma * vAll[:,

i + NHistory - 1] * np.sum(

dW[R * i:R * (i + 1)])
<img src="/tex/787e09f17da1bf1ed2e0885c48397f2c626546fa1121e02f1d039b4dad596f40a869c30fe415810a2c37fae220ca6d72754ed120b2455.svg?invert_in_darkmode&sanitize=true" align=middle width=31.96356239999999288.6825149919.0021866184.213342439.14952689999979584pt height=24.6575332.8310565999999986pt/>, and estimate the Lyapunov Spectrum by choosing a large positive t .

# 5. Invariant Measure

The invariant measure could be understood as a (probability) distribution of particles in the dynamical system when time tends to infinity. Two methods are used to cross-validate the calculation of the invariant measure.

## 5.1 Ulam's Method

The Ulam's method uses Markov Chain.

  Suppose  we have a discrete finite state process, then the invariant measure could be calculated by computing the eigenvalues and eigenvectors of the corresponding Markov matrix. The normalised eigenvector corresponding to eigenvalue 1, is the required invariant measure.

If  we have <img src="/tex/000539e7b945764a553fc8c5317c72b2.svg?invert_in_darkmode&sanitize=true" align=middle width=118.54225185pt height=24.65753399999998pt/>, we could partition the interval <img src="/tex/fe477a2781d275b4481790690fccd15f.svg?invert_in_darkmode&sanitize=true" align=middle width=32.18228144999999pt height=24.65753399999998pt/> into N equally spaced subintervals. Then, considering each subinterval as a single state, we have a N-state Markov process.

And the Markov matrix of this process is estimated by:

<img src="/tex/0e3d73c5cdb36fc4ce1327d3097c921b227b517573842382ed7e695dca301f7420e848d96b3bb30942eeea2b2121d29c11d72c042c143b505dd1fc0d9347a295c7361f1b9be348411a9071d3b078.svg?invert_in_darkmode&sanitize=true" align=middle width=350.470370249527.0710203440.036374292.4629966454.53721569999998pt height=33.20539850.648287999999997pt/>

Of course, the expectation is es

However the process we are dealing with is 

<img src="/tex/d354468c96a6bb401a5480a7c87f2691812fb1b87ae5a39b0a622daf03ce876d.svg?invert_in_darkmode&sanitize=true" align=middle width=290.626248163.32621855pt height=24.65753399999998pt/>

What should we do?

#### Step1:

We have to map <img src="/tex/d64e8fdb66c8e9e5c12535b6b40a9a3b.svg?invert_in_darkmode&sanitize=true" align=middle width=49.751113499999995pt height=24.65753399999998pt/> onto <img src="/tex/f3e711926cecfed3003f9ae341f3d92b.svg?invert_in_darkmode&sanitize=true" align=middle width=11.87217899999999pt height=22.648391699999998pt/>, so we can define a Mmarkov process that we can work with. A natural choice isWe define an equivalence relation on <img src="/tex/d64e8fdb66c8e9e5c12535b6b40a9a3b.svg?invert_in_darkmode&sanitize=true" align=middle width=49.751113499999995pt height=24.65753399999998pt/>:

<img src="/tex/02f3befc8627c8129b0b56479ca81c6d.svg?invert_in_darkmode&sanitize=true" align=middle width=190.13939009999999pt height=24.65753399999998pt/>

Thus <img src="/tex/2579988a07702d44ea90a6cf5484e0d1.svg?invert_in_darkmode&sanitize=true" align=middle width=111.40096934999998pt height=24.65753399999998pt/> i## 5.2 Birkhoff's tThe orequiredm









## 5.2 Birkhoff's tThe mapping which maps <img src="/tex/d64e8fdb66c8e9e5c12535b6b40a9a3b.svg?invert_in_darkmode&sanitize=true" align=middle width=49.751113499999995pt height=24.65753399999998pt/> onto <img src="/tex/f3e711926cecfed3003f9ae341f3d92b.svg?invert_in_darkmode&sanitize=true" align=middle width=11.87217899999999pt height=22.648391699999998pt/>

#### Step2:

We have to restrict <img src="/tex/f3e711926cecfed3003f9ae341f3d92b.svg?invert_in_darkmode&sanitize=true" align=middle width=11.87217899999999pt height=22.648391699999998pt/> to a finite interval <img src="/tex/fe477a2781d275b4481790690fccd15f.svg?invert_in_darkmode&sanitize=true" align=middle width=32.18228144999999pt height=24.65753399999998pt/>. This is totally for the purpose of numerical simulation. We choose:

b = the largest  value a path has ever visited

a = the smallest value a path has ever visited


Now, we have reduced our problem to what we have discussed at the very beginning of this section. So then we estimate the the Markov matrix and compute its normalised eigenvector with  eigenvalue 1.



## 5.2 Birkhoff's Theorem
The Birkhoff's Theorem states that given <img src="/tex/d71b2d59eef000a071a8f0a6b220e28e.svg?invert_in_darkmode&sanitize=true" align=middle width=44.13374129999999pt height=24.65753399999998pt/>



















## 5.2 Birkhoff's Theorem








## 5.2 Birkhoff's Theorem










## 5.2 Birkhoff's Theorem








orem








## 5.2 Birkhoff's Theorem








## 5.2 Birkhoff's Theorem








$X:[a,b]
## 5.2 Birkhoff's Theorem









## 5.2 Birkhoff's Theorem









## 5.2 Birkhoff's Theorem









## 5.1 Ulam's Method
## 5.2 Birkhoff's Theorem









# 6.







.

















X_\tau )dt+\sigma X dWis 


 
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTgyMjM5OTAyNSwtMTUzNzYyNzI5LDg0OD
k1NDA0OCwxMTczNDUwNzYwLDE1NzU4MDMyNzIsODYyNTI1MTE4
LC0yNzQ5NzgwNjYsMTkxMTYzMDk1OCwtMTAxMzgzNzk1MCwtNj
A4ODM1MzQyLC04Njc5NTE2NSwxMzQyNjcxODY0LDI2NTg3NDE0
MCwxNDQ2MjAzNDUxLC02MjE3MDIwMzUsLTI1OTIwODQzMiwtMj
EzMjE2MDM0NSwtMTQwNTA4MzcxMSwtMTM2NzgxNzc3MSwtODAy
NTg1MjcxXX0=
-->