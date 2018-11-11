
# Introduction  
This is a summer project I did in 2018 summer in Imperial College London. The main task is to simulate the following logistic stochastic delay differential equation. 

<img src="/tex/d54e858e85ff6acf67687130083eb796.svg?invert_in_darkmode&sanitize=true" align=middle width=221.58144525pt height=24.65753399999998pt/>

where <img src="/tex/c745b9b57c145ec5577b82542b2df546.svg?invert_in_darkmode&sanitize=true" align=middle width=10.57650494999999pt height=14.15524440000002pt/> and <img src="/tex/8217ed3c32a785f0b5aad4055f432ad8.svg?invert_in_darkmode&sanitize=true" align=middle width=10.16555099999999pt height=22.831056599999986pt/> are parameters of the deterministic delay logistic equation, and <img src="/tex/647ddedd0d2f600c40dbbe8108056d5d.svg?invert_in_darkmode&sanitize=true" align=middle width=125.24022225pt height=24.65753399999998pt/>.

We were also askd to simulate the system when a pullback in time is set, in order to find pullback attractors.

I am going to explain some of the algorithms I used in the scripts.

# How Noise Realisation Is Implemented When Pullback Exists

A pullback means we have  t taking negative values. For example, if pullback is 200, then the system starts at t = -200.  We want to make sure that the noise realisation should be the same at each time point, so we can adjust the value of pullback and explore pullback attractors.  

Thus, in order to allow pullback to vary without altering the noise realisation, we should not generate noise forward from the very initial time point (e.g. t = -200). Instead, we should generate **two** streams of noise, starting from t = 0. 

Stream A  runs forward from t = 0 to <img src="/tex/701fa44621fd283e3f2c5468958859d8.svg?invert_in_darkmode&sanitize=true" align=middle width=29.223836399999986pt height=19.1781018pt/>.

Stream B runs backward from t = 0 to <img src="/tex/1d5ba78bbbafd3226f371146bc348363.svg?invert_in_darkmode&sanitize=true" align=middle width=29.223836399999986pt height=19.1781018pt/>.

**NOTICE** : we need two different seeds to generate those streams, otherwise there is a symmetry in noise.

This makes sure that once the random seed is fixed, we can adjust pullback and the final time as we wish without altering the noise.

# Stochastic Integration Schemes
The philosophy of numeric integration is to discretise time, and use summation to replace integration.
 
Two integration schemes are used for integrations. In most of the scripts, Euler-Maruyama method are used to save computing time. Heun's Method  is only used when stated in the title of the scripts. More sophisticated integration schemes like Runge-Kutta, requires fractional time step, which I found infeasible for stochastic delay differential equation. 

##  Euler-Maruyama

The Euler-Maruyama method is basically a stochastic version of the Euler's method for deterministic equation. Under Euler-Maruyama method, our equation becomes

<img src="/tex/4c6fda2b84bf19d1ef8c7d29f1c50d75.svg?invert_in_darkmode&sanitize=true" align=middle width=488.5447545pt height=24.65753399999998pt/>

<img src="/tex/f11145648cff3ba9c4465b8461448c77.svg?invert_in_darkmode&sanitize=true" align=middle width=127.73402069999999pt height=24.65753399999998pt/> follows a normal distribution with variance <img src="/tex/5a8af6f173febd968ef4c52695efcf85.svg?invert_in_darkmode&sanitize=true" align=middle width=14.492060549999989pt height=22.831056599999986pt/>. Thus, <img src="/tex/f11145648cff3ba9c4465b8461448c77.svg?invert_in_darkmode&sanitize=true" align=middle width=127.73402069999999pt height=24.65753399999998pt/> is realised by drawing a sample from the normal distribution with variance <img src="/tex/5a8af6f173febd968ef4c52695efcf85.svg?invert_in_darkmode&sanitize=true" align=middle width=14.492060549999989pt height=22.831056599999986pt/>. In fact, in my implementation, a smaller times step called **tDelta** is set, and <img src="/tex/e65c57c88cc602403a9760a73adca1ec.svg?invert_in_darkmode&sanitize=true" align=middle width=112.05319124999998pt height=22.831056599999986pt/>, where <img src="/tex/1e438235ef9ec72fc51ac5025516017c.svg?invert_in_darkmode&sanitize=true" align=middle width=12.60847334999999pt height=22.465723500000017pt/> is an integer. Now

<img src="/tex/c5b43f0282cd8e2163c09c50dbc558a4.svg?invert_in_darkmode&sanitize=true" align=middle width=510.09202365pt height=29.789954700000024pt/>

##   Heun's Method

Heun's Method is supposed to be more accurate than Euler's Method for integrating the deterministic equation, but it is more time- consuming. The scheme for integrating the random variable is the same as Euler-Maruyama, as Heun's Method only improves evaluation of the deterministic gradient. 

For the sake of simplicity, suppose our equation is 

<img src="/tex/331b59c0123e3ea1aa2c0f4710955fc5.svg?invert_in_darkmode&sanitize=true" align=middle width=217.12907039999996pt height=24.65753399999998pt/>

Then, under Heun's method, we first use the same technique as Euler's method to find the value of X(t+dt) by

<img src="/tex/09cf12e61c237b88b551100a84c3023c.svg?invert_in_darkmode&sanitize=true" align=middle width=315.59927849999997pt height=24.65753399999998pt/>

However, this <img src="/tex/2d5e9e9001f4057fdc75304f69d1b973.svg?invert_in_darkmode&sanitize=true" align=middle width=68.21345024999998pt height=24.65753399999998pt/> is only an intermediate value. The purpose is to use this to evaluate <img src="/tex/18da5344aa0718fc0560cf835cbbb5ea.svg?invert_in_darkmode&sanitize=true" align=middle width=90.79342634999999pt height=24.65753399999998pt/> and then use the average gradient to evaluate <img src="/tex/2d5e9e9001f4057fdc75304f69d1b973.svg?invert_in_darkmode&sanitize=true" align=middle width=68.21345024999998pt height=24.65753399999998pt/> again.

<img src="/tex/f0ab9226dc1b7ec0ca0e2a9cbdc2d5b3.svg?invert_in_darkmode&sanitize=true" align=middle width=446.11408215000006pt height=27.77565449999998pt/>

# Lyapunov spectrum

The Lyapunov specturm is very useful in determining types of attractors. Negative Lyapunov spectra mean **stable** attractors, and positive Lyapunov spectra mean **unstable** attractors, and Lyapunov spectra with both positive and negative values mean **strange** attractors. 

Although the delay function <img src="/tex/057fde3677e10e4628746048e05a0584.svg?invert_in_darkmode&sanitize=true" align=middle width=111.62643029999998pt height=24.65753399999998pt/> is infinite dimensional, we estimate lyapunov spectrum using finite points to approximate. 

For the purpose of more robust computation, we use <img src="/tex/00e90768cd54cf7237b1c509f34fd44f.svg?invert_in_darkmode&sanitize=true" align=middle width=127.47068729999998pt height=24.65753399999998pt/> points, which are evenly spaced in the time domain, to represent the function, where <img src="/tex/5a8af6f173febd968ef4c52695efcf85.svg?invert_in_darkmode&sanitize=true" align=middle width=14.492060549999989pt height=22.831056599999986pt/> is set to guarantee that <img src="/tex/c7f4bd27295dced6069c17f5d8b3a740.svg?invert_in_darkmode&sanitize=true" align=middle width=73.79497124999999pt height=22.465723500000017pt/> is an integer. Now the infinite dimensional functions is approximated by a <img src="/tex/c7f4bd27295dced6069c17f5d8b3a740.svg?invert_in_darkmode&sanitize=true" align=middle width=73.79497124999999pt height=22.465723500000017pt/>-dimensional object.  Thus, this <img src="/tex/c7f4bd27295dced6069c17f5d8b3a740.svg?invert_in_darkmode&sanitize=true" align=middle width=73.79497124999999pt height=22.465723500000017pt/>-dimensional object can be represented by a <img src="/tex/c7f4bd27295dced6069c17f5d8b3a740.svg?invert_in_darkmode&sanitize=true" align=middle width=73.79497124999999pt height=22.465723500000017pt/>-dimensional vector, indicating the values the object takes at each time point. 

Since we are dealing with <img src="/tex/c7f4bd27295dced6069c17f5d8b3a740.svg?invert_in_darkmode&sanitize=true" align=middle width=73.79497124999999pt height=22.465723500000017pt/> -dimension, we choose the canonical basis <img src="/tex/94e46ff3acf48b47afcef8b7dc467385.svg?invert_in_darkmode&sanitize=true" align=middle width=153.99570945pt height=24.65753399999998pt/>.  Aligning them together gives us the identity matrix. To estimate the Lyapunov spectrum,  we are going to see how this canonical basis develops when time span is very large, under the linearisation of the of stochastic delay differential equation. Since this system is autonomous, the dynamics of the tangent equation can be written in the form 

<img src="/tex/0e797463df34a1d249fa77d9a0681d94.svg?invert_in_darkmode&sanitize=true" align=middle width=290.76076424999997pt height=27.6567522pt/>

Thus starting with <img src="/tex/d6ca07c8f420b618e165faa2b0de3548.svg?invert_in_darkmode&sanitize=true" align=middle width=47.333168849999986pt height=22.465723500000017pt/>, we have <img src="/tex/9e973a531e7be3193ea1c99e5e69c0ef.svg?invert_in_darkmode&sanitize=true" align=middle width=76.38703379999998pt height=27.6567522pt/>. The Lyapunov spectrum is a vector :

<img src="/tex/661a63e5988afc00d5f43115086f5be3.svg?invert_in_darkmode&sanitize=true" align=middle width=220.08687195000002pt height=27.77565449999998pt/>

where <img src="/tex/b7680b03af3ad50555b9981995ff1ae5.svg?invert_in_darkmode&sanitize=true" align=middle width=40.80723569999999pt height=24.65753399999998pt/> is the vector of all real parts of eigenvalues of <img src="/tex/02fa8b922302d21626546fa1121e02f1.svg?invert_in_darkmode&sanitize=true" align=middle width=31.96356239999999pt height=24.65753399999998pt/>. 

Working from first principles, the tangent equation is found to be

<img src="/tex/306dc770a69a0e1b8243fe9cf150bf83.svg?invert_in_darkmode&sanitize=true" align=middle width=288.68251499999997pt height=24.65753399999998pt/>

Now, implementing the integration scheme we can find <img src="/tex/02fa8b922302d21626546fa1121e02f1.svg?invert_in_darkmode&sanitize=true" align=middle width=31.96356239999999pt height=24.65753399999998pt/>, and estimate the Lyapunov Spectrum by choosing a large positive t .

#  Invariant Measure

The invariant measure could be understood as a (probability) distribution of particles in the dynamical system when time tends to infinity. Two methods are used to cross-validate the calculation of the invariant measure.

##  Ulam's Method

The Ulam's method uses Markov Chain.

 Suppose  we have a discrete finite state process, then the invariant measure could be calculated by computing the eigenvalues and eigenvectors of the corresponding Markov matrix. The normalised eigenvector corresponding to eigenvalue 1, is the required invariant measure.

If  we have <img src="/tex/000539e7b945764a553fc8c5317c72b2.svg?invert_in_darkmode&sanitize=true" align=middle width=118.54225185pt height=24.65753399999998pt/>, we could partition the interval <img src="/tex/fe477a2781d275b4481790690fccd15f.svg?invert_in_darkmode&sanitize=true" align=middle width=32.18228144999999pt height=24.65753399999998pt/> into N equally spaced subintervals. Then, considering each subinterval as a single state, we have a N-state Markov process.

And the Markov matrix of this process is estimated by:

<img src="/tex/0e3d73c5cdb36fc4ce1327d3097c921b.svg?invert_in_darkmode&sanitize=true" align=middle width=350.4703702499999pt height=33.20539859999999pt/>



However the process we are dealing with is 

<img src="/tex/d354468c96a6bb401a5480a7c87f2691.svg?invert_in_darkmode&sanitize=true" align=middle width=290.626248pt height=24.65753399999998pt/>

What should we do?

#### Step1:

We have to map <img src="/tex/d64e8fdb66c8e9e5c12535b6b40a9a3b.svg?invert_in_darkmode&sanitize=true" align=middle width=49.751113499999995pt height=24.65753399999998pt/> onto <img src="/tex/f3e711926cecfed3003f9ae341f3d92b.svg?invert_in_darkmode&sanitize=true" align=middle width=11.87217899999999pt height=22.648391699999998pt/>, so we can define a Markov process that we can work with. A natural choice is an equivalence relation on <img src="/tex/d64e8fdb66c8e9e5c12535b6b40a9a3b.svg?invert_in_darkmode&sanitize=true" align=middle width=49.751113499999995pt height=24.65753399999998pt/>:

<img src="/tex/02f3befc8627c8129b0b56479ca81c6d.svg?invert_in_darkmode&sanitize=true" align=middle width=190.13939009999999pt height=24.65753399999998pt/>

Thus <img src="/tex/2579988a07702d44ea90a6cf5484e0d1.svg?invert_in_darkmode&sanitize=true" align=middle width=111.40096934999998pt height=24.65753399999998pt/> is the required mapping which maps <img src="/tex/d64e8fdb66c8e9e5c12535b6b40a9a3b.svg?invert_in_darkmode&sanitize=true" align=middle width=49.751113499999995pt height=24.65753399999998pt/> onto <img src="/tex/f3e711926cecfed3003f9ae341f3d92b.svg?invert_in_darkmode&sanitize=true" align=middle width=11.87217899999999pt height=22.648391699999998pt/>

#### Step2:

We have to restrict <img src="/tex/f3e711926cecfed3003f9ae341f3d92b.svg?invert_in_darkmode&sanitize=true" align=middle width=11.87217899999999pt height=22.648391699999998pt/> to a finite interval <img src="/tex/fe477a2781d275b4481790690fccd15f.svg?invert_in_darkmode&sanitize=true" align=middle width=32.18228144999999pt height=24.65753399999998pt/>. This is totally for the purpose of numerical simulation. We choose:

b = the largest  value a path has ever visited

a = the smallest value a path has ever visited


Now, we have reduced our problem to what we have discussed at the very beginning of this section. So then we estimate the the Markov matrix and compute its normalised eigenvector with  eigenvalue 1.



##  Birkhoff's Theorem

The Birkhoff's Theorem states that given <img src="/tex/b5723c8b7457075f9240a1bbd057a494.svg?invert_in_darkmode&sanitize=true" align=middle width=62.32495334999999pt height=24.65753399999998pt/> is a probability space and <img src="/tex/b46ab6a9a6eb1f2e11704d7cee79338e.svg?invert_in_darkmode&sanitize=true" align=middle width=78.90379034999998pt height=22.831056599999986pt/> a measurable function such that <img src="/tex/07617f9d8fe48b4a7b3f523d6730eef0.svg?invert_in_darkmode&sanitize=true" align=middle width=9.90492359999999pt height=14.15524440000002pt/> is invariant with respect to f then we have:

<img src="/tex/1a9f3fad8f92226ac0fb3fcff047730f.svg?invert_in_darkmode&sanitize=true" align=middle width=264.83056545pt height=28.894955100000008pt/>

This means that if <img src="/tex/fed0272dc708799a525db6cf5d8e9450.svg?invert_in_darkmode&sanitize=true" align=middle width=148.98183629999997pt height=24.65753399999998pt/> is finite, then

<img src="/tex/ad07fed076d905f183d732d1a6cce498.svg?invert_in_darkmode&sanitize=true" align=middle width=304.30950659999996pt height=30.648287999999997pt/>

We use the same mapping to reduce our infinite dimensional object to one dimension with finite states as what we have done in Ulam's method.

# Another Two Equations

Apart from the Logistic equation, the same scripts were also written for other two equations, which are the Predator-Pray equation, and Wreight's conjecture.

Predator-Pray:

<img src="/tex/1e1fb597f72d84098daaff10e936cf68.svg?invert_in_darkmode&sanitize=true" align=middle width=268.83810525pt height=24.65753399999998pt/>

<img src="/tex/361ec937520e1899ad9df6225b2910c3.svg?invert_in_darkmode&sanitize=true" align=middle width=278.7964179pt height=24.65753399999998pt/>


Wreight's conjecture:

<img src="/tex/d32ff6994de892241f4f4bd0eb8bb07f.svg?invert_in_darkmode&sanitize=true" align=middle width=232.42055595000002pt height=24.65753399999998pt/>









 
<!--stackedit_data:
eyJoaXN0b3J5IjpbMTAyOTY1NDcxMCwtMjIwMzIxNTUzLDEzOT
kzMzMwMDgsLTg0MDU1Njk1OSwxNjQ2NDQwMTQzLDE0MTM4NDYx
MCwtMTUzNzYyNzI5LDg0ODk1NDA0OCwxMTczNDUwNzYwLDE1Nz
U4MDMyNzIsODYyNTI1MTE4LC0yNzQ5NzgwNjYsMTkxMTYzMDk1
OCwtMTAxMzgzNzk1MCwtNjA4ODM1MzQyLC04Njc5NTE2NSwxMz
QyNjcxODY0LDI2NTg3NDE0MCwxNDQ2MjAzNDUxLC02MjE3MDIw
MzVdfQ==
-->