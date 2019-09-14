# Spatial Model Solver

This is code that solves for the equilibrium of "my" spatial trade model. It bases upon Allen, Arkolakis (2014). I start with a version of my model without agglomeration/congestion effects on amenity utility. This is a "economy on a line" version of the model.

## Solving algorithm

For simplicity, I start with a version of the model in which <img src="/tex/cf2889bdb114fa1a6a34972c325d36bd.svg?invert_in_darkmode&sanitize=true" align=middle width=40.302373649999986pt height=22.831056599999986pt/> and <img src="/tex/3b3ab4c7f8eddc6096361d7d5269bc6d.svg?invert_in_darkmode&sanitize=true" align=middle width=45.01990349999999pt height=21.18721440000001pt/> <img src="/tex/60cd0f191d74161a0afffe3bb569399d.svg?invert_in_darkmode&sanitize=true" align=middle width=14.79567374999999pt height=22.831056599999986pt/>. Then, the equilibrium conditions would imply a system of <img src="/tex/6251b0c7ae3c532d6b603415e5a0ce66.svg?invert_in_darkmode&sanitize=true" align=middle width=43.310391299999985pt height=22.465723500000017pt/> equations as follows
<p align="center"><img src="/tex/5057afc656beeee273d7a6dec5e89a5b.svg?invert_in_darkmode&sanitize=true" align=middle width=473.86532654999996pt height=168.81070305pt/></p>
Importantly, (2) and (3) all depend on wages, so that they can be inserted into (1) so to obtain a system of N equations, which can be solved for N wages. To do so as efficient as possible, I build all functions in matrix format so that I can represent the system with a vector of equations that must equal to zero in equilibrium.
