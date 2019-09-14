# Spatial Model Solver

This is code that solves for the equilibrium of "my" spatial trade model. It bases upon Allen, Arkolakis (2014). I start with a version of my model without agglomeration/congestion effects on amenity utility. This is a "economy on a line" version of the model.

## Solving algorithm

For simplicity, I start with a version of the model in which $\beta=0$ and $\bar{u}_i=1$ $\forall i$. Then, the equilibrium conditions would imply a system of $3 \times N$ equations as follows
\begin{align}
w_i L_i 	&= \sum_{j \in S} \bigg( \frac{w_i \tau_{ij}}{A_i} \bigg)^{1-\sigma} P_j^{\sigma-1} w_j L_j	\label{eq:ap1} \\		
L_i &= \sum_{j \in S} \frac{(w_i / P_i \mu_{ij})^\theta u_i}{\sum\limits_{s \in S} (w_s / P_s \mu_{sj})^\theta u_s} L_j^0	\label{eq:ap2} \\
P_i &= \bigg( \sum_{j \in S} (w_i \tau_{ij} / A_i)^{1-\sigma} \bigg)^{\frac{1}{1-\sigma}}.	\label{eq:ap3}
\end{align}
Importantly, (2) and (3) all depend on wages, so that they can be inserted into (1) so to obtain a system of N equations, which can be solved for N wages. To do so as efficient as possible, I build all functions in matrix format so that I can represent the system with a vector of equations that must equal to zero in equilibrium. As of notation, I define circled operations as element-wise ones, i.e. circled dot (slash) for element-wise multiplication (division). Element-wise power is denoted with diamond, and in line dot stands for standard matrix multiplication.

First, I calculate a $N \times 1$ vector of price indexes, given by \cref{eq:ap3}, as
\begin{equation}
\mathbf{P}_{N\times 1} = \bigg[ \{ \mathbf{T}_{N\times N}\diamond (1-\sigma) \} \cdot \{ (\mathbf{w} \oslash \mathbf{A}) \diamond (1-\sigma)\}_{N\times 1} \bigg] \diamond (1/1-\sigma).	\label{eq:prices_vector}
\end{equation}
Note that $\mathbf{w} \oslash \mathbf{A} = [w_1/A_1, w_2/A_2, ...]'$, so that the first row of $\mathbf{T} \cdot (\mathbf{w} \oslash \mathbf{A}) =  (w_1 / A_1) \times \tau_{11} + ... + (w_N / A_N) \times \tau_{N}$. Thus, by element-wise powering both matrices to $1-\sigma$ before mutiplying them and then to $1/1-\sigma$ after doing so, one obtains $\mathbf{P} = [P_1, ..., P_N]$. The second step consists in calculating the equilibrium allocation of labor. The challenge there is to calculate $\mathbf{\Pi}$, a symmetric matrix of bilateral probabilities $\{ \Pi_{ij} \}_{i,j \in S}$. To do so, I define
\begin{align*}
\phi_{N \times 1} &= \{\mathbf{w} \oslash \mathbf{P} \diamond \theta \} \odot u_{N \times 1}, \quad \text{and} \\
\Phi_{N \times 1} &= [\mathbf{M}_{N \times N} \diamond (-\theta) ] \cdot \phi.
\end{align*}
Importantly, the first row of $\Phi = (w_1/P_1 \mu_{11})^{\theta} u_1 + (w_2/P_2 \mu_{21})^{\theta} u_2 + ... + (w_N/P_N \mu_{N1})^{\theta} u_N = \sum_s (w_s/P_s \mu_{s1})^\theta u_s$, which is exactly the denominator of $\Pi_{i1}$ $\forall i$. Moreover, $\phi$ is almost the vector of the nominators of the same function -- each of its elements need only to be multiplied by the respective $\mu_{ij}^{-\theta}$. Thus, $\mathbf{\Pi}$ and $\mathbf{L}$ can be obtained as
\begin{align}
\mathbf{\Pi}_{N \times N} &= \big[ \phi \cdot \{ \Phi \diamond (-1) \}' \big]_{N \times N} \odot \{ \mathbf{M}_{N \times N} \diamond (-\theta) \}_{N \times N}, \quad \text{and} \label{eq:transitions} \\
\mathbf{L}_{N \times 1} &= \mathbf{\Pi} \cdot \mathbf{L}^0.	\label{eq:labormkt}
\end{align}
The final step is to write the market clearing condition -- \cref{eq:ap1} -- in vector format. One can do so by noting that
\begin{align}
w_i L_i 	&= \sum_{j \in S} \bigg( \frac{w_i \tau_{ij}}{A_i} \bigg)^{1-\sigma} P_j^{\sigma-1} w_j L_j	\nonumber		 \\
		&= \bigg( \frac{w_i}{A_i} \bigg)^{1-\sigma} \sum_{j \in S} \bigg( \frac{\tau_{ij}}{P_j} \bigg)^{1-\sigma} w_j L_j \rightarrow	\nonumber		 \\
w_i^\sigma	&= \frac{1}{L_i A_i^{1-\sigma}} \sum_{j \in S} \bigg( \frac{\tau_{ij}}{P_j} \bigg)^{1-\sigma} w_j L_j. \label{eq:eq_cond}
\end{align}
By defining
\begin{align*}
\mathbf{\Psi}_{N \times 1} 	&= \{ \mathbf{T} \diamond (1 - \sigma) \}_{N \times N} \cdot \{ \mathbf{w} \odot \mathbf{L} \oslash [ \mathbf{P} \diamond (1 - \sigma) ] \}_{N \times 1}, \quad \text{and} \\
\Omega_{N \times 1}		&= \{ \mathbf{L} \odot [\mathbf{A} \diamond (1-\sigma)] \} \diamond (-1)
\end{align*}
one would note that its first row of $\Omega \cdot \mathbf{\Psi} = (L_1 A_1^{1-\sigma})^{-1} \sum_{j \in S} ( \tau_{1j}/P_j )^{1-\sigma} w_j L_j$, which is the RHS of \cref{eq:eq_cond} for $i=1$. Then, given a certain geography -- vector of productivities, trade and migration bilateral frictions and initial labor distribution -- in equilibruim the following condition should hold:
\begin{equation}
f(\mathbf{w}, \mathcal{G}(S)) = \mathbf{w}_{N \times 1} \diamond \sigma - [ \Omega \cdot \mathbf{\Psi} ]_{N \times 1} = \mathbf{0}_{N \times 1}.	\label{eq:obj_function}
\end{equation}
To find equilibrium wages, I search for a $N\times 1$ vector that minimizes the square of \cref{eq:obj_function}; thus
\begin{equation}
\mathbf{w}^* = \text{arg}\min\limits_{w \in R_{+}^N} f(w, \mathcal{G}(S))' \cdot f(w, \mathcal{G}(S)).
\end{equation}
To close the solution, the equilibrium prices and labor allocations can be found by inserting $\mathbf{w}^*$ into \cref{eq:prices_vector,eq:transitions,eq:labormkt}.
