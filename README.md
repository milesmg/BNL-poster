# Accelerated Neural System Identification via Reduced-Order Modeling

Miles M. Gantcher$^1$, Nathan M. Urban$^1$


$^1$Department of Applied Mathematics, Brookhaven National Laboratory, Upton, NY 11973.


## Figures

- [Animations](#figures)
  - [Allen-Cahn Trajectory on a FOM and ROM](#figure=Figures/Animations/Allen-Cahn_trajectory.gif)
    - Galerkin projection with r = 20 spatial modes
    - DEIM hyperreduction, m = 20 DEIM points
  - [Reaction-Diffusion Trajectory on a FOM and ROM](#figure=Figures/Animations/Reaction-Diffusion_trajectory.gif)
    - Galerkin projection, r = 20 spatial modes
    - DEIM hyperreduction, m = 20 DEIM points split evenly between $s_1$ and $s_2$
  - [Cahn-Hilliard Trajectory on a FOM and ROM](#figure=Figures/Animations/Cahn-Hilliard_Trajectory.gif)
    - Petrov-Galerkin projection with a $(-\Delta)^{-1}$ trial basis, r = 5 spatial modes
    - ECSW hyperreduction, m = 30 maximum non-zero ECSW elements
<!-- - $\kappa = 0.1$, $\bar{c} = 0.0$ -->
  - **[Key Result: Learned (neural nonlinearity) Allen-Cahn Trajectory on a FOM and ROM](#figure=Figures/Animations/Allen-Cahn_trajectory.gif)**
  - [1-dimensional Allen-Cahn trajectory](#figure=Figures/Animations/1d_allen_cahn.gif)
    - [FOM vs. ROM](#figure=Figures/Animations/1d_allen_cahn.gif)
    - [Evolution of reduced coefficients over time (reduced trajectory)](#figure=Figures/Animations/1d_rom_reduced_state.gif)
    - [Evolution of spatial modes over time](#figure=Figures/Animations/1d_rom_mode_contributions.gif)
- [Images](#folder=Figures/Images)
  - [2-D reaction diffusion target function plot](#figure=Figures/Images/Nonlinearities/rd_nonlinearity.png)
  - **[Key Result: 1-D Allen-Cahn target function with learned neural approximation](#figure=Figures/Images/Poster_figures/learned_function_AC.png)**
  - [Spatial modes, A-C trajactory](#folder=Figures/Images/2D_Allen-Cahn_Spatial_Modes)
  - [Spatial modes, R-D trajectory](#folder=Figures/Images/Reaction-diffusion_spatial_modes)
  - [Spatial modes, C-H trajectory](#folder=Figures/Images/Cahn-Hilliard_Spatial_Modes)

## A Note on Gradient Flows and the Allen-Cahn and Cahn-Hilliard Equations

Let $c\in [-1,1]^{N^d}$ denote a scalar relative concentration field; that is, each point in our field, $c(x_1,...,x_d)$, denotes the relative concentration $[A]-[B]$ of two solutes. Both the Allen-Cahn and Cahn-Hilliard equation are formulated as gradient flows of the Ginzburg-Landau free energy functional on this field. That functional takes the form 
$$E[c] = \int_\Omega \frac{\kappa}{2} |\nabla c|^2 + F(c)d\Omega,$$ 
where $\Omega$ is the spatial domain. Here, the gradient term penalizes rapid changes in concentration; $F$ is the free energy term, often represented by a double well potential $F(c) = \frac{(c^2-1)^2}{4}$

The Allen-Cahn equation is an $L^2-$gradient flow; under the Allen-Cahn equation, $c_t$ is such that $E[c]$ is minimized as much as possible at each time. The Cahn-Hilliard equation is similar, but $c_t$ must preserve mass. The *nonlocal* Cahn-Hilliard equation describes the dissolution of species $A$ and $B$ when they are attached to a polymer, known as a block copolymer. This equation adds a term to the free energy which represents the polymer constraints, giving 
$$E_{nonlocal}[c] = \int_\Omega \frac {\kappa}{2}|\nabla c|^2 + F(c) + \frac{\sigma}{2}(c - \bar{c})(-\Delta ^{-1})(c-\bar{c})d\Omega,$$
where $\bar{c}$ is the mean concentration (average mass) and $-\Delta ^{-1}$ is the negative inverse of the Laplacian operator, which is defined on the zero mean subspace in which $c - \bar{c}$ lives. 

## A Note on Projection-Based Model-Order Reduction

As described on the poster, projection-based model order reduction approximates a full state vector $u \in \mathbb R^{N^d}$, where $d$ is the dimension, as a linear combination of $r$ basis elements. In proper orthogonal decomposition (POD), this basis $U$ is chosen from data by taking the first $r$ left singular vectors of a matrix of data snapshots. An initial condition $u_0$ is projected into reduced space via $\tilde{u_0} = U^Tu_0$; a state is projected up from the reduced space to the full space via $U \in R^{N^d \times r}.$ 

How do we evolve our reduced state so that its dynamics respect the dynamics of the full order equation, albeit in reduced form, and so that we can extract the resultant true state at any point? We project our equation down; specifically, if our PDE is of the form
 $$u_t = L(u) \iff u_t - L(u) = 0,$$ 
 then we would *like* to require 
$$U\tilde{u}_t = L(U \tilde{u}) \iff U\tilde{u}_t - L(U \tilde{u}) = 0$$
but this equation is over-constrained; we're projecting up from a low dimensional space and making $N^d$ requirements of the result of that projection. Thus, we instead define the residual 
$$U\tilde{u}_t - L(U \tilde{u}) =: r(\tilde{u}) \in \mathbb R^{N^d}.$$ 
We can force our residual to be orthogonal to whatever $r-$dimensional space we choose. **Galerkin projection** says: force our residual to be orthogonal to the reduced space, the span of the first $r$ singular vectors of our data matrix $\equiv$ the columns of $U$. That is, enforce $U^Tr(\tilde{u}) = 0$. **Petrov-Galerkin projection** chooses a distinct basis. This is may be helpful in the context of the Cahn-Hilliard equation, when we want force our ROM to respect the mass-preserving gradient flow construction of the Cahn-Hilliard equation. The brief sketch shown on the paper describes the process of pre-computing the projected operators of our PDE. In the case of the above equation, applying Galerkin projection leads to 
$$U^TU\tilde{u}_t = U^TL(U(\tilde{u})).$$ 
If $L$ is linear then we can precompute $\tilde{L} = U^TLU$; if, as in POD, the columns of $U$ are orthonormal, then $U^TU = I$, and we get the reduced differential equation
$$\tilde{u}_t = \tilde{L}(\tilde{u}). $$

## A Note on Hyperreduction

The two methods of hyperreduction employed in this paper are the **discrete empirical interpolation method (DEIM)** and **energy-conserving sampling and weighting (ECSW)**. Both require function 'snapshots,' evaluations of the function over the course of a trajectory. In the context of system identification, this seems unappealing, but in many cases it may be possible to build these function snapshots from a concatenation of of pre-computed trajectories from various positions in parameter space. 

**DEIM** is an approximate-then-project method. Similar to POD, the function $f$ is approximated as a linear combination of basis vectors; this basis is the first $m$ left singular vectors of $F$, the function snapshot matrix. Note that we are not approximating $f(\cdot)$ via a function basis; we're approximating each function evaluation vector $f(x) \in \mathbb R^{N^d}$ as a linear combination of our evaluation basis. We then select $m$ points out of our $N^d-$length state vector, and require that our $m-$dimensional approximation is correct at those $m$ spatial (DEIM) points for all time. If our function snapshot basis is $\Phi$, then this is equivalent to enforcing 
$$Z_m f(U\tilde{u}) = Z_m \Phi \alpha (\tilde{u}),$$
where $Z_m$ is the matrix that selects our $m$ points and $\alpha(\tilde{u})$ is the vector of coefficients with which we interpolate our function basis for a given reduced state. Thus, our approximation of $f$ is given by 
$$\Phi\alpha (\tilde {u}) = \Phi(Z_m \Phi)^{-1}Z_m f(U\tilde{u}),$$
which we then project down again via our chosen projection method. For example, in a Galerkin projection with state basis matrix $U$ of the equation $u_t = L(u) + N(u)$ where $L$ is linear and $N$ is nonlinear, we get 
$$\tilde{u}_t = \tilde{L}\tilde{u} + \tilde{V}f(U_m\tilde{u}),$$
where $\tilde{L} = U^T L U$, $V = U^T\Phi(Z_m \Phi)^{-1}Z_m,$ and $U_m$ projects $\tilde{u}$ into only the $m$ DEIM points.

The choice of DEIM points is nontrival. In many cases, they are chosen via a greedy algorithm: for each basis vector (column of $\Phi$) $\phi_i$, we include the point $m$ where our approximation of $\phi_i$ via the first ${\phi_1,...,\phi_{i-1}}$ basis vectors is the worst. This method is easy to compute, but may be suboptimal in many cases. 

**ECSW** is a project-then-approximate method. It was originally developed for a finite element context; here, we treat each grid point as an element and preserve the standard notation. ECSW computes its approximation of the nonlinearity in reduced space, unlike DEIM, and because of its construction preserves the Lagrangian structure (least action princple) of certain systems to which it is applied. This makes it particularly appealing in a Cahn-Hilliard context, where we would like to preserve the gradient flow character of the FOM in our ROM. 

Let $U$ denote our Galerkin projector, as above. ECSW is effectively a quadrature method: we would like the approximate the function 
$$f(u) = \sum_{e \in E} (L^e U)^TL^e f(Uu)$$
where $L^e$ selects the $e-$th element and $E$ is the set of elements, via 
$$f(u) \approx \tilde{f}(u) =  \sum_{e \in \tilde{E} \subset E} \xi_e (L^e U)^TL^e f(Uu).$$ 
Enforcing the nonnegativity of the $\{\xi_e\}$ makes ECSW energy conserving in the right contexts. We let $\xi$ denote the vector of these 'quadrature weights.' Let $G$ denote a matrix of per-element reduced force snapshots generated from data; let $b$ denote the matrix of their sums across elements (that is, the matrix of reduced force snapshots). Then the offline 'training' of ECSW, analogous to the SVD and greedy evaluation point selection of the DEIM algorithm, involves adjusting $\xi$ to minimize 
$$||G\xi = b||$$
 while keeping all elements nonnegative, and keeping as many **zero** as possible. 
