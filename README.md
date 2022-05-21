# Spin-dynamics via dynamic message-passing </br> (dynamic cavity-method)
This repo contains the codes used to compute out-of-equilibrium spin probability distribution for the Ising model at any given temperature. 
The codes implement a novel algorithm introduced in this paper: </br>
- <a href="https://journals.aps.org/pre/abstract/10.1103/PhysRevE.92.010102"> Dynamic message-passing approach for kinetic spin models with reversible dynamics</a> </br>
  Gino Del Ferraro and Erik Aurell, Phys. Rev. E 92, 010102(R) 
  
The code computes the single-sping probability distribution $P_i(\sigma_i(t))$ at any given time and for any given temperature. This distribution is used to compute the out-of-equilibrium magnetization $m_i(t)$. Results are compared with Monte Carlo out-of-equilibrium simulations. 

The code can, in principle, compute the single-spin off-equilibrium $P_i(\sigma_i(t))$ on any network architectures. Message-passing algorithms, though, are known to converge only on tree-like network structures. Therefore, the algorithms employed here is guaranteed to converge only on such structures. 
