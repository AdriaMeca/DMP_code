# Dynamic message-passing code for time-varying networks

This repository contains the Fortran code I used to simulate epidemics on time-varying networks for my bachelor's thesis, which was supervised by Dr. Matteo Palassini. Now, it has been expanded to contain the code that I am currently using for my master's thesis.

The folder *modules* has several Fortran modules that provide procedures for creating and manipulating different networks. It also contains **dynamic message-passing** algorithms for different epidemiological models (SIR and SEIR) that provide the marginal probabilities that each node on a network is in a given state at time $t$.

The **dynamic message-passing** algorithm gives exact results on tree-like networks and good approximations on real networks. Its complexity is $\mathcal{O}(cNt)$, where $N$ is the number of nodes and $c$ is the average number of neighbors each node has.

Furthermore, the **dynamic message-passing** algorithm provides the marginal probabilities of interest with only one run per network, whereas the **Monte Carlo** approach requires a large number of runs. Therefore, the **dynamic message-passing** method is more efficient.
