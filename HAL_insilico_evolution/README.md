# HAL *in silico* evolution

This script performs *in silico* evolution of hallucinated protein homooligomers (HALs) using AF2 confidence metrics (=fitness) as selection pressure. 

To maximize the exploration of the design space while minimizing use of computational resources, we used the following evolution-based computational strategy: many short MCMC trajectories (< 50 steps) outputs are clustered by fitness score, and then used to seed new trajectories (=generations).

This approach was used to generate HALs of greater complexities across longer length-scales by extending the design specifications to structures of higher symmetry (up to C42) and longer oligomeric assembly sequence lengths (up to 1800 residues). To generate multiple possible oligomers from a single structure, we specified the MCMC trajectories as single-chains with internal sequence symmetry; the resulting structure-symmetric repeat proteins can be split into any desired oligomeric assembly compatible with factorization (e.g. C15 into a pentamer, shorthanded as C15-5). 

Using this approach, we hallucinated cyclic homo-oligomers from C5 to C42 with their largest dimension ranging from 7 to 14 nm (median: 10 nm), which were then divided into homo-trimers, tetramers, pentamers, hexamers, heptamers, octamers, and dodecamer, and the backbones were re-designed with [ProteinMPNN](https://doi.org/10.1126/science.add2187).


