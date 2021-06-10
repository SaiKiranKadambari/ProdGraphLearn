In product graph learning or clustering, a practical question that arises is, how to choose P, Q, KP, KQ? 

For many applications (a few of them are discussed next), the values of P, Q, KP , and KQ might have some physical meaning and are naturally available. 

However, in general, there is no standard recipe to choose P and Q as for a given product graph with N nodes there are many ways of choosing P and Q so that PQ = N. Nonetheless, for the clustering, we can choose KP and KQ by increasing them till the objective function value increases as illustrated in Figs. 4(d)-(g). 

In these figures, the objective function value increase after KP = 3 and KQ = 4, which correspond to the actual number of components of the
graph factors.

Use the code Figure_OptKval.m to generate the figures 4 d,e,g,f as in paper, and use the codes OptimalVal_k1 and OptimalVal_K2 to generate the corresponding data files. 