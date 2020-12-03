# Product Graph Learning and Clustering

This work amis, learing the fators of the product graph from multi-domian graph data. We assume that, the product graph can be factorized as the Cartesian product of two smaller graphs and assume that the available graph data is smooth of the product graph.
We pose the product graph learning problem as a convex optimization problem and solve efficiently using the waterfilling like solution.
We extend this solver to infer multi-component graph factors with applications to product graph clustering by introducing rank constraints on the graph Laplacian matrices. 

Although working with smaller graph factors
is computationally more attractive, not all graphs may readily
admit an exact Cartesian product factorization. To this end, we
propose efficient algorithms to approximate a graph by a nearest
Cartesian product of two smaller graphs. The efficacy of the
developed framework is demonstrated using several numerical
experiments on synthetic data and real data.
