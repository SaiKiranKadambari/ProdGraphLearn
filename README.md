# Product Graph Learning and Clustering

We porpose an efficient graph learning algorithm, to learn the factor the product graph from multi-domain graph data.
Product graphs are useful to explain the complex realtions in the multidomain graph data.

Cartesian product
graphs are useful to explain complex relationships in multidomain
graph data. For instance, consider a rating matrix in a
recommender system such as Netflix or Amazon. It has a user
dimension as well as an item or a movie dimension. The graph
underlying such multi-domain data can often be factorized
into a user graph and an item or a movie graph.
This work aims, learing the fators of the product graph from multi-domian graph data. We assume that, the product graph can be factorized as the Cartesian product of two smaller graphs and assume that the available graph data is smooth of the product graph.
We pose the product graph learning problem as a convex optimization problem and solve efficiently using the waterfilling like solution.
We extend this solver to infer multi-component graph factors with applications to product graph clustering by introducing rank constraints on the graph Laplacian matrices. 

Although working with smaller graph factors
is computationally more attractive, not all graphs may readily
admit an exact Cartesian product factorization. To this end, we
propose efficient algorithms to approximate a graph by a nearest
Cartesian product of two smaller graphs. The efficacy of the
developed framework is demonstrated using several numerical
experiments on synthetic data and real data.
