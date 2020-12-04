This folder has the files and data required for Learning graph factors with rank constraints. 

Use the following codes for generating the results. 


1. RPGL_method.m is the proposed RPGL method. We save the results from this code as ResultsRPGL.mat

2. PGL_method.m is for obtained the factor graphs without rank constraints.

2. BigLASSO.m is for obtaining the graph learning performance using the method
 A. Kalaitzis, J. Lafferty, N. D. Lawrence, and S. Zhou, “The bigraphicallasso,” in Proc.  of the 30th  Int.  Conf.  on  Machine  Learning,  vol.  28, no. 3, Atlanta, Georgia,    USA, June 2013.

Since the precision matrix estimates from BiGLasso are not valid Laplacian matrices, we project them onto the set of all the valid graph Laplacian matrices (ProjectedBiGLasso).

3. DongsMethod.m is the code is for obtaining the graph learning performance using the method:
 X.Dong,  D.Thanou,  M.Rabbat,  and  P.Frossard, “Learning graphs from data: A signal representation perspective,” IEEE Signal Process.Mag., vol. 36, no. 3, pp. 44–63, May 2019.
 Here, we obtain the product graph Ln with N nodes. To get the graph factors Lp and Lq with rank constraints, we factorize Ln as Lp and Lq using the water-filling based KronSum  Factorization solver with spectral constraints (waterfill_RKF.m). 

4. NoiseFree.m is for obtaining the graph factors Lp and Lq from the noise-free Ln.

5. NavieMethod.m- In this code, we will test the clustering performance of the graph data’s navie clustering approaches.
 We evaluate the clustering performance on the k-means and spectral clustering algorithms.
 
 For the k-means clustering algorithm, we refer to
 S. Lloyd, “Least squares quantization in PCM,” IEEE transactions on information theory, vol. 28, no. 2, pp. 129–137, 1982.

 For the spectral clustering algorithm, we refer to
 A. Y. Ng, M. I. Jordan, and Y. Weiss, “On spectral clustering: Analysis and an algorithm,” 
 in Proc. of the Adv. in neural inf. Proc. sys. (NIPS),Vancouver, Canada, Dec. 2002

6. The folder Fig3 Data has the Matlab codes and mat files for obtaining Fig 3 of the paper.

