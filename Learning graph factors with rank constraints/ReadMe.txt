This folder has the files and data required for Learning graph factors with rank constraints 

Use the following codes for generating the results 


1. RPGL_method.m is the proposed RPGL method. We save the results from this code as ResultsRPGL.mat

2. PGL_method.m is for obtained the fatcor graphs with out rank constraints.

2. BigLASSO.m is for obtaining the graph learing performance using the method
 A. Kalaitzis, J. Lafferty, N. D. Lawrence, and S. Zhou, “The bigraphicallasso,” 
 inProc.  of  the  30th  Int.  Conf.  on  Machine  Learning,  vol.  28,no. 3, 
 Atlanta, Georgia, USA, June 2013.


 Since the precision matrix estimates from BiGLasso are not valid Laplacian matrices,
 we project them onto a set of all the valid graph Laplacian matrices (Projected BiGLasso).

 We save the results as ResultsBIGLASSO.mat

3. DongsMethod.m is the code is for obtaining the graph learing performance using the method:
 X.Dong,  D.Thanou,  M.Rabbat,  and  P.Frossard, “Learning graphs from data: A signal representation 
 perspective,” IEEE Signal Process.Mag., vol. 36, no. 3, pp. 44–63, May 2019.

4. NoiseFree.m is for obtaining the graph facors Lp and Lq from the noise free Ln.
 We also report the clustering accuracy of the obtained graph factors

5. In this code, we will test the clustering performace of the navie clustering approaches given the graph data.
 We evaluate the clustering performance on the kmeans and spectral clustering algotithms.
 For k-means clustering algorithm,we refer to
 S. Lloyd, “Least squares quantization in pcm,” IEEE transactions on
 information theory, vol. 28, no. 2, pp. 129–137, 1982.

 For spectral clustering algorithm we refer to
 A. Y. Ng, M. I. Jordan, and Y. Weiss, “On spectral clustering: Analysisand an algorithm,” 
 in Proc. of the Adv. in neural inf. proc. sys. (NIPS),Vancouver, Canada, Dec. 2002



6. The folder Fig3 Data has the matlab codes and mat files for obtaining the Fig 3 of the paper.

