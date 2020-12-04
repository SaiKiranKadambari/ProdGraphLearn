This folder has the files and data required for Learning sparse graph factors.

Use the following codes for generating the results. 

1. PGL_graphs.m is the file for the proposed product graph learning method.

2. DongsMethod.m is the file for the Dongs method:
X.Dong,  D.Thanou,  M.Rabbat,  and  P.Frossard, “Learning graphs from data: A signal representation perspective,” IEEE Signal Process.Mag., vol. 36, no. 3, pp. 44–63, May 2019.

Where we obtain the product graph Ln. We then factorize Ln as Lp and Lq using the unconstrained water-filling based KronSum Factorization solver (waterfill_KF.m). 

3. BIGLasso.m is the Matlab file for the BIGLasso method of obtaining the product graphs.
We refer to,
A. Kalaitzis, J. Lafferty, N. D. Lawrence, and S. Zhou, “The bigraphicallasso,” 
 inProc.  of the  30th  Int.  Conf.  on  Machine  Learning,  vol.  28, no. 3, 
 Atlanta, Georgia, USA, June 2013.

Since the precision matrix estimates from BiGLasso are not valid Laplacian matrices,
we project them onto a set of all the valid graph Laplacian matrices (Projected BiGLasso).
We use the BIGLasso papers code for the implementation of inverse covariance matrices. 
We then project them on to valid Laplacian matrice(LapFactorize.m is the pseudo-code).


%%%%%%%%%%%%% Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

In the folder Fig2 Data, we have the Matlab files required for generating the Fig. 2 of the paper. 
