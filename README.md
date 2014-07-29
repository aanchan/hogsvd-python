hogsvd-python
=============



This script calculates the generalized SVD(GSVD) for N>2 matrices. This is referred to as the higher-order
GSVD(HOGSVD). This is based upon work done by Ponnapalli, Alter et al in the paper: 
Ponnapalli, Sri Priya, et al. "A higher-order generalized singular value decomposition for comparison of global mRNA expression 
from multiple organisms." PloS one 6.12 (2011): e28072.

Usage : "hogsvd.py -m textFileWithMatrixD1 -m textFileWithMatrixD2 -m textFileWithMatrixD3 -o outputDir"

Notes:

1.outputDir will contain the text files, with U and Sigma matrices for each of the input matrices specified by the -m switch. 
A single file in the output directory is created for the shared subspace V.

2. When only two matrices are specified, the QR-decomposition based algorithm is called from LAPACK.

3. When three or more matrices are specified the algorithm specified in the paper is called. This algorithm in *not* based on the
   QR decomposition.
   a. The algorithm is as follows(ref. Ponnapalli et al.):

      i. For each input matrix Di form Ai=Di^T*Di

      ii. Using all Ai,Aj pairs calculate the balanced sum Sij

      iii. Sum over all i,j ~ Si,j to get S

      iv. Do eigen-decomposition on S --> get V and L(Lambda).

      v. Solve a linear system of equations for each Di and V to get Bi

      vi. Normalize columns of Bi to get Ui, norm of the columns form the elements of elements of the diagonal matrix Sigma_i

4. To generate the sggsvd binding use f2py like this:

   Step 1: Generate Signature file:

   f2py ../lapack-3.5.0/SRC/sggsvd.f -m sggsvd2 -h sig.pyf

   Step 2: Modify sig.pyf to generate sig2.pyf (Trick is in declaring the input and output variables using input=intent(in),and output=intent(in)+intent(out)

   Step 3: Use sig2.pyf to generate the ssgsvd2 module

   f2py -L../lapack-3.5.0/liblapack.a -llapack -c sig2.pyf ../lapack-3.5.0/SRC/sggsvd.f

   Step 4: import ssgsvd2 in the python script as in hogsvd.py 

```
   #############################################################
   Note: Comparing generalized-SVD N=2 between LAPACK and MATLAB
   #############################################################
   A = [ 0.2288177 ,  0.87973283,  0.66893635,  0.97611653,  0.18300222],
       [ 0.83183603,  0.9563369 ,  0.49098867,  0.62951634,  0.06390671],
       [ 0.56031347,  0.33504639,  0.94928925,  0.62597671,  0.26223665],
       [ 0.94485407,  0.30746061,  0.97271252,  0.69691154,  0.18479945],
       [ 0.52239166,  0.72090532,  0.39089873,  0.88321032,  0.34586975]

   B = [ 0.69928475,  0.97095874,  0.3716468 ,  0.1987496 ,  0.85663137],
       [ 0.89528486,  0.67339985,  0.32713991,  0.09805244,  0.39415708],
       [ 0.92936619,  0.55483178,  0.25689855,  0.29377181,  0.88427309],
       [ 0.53662783,  0.5599302 ,  0.80107268,  0.07220077,  0.5653766 ],
       [ 0.1618458 ,  0.27674447,  0.65759925,  0.04619592,  0.32451514]


   #######
   LAPACK
   #######
   U - Left basis for Matrix A

   [-0.59489739  0.26750204  0.45597154  0.45779526  0.39629957]
   [-0.39530799 -0.8765651   0.04360846  0.15031476 -0.22554204]
   [-0.34707662  0.2788839   0.39976519 -0.42321956 -0.68031949]
   [-0.38386795 -0.09590145 -0.15791212 -0.73824096  0.52298397]
   [-0.47118184  0.27037612 -0.77810037  0.2090802  -0.23607211]

   V - Left basis for Matrix B

   [-0.02534155  0.0219933  -0.05234037  0.78979003 -0.61021781]
   [-0.21906021 -0.88348252 -0.10359917 -0.23836443 -0.32236814]
   [-0.05505875  0.37917852 -0.6710071  -0.42269319 -0.47357407]
   [ 0.67014527  0.03052208  0.48630968 -0.30824611 -0.46739778]
   [-0.70657188  0.27252084  0.54752225 -0.21384235 -0.28456864]
   
   X - Common right basis for both A and B - Can be composed (I think...) from Q and the two R 
   matrices (upper triangular) corresponding R_A and R_B
   
   R_A

   [ 1.67658412  0.25851229  0.2992408   0.32308954 -2.47159886]
   [ 0.          0.58611751 -0.21575601 -0.03035932 -0.42836428]
   [ 0.          0.         -0.72558028 -0.24683948 -0.21623912]
   [ 0.          0.          0.          0.91359884 -0.48495629]
   [ 0.          0.          0.          0.         -2.70593882]
   
   R_B

   [ 0.02610406  0.00402498  0.00465912  0.00503044 -0.03848226]
   [ 0.          0.36239293 -0.13340063 -0.01877097 -0.26485506]
   [ 0.          0.         -0.72558028 -0.24683948 -0.21623912]
   [ 0.          0.          0.          0.35832787 -0.19020751]
   [ 0.          0.          0.          0.         -2.70593882]
   
   Q

   [ 0.07457899 -0.40906438  0.59870666 -0.42310244  0.53818101]
   [-0.12090214 -0.4321093  -0.27820617  0.66871464  0.52353233]
   [-0.09388376  0.03826843 -0.71264189 -0.58332396  0.37629262]
   [-0.92798996  0.28727585  0.20130633  0.00328219  0.12558617]
   [ 0.33142242  0.74963719  0.12557456  0.18310411  0.52811635]  

   Singular Values

   Alpha - Corresponding to matrix A

   [ 0.99987876,  0.78594691,  0.4280833 ,  0.91987318,  0.07368365]

   Beta - Corresponding to matrix B

   [ 0.01556979,  0.618294  ,  0.90373927,  0.39221576,  0.99728167]

   #################
   MATLAB
   #################
   U =

   -0.3963    0.4560   -0.2675   -0.4578    0.5949
    0.2255    0.0436    0.8766   -0.1503    0.3953
    0.6803    0.3998   -0.2789    0.4232    0.3471
   -0.5230   -0.1579    0.0959    0.7382    0.3839
    0.2361   -0.7781   -0.2704   -0.2091    0.4712


    V =

    0.6102   -0.0523   -0.0220   -0.7898    0.0253
    0.3224   -0.1036    0.8835    0.2384    0.2191
    0.4736   -0.6710   -0.3792    0.4227    0.0551
    0.4674    0.4863   -0.0305    0.3082   -0.6701
    0.2846    0.5475   -0.2725    0.2138    0.7066


    X =

    1.4563   -0.4463    0.5866    0.6475    1.2684
    1.4166   -0.0764    0.4378   -0.3570    1.4756
    1.0182    0.5797   -0.0327    0.7154    1.4793
    0.3398   -0.1740   -0.0710    0.0579    1.7307
    1.4291   -0.2505   -0.1805    0.0888    0.4591


    C =

    0.0737         0         0         0         0
         0    0.4281         0         0         0
         0         0    0.7859         0         0
         0         0         0    0.9199         0
         0         0         0         0    0.9999


    S =

    0.9973         0         0         0         0
         0    0.9037         0         0         0
         0         0    0.6183         0         0
         0         0         0    0.3922         0
         0         0         0         0    0.0156

    ######################
    Sanity Check: Singular Values are unique for both
    implementations, though the order is of the singular
    values is not the same. Also the values of Left and 
    Right bases need not be unique just as in the SVD.	   
    #####################

```