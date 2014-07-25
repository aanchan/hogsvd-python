#!/usr/bin/python

"""
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
4. sggsvd(jobu,jobv,jobq,m,n,p,k,l,a,b,alpha,beta,u,v,q,work,iwork,info,lda=shape(a,0),ldb=shape(b,0),ldu=shape(u,0),ldv=shape(v,0),ldq=shape(q,0))

"""

import numpy as np
import sys
#from cmdLineTool import processCommandLine
import sggsvd2 as sggsvd
usage = "hogsvd.py -m textFileWithMatrixA -m textFileWithMatrixB -o outputDir"


def sggsvdOpts(A,B):
    #set options as per the LAPACK documentation of sggsvd
    jobu="U"
    jobv="V"
    jobq="Q"
    m=A.shape[0] #number of rows of the matrix A
    n=A.shape[1] #number of columns of the matrix A
    p=B.shape[0] #number of rows of the matrix B
    k=0 # (output) Integer
    l=0 # (output) Integer, k,l specify the the dimension of the sub-blocks. k+l = effective numerical rank of (A^H,B^H)^H
    lda=max(m,1) # leading dimension of array A
    ldb=max(p,1) # leading dimension of array B
    alpha=np.zeros(n) #(output)
    beta=np.zeros(n) #(output) alpha and beta contain generalized singular values of the matrix A and B
    if jobu=="U":
        ldu=max(m,1)
        U=np.zeros(ldu*m) #contains m-by-m orthogonal matrix U
        U=U.reshape(ldu,m) 
    if jobv=="V":
        ldv=max(1,p)
        V=np.zeros(ldv*p) #contains p-by-p orthogonal matrix V
        V=V.reshape(ldv,p)
    if jobq=="Q":
        ldq=max(1,n)
        Q=np.zeros(n*n) # contains n-by-n orthogonal matrix Q
        Q=Q.reshape(n,n)
    work=np.zeros(max(3*n,m,p)+n) #
    iwork=np.zeros(n,dtype=np.int)
    info=0
    #return(jobu,jobv,jobq,m,n,p,lda,ldb,ldu,ldv,ldq)
    return(jobu,jobv,jobq,m,n,p,k,l,lda,ldb,alpha,beta,ldu,U,ldv,V,ldq,Q,work,iwork,info)
    
def main(argv):
    A=np.random.rand(25)
    A=A.reshape(5,5)
    B=np.random.rand(25)
    B=B.reshape(5,5)
    #(jobu,jobv,jobq,m,n,p,lda,ldb,ldu,ldv,ldq)=sggsvdOpts(A,B)
    (jobu,jobv,jobq,m,n,p,k,l,lda,ldb,alpha,beta,ldu,U,ldv,V,ldq,Q,work,iwork,info)=sggsvdOpts(A,B)
    #k,l,a,b,alpha,beta,u,v,q,work,iwork,info = sggsvd.sggsvd(jobu,jobv,jobq,m,n,p,A,B,lda,ldb,ldu,ldv,ldq)
    #sggsvd.sggsvd(jobu=jobu,jobv=jobv,jobq=jobq,m=m,n=n,p=p,k=k,l=l,a=A,b=B,alpha=alpha,beta=beta,u=U,v=V,q=Q,work=work,iwork=iwork,info=info,lda=lda,ldb=ldb,ldu=ldu,ldv=ldv,ldq=ldq)
    k,l,a,b,alpha,beta,u,v,q,work,iwork,info = sggsvd.sggsvd(jobu=jobu,jobv=jobv,jobq=jobq,m=m,n=n,p=p,k=k,l=l,a=A,b=B,alpha=alpha,beta=beta,u=U,v=V,q=Q,work=work,iwork=iwork,info=info,lda=lda,ldb=ldb,ldu=ldu,ldv=ldv,ldq=ldq)
    print usage






if __name__=="__main__":
    main(sys.argv)
