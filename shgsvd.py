import sggsvd2 as sggsvd
import numpy as np


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
    return(jobu,jobv,jobq,m,n,p,k,l,lda,ldb,alpha,beta,ldu,U,ldv,V,ldq,Q,work,iwork,info)

def gsvd(A,B):
    #Generate command line options for calling LAPACK Fortran Routine for two matrices
    (jobu,jobv,jobq,m,n,p,k,l,lda,ldb,alpha,beta,ldu,U,ldv,V,ldq,Q,work,iwork,info)=sggsvdOpts(A,B)
    #Call the Fortran wrapped generalized SVD routine for real-numbered matrices
    k,l,a,b,alpha,beta,u,v,q,work,iwork,info = sggsvd.sggsvd(jobu=jobu,jobv=jobv,jobq=jobq,m=m,n=n,p=p,
                                                             k=k,l=l,a=A,b=B,alpha=alpha,beta=beta,u=U,
                                                             v=V,q=Q,work=work,iwork=iwork,info=info,lda=lda,
                                                             ldb=ldb,ldu=ldu,ldv=ldv,ldq=ldq)
    #a:contains a part the upper-triangular matrix R(look at doc of the lapack function)
    #b:contains a part the upper-triangular matrix R(look at doc of the lapack function)
    #If col-rank A == col-rank B then a and b contain the same matrix R 
    #u contains the left singular vectors of A
    #v contains the left singular vectors of B
    #q contains matrix Q
    #alpha, beta -> singular values of matrix A and B respectively.
    
    #get the ratio of the singular values
    
    return(a,b,alpha,beta,u,v,q)
    
