#!/usr/bin/python


import numpy as np
def calcHOGSVD(matList):
    #assert len(matList)==3
    #Step 1: calculate Di^T*Di
    aMat=list()
    sMat=list()
    N=len(matList)
    for mat in matList:
        D=mat
        aMatTemp=np.dot(D.T,D)
        aMat.append(aMatTemp)

    #calculate pairwise Sij sums
    for i in range(len(aMat)):
        for j in range(i+1,len(aMat)):
            print i,j
            sMatTemp=0.5 * ( np.dot(aMat[i],np.linalg.inv(aMat[j])) + np.dot(aMat[j],np.linalg.inv(aMat[i])) )
            sMat.append(sMatTemp)
    
    S=np.zeros(sMat[0].shape)

    #sum all Si,j matrices to get S
    for s in sMat:
        S=S+s

    #normalize S
    S=(2.0/(N*(N-1)))*S
    

    #Do the eigen-decomposition on S
    L,V=np.linalg.eig(S)

    #Calculate the Bi matrices
    bMatList=list()
    for mat in matList:
        D=mat
        BT=np.dot(np.linalg.inv(V),D.T)
        bMatList.append(BT.T)


    #from Bi matrices calculate U and Sigma matrices
    sigList=list()
    uMatList=list()
    
    for B in bMatList:
        sig=np.sqrt(np.sum(B*B, axis=0))
        U=B/sig
        SIG=np.diag(sig)
        sigList.append(SIG)
        uMatList.append(U)
        

    return(uMatList,sigList,V)
