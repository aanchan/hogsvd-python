#!/usr/bin/python


import numpy as np
def calcHOGSVD(matList):
    assert len(matList)==3
    #Step 1: calculate Di^T*Di
    aMat=list()
    sMat=list()
    for mat in matList:
        D=mat
        aMatTemp=np.dot(D.T,D)
        aMat.append(aMatTemp)

    #calculate pairwise Sij sums
    for i in range(len(aMat)):
        for j in range(i+1,len(aMat)):
            print i,j
            #sMatTemp=0.5 * ( np.dot(aMat[i],np.linalg.inv(aMat[j])) + np.dot(aMat[j],np.linalg.inv(aMat[i])) )
            #sMat.append(sMatTemp)


