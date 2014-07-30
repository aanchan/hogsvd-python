import numpy as np

def writeMat(saveDir,matName,mat,mode="T"):
    matName="%s/%s"%(saveDir,matName)
    if mode=="T":
        np.savetxt(matName,mat)
    elif mode=="B":
        np.save(matName,mat)
