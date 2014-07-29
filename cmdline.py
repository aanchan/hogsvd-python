import sys,getopt
import numpy as np

def processCommandLine(argv):
    
    scriptName=argv[0]
    
    outDir=None
    matList=None
    useLapack=False
    usage = "%s -m textFileWithMatrixA -m textFileWithMatrixB -m textFileWithMatrixD3 -o outputDir -l lapack"%(scriptName)
    try:
        opts, args = getopt.getopt(argv[1:], "ho:l:m:",["outputDir=","mat=","lapack="])
    except getopt.GetoptError:
        print usage
        sys.exit(2)
    if len(argv[1:]) % 2 != 0:
        print "Incorrect Usage"
        print "Usage:%s"%(usage)
        sys.exit()
    if not opts:
        print "No command line arguments specified"
        print usage
        sys.exit()
    for opt, arg in opts:
        if opt=='-h':
           print usage
           sys.exit()
        elif opt in ("-m","--mat"):
            if arg:
                if matList==None:
                    matList=list()
                try:
                    matTemp=np.loadtxt(arg)
                    matList.append(matTemp)
                except:
                    print "Unable to open file specified by:%s"%arg
            else:
                print "Matrix file name not specified\n"
                sys.exit()
        elif opt in ("-l","--lapack"):
            useLapack=True
        elif opt in ("-o","--outputDir"):
            if arg:
                outDir=arg
            else:
                print "output-directory not specified\n"
                sys.exit()
    return matList,outDir,useLapack
        
