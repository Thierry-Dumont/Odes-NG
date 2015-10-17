nrecf=4382
nfpa=50*6
nfpb=50*4
nfpbe=50*5
MS=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,	22,24,26,28,30,32,34,36,38,41,44,47,50,53,56,59,63,67,71,76,81,86,92,98,105,112,120,129,138,148]
#
allitems=[]
file = open('coeffs', 'r')
for line in file:
    items=line.replace("\n","").split(",")
    for x in items:
        if x!="":
            allitems.append(float(x))
RECF=[allitems[i] for i in range(0,nrecf)]
nn=nrecf+nfpa
FPA=[allitems[i] for i in range(nrecf,nn)]
n1=nn+nfpb
FPB=[allitems[i] for i in range(nn,n1)]
n2=n1+nfpbe
FPBE=[allitems[i] for i in range(n1,n2)]
#
# access to arrays:
#
def recf(i):
    return RECF[i-1];
def fpa(i,j):
    return FPA[(j-1)*50+i-1];
def fpb(i,j):
    return FPB[(j-1)*50+i-1];
def fpbe(i,j):
    return FPBE[(j-1)*50+i-1];
def ms(i):
    return MS[i-1];
