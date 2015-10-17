#
# Compute polynomial associated to Rock4
#
R=RDF
#var("x")
load("coeffs-rock4.sage")
ac=Matrix(R,5,4)
bc=vector(R,5)
def recf(i):
    return RECF[i-1]
def fpa(i,j):
    return FPA[(j-1)*50+i-1]
def fpb(i,j):
    return FPB[(j-1)*50+i-1]
def fpbe(i,j):
    return FPBE[(j-1)*50+i-1]
def ms(i):
    return MS[i-1]
def R1(p2,p1,a,b):
    return a*x*p2+b*p1
def R2(p2,p1,a,b,c):
    return a*x*p2+b*p2+c*p1
mdeg=3
mp=2
mz=mp
mr=2
P=R["x"]
pols=[1]
pol=1
pols.append(R1(pol,pol,recf(mr),1))
for i in [2..mdeg]:
    a1=recf(mr+2*(i-2)+1)
    a3=recf(mr+2*(i-2)+2)
    a2=1.-a3
    pols.append(R2(pols[i-1],pols[i-2],a1,a2,a3))
#
#
#
ac[2,1]=fpa(mz,1)
ac[3,1]=fpa(mz,2)
ac[3,2]=fpa(mz,3)
ac[4,1]=fpa(mz,4)
ac[4,2]=fpa(mz,5)
ac[4,3]=fpa(mz,6)
bc[1]=fpb(mz,1)
bc[2]=fpb(mz,2)
bc[3]=fpb(mz,3)
bc[4]=fpb(mz,4)
#
#
#
polrec=pols[mdeg]
k1=x*polrec
k2=x*(polrec+ac[2,1]*k1)
k3=x*(polrec+ac[3,1]*k1+ac[3,2]*k2)
k4=x*(polrec+ac[4,1]*k1+ac[4,2]*k2+ac[4,3]*k3)
polf=(1+bc[1]*k1+bc[2]*k2+bc[3]*k3+bc[4]*k4).factor(x)
#
racs1=solve([polf==1],[x])
racs=[i.rhs().n() for i in racs1]
