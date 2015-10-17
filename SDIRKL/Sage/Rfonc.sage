from sage.symbolic.expression_conversions import AlgebraicConverter
def doracf(n,A,b,F):
    R.<z>=PolynomialRing(F)
    II=vector(F,[1 for i in range(0,n)])
    AA=identity_matrix(F,n)-z*A
    AAi=AA.inverse()
    v=AAi*b
    q1=(z*b)*v
    return 1-q1,z
print("SDIRK,ordre 4")
n=5
l1=[1/4,0,0,0,0]
l2=[1/2,1/4,0,0,0]
l3=[17/50,-1/25,1/4,0,0]
l4=[371/1360,-137/2720,15/544,1/4,0]
l5=[25/24,-49/48,125/16,-85/12,1/4]
A=matrix(QQ,[l1,l2,l3,l4,l5])
b=vector(l5)
Racf4,z=doracf(n,A,b,QQ)
#Ropp et Shadid
n=2
gamma=1+1/sqrt(2)
l1=[gamma,0]
l2=[1-2*gamma,gamma]
A=matrix(RDF,[l1,l2])
b=vector(RDF,[1/2,1/2])
RacfRS,z=doracf(n,A,b,RDF)

