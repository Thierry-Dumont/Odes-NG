def doracf(n,A,b):
    z=var("z")
    II=vector([1 for i in range(0,n)])
    AA=identity_matrix(n)-z*A
    AAi=AA.inverse()
    v=AAi*b
    q1=(z*b)*v
    return 1+q1,z
# R & S.
n=2
g=1+1/sqrt(2)
gamma=var("gamma")
l1=[gamma,0]
l2=[1-2*gamma,gamma]
A=matrix([l1,l2])
b=vector([1/2,1/2])
R=gamma.parent()
RacfC,z=doracf(n,A,b)
R1=RacfC.substitute(gamma=g.n())
R2=R1.collect(z)
RF=FractionField(PolynomialRing(RDF,"z"))
RC=RF(R2)
 
