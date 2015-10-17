def doracf(n,A,b):
    z=var("z")
    II=vector([1 for i in range(0,n)])
    AA=identity_matrix(n)-z*A
    AAi=AA.inverse()
    v=AAi*b
    q1=(z*b)*v
    return 1+q1,z
# Crouzeix
n=3
gamma=var("gamma")
l1=[gamma,0,0]
l2=[1/2-gamma,gamma,0]
l3=[2*gamma,1-4*gamma,gamma]
A=matrix([l1,l2,l3])
delta=1/6/(2*gamma-1)^2
b=vector([delta,1-2*delta,delta])
R=gamma.parent()
RacfC,z=doracf(n,A,b)
g=(cos(pi/18)/sqrt(3)+1/2)
R1=RacfC.substitute(gamma=g.n())
R2=R1.collect(z)
RF=FractionField(PolynomialRing(RDF,"z"))
RCc=RF(R2)