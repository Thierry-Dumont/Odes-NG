from readCoeffs import *
#
# fonction mdegre comme dans le code Fortran/C++
def mdegre(mdeg):
    mp2=1;
    mdegr=mdeg
    for i in range(1,51):
        if ms(i)/mdegr>=1:
            mdegr=ms(i);
            mp1=i;
            break;
        else:
	    mp2=mp2+ms(i)*2-1;
    return mdegr,mp1,mp2
#
# indexation des tableaux de coefficients.
#
def imat(i,j):
    return i-1,j-1
#
# calcul du polynome pour mdeg0 >=5.
#
def pol(mdeg0):
    mdeg,mz,mr=mdegre(mdeg0-4)
    print(mdeg,mz,mr)
    R.<z>=PolynomialRing(RDF)
    Polid=R(1)
    g=[Polid]
    g.append(1+recf(mr)*z)
    for i in range(2,mdeg+1):
        a=recf(mr+2*(i-2)+1)
        c=-recf(mr+2*(i-2)+2)
        b=1-c
        g.append((a*z+b)*g[i-1]+c*g[i-2])
    polrec=g[mdeg]

    Ac=matrix(RDF,4,4)
    Ac[imat(2,1)]=fpa(mz,1)
    Ac[imat(3,1)]=fpa(mz,2)
    Ac[imat(3,2)]=fpa(mz,3)
    Ac[imat(4,1)]=fpa(mz,4)
    Ac[imat(4,2)]=fpa(mz,5)
    Ac[imat(4,3)]=fpa(mz,6)
    bc=vector(RDF,[fpb(mz,1),fpb(mz,2),fpb(mz,3),fpb(mz,4)])

    k=[]
    k.append(Polid)
    for i in range(1,4):
        q=Polid
        for j in range(0,i):
            q+=z*Ac[i,j]*k[j]
        k.append(q)

    p=Polid+sum(bc[i]*z*k[i] for i in range(0,4))

    pol=p*polrec
    return pol
#
# main
#
mdeg0max=Integer(input("mdeg0 max:"))
pols={}
for i in range(5,mdeg0max+1):
    pols[i]= pol(i)
#
#find interesting interval for plotting.
#
pmax=pols[mdeg0max]-1
r=pmax.roots()
r.reverse()
zmin=r[1][0]

p=sum(plot(pols[i],zmin,0) for i in range(5,mdeg0max+1))

R.<z>=PolynomialRing(RDF)
q=sum(plot(pols[i]-exp(z),zmin,0) for i in range(5,mdeg0max+1))
