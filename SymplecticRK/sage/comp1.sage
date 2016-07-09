
#-----------------------------------------------------------------------------
# Compute the coefficients of the Gauss symplectic formula.
#-----------------------------------------------------------------------------
# We fist compute the roots of the polynomial d/dx(x^s*(x-1)^s) in QQbar
# (field of algebraic numbers). So roots are computed exactly.
# We convert these roots in the field Rpivot, which is supposed to be
# an extremly precise field of real floatting points numbers.
# At the end the fnal result is converted in the field Rout which can be RDF
# (corresponding to C language "double" or RealField(112) which corresponds
# to C "long double"
# Then we generate the specialized templates used by the C++ code.
#
# If we use directly RDF or RealField(112) a Rpivot, the precision of the
# Cooper test is restrained, which possibly indicates a loss of precision.
#

#-----------------------------------------------------------------------------

R.<x> = QQbar[] #the polynomial ring.

def Cooper(s,A,B):
#
#  Verify; theorem of Cooper (see Hairer et al.).
#  All the   [B[i][0]*A[i,j]+B[j][0]*A[j,i]-B[i][0]*B[j][0]must be (about)
#  equal zero.
    return sum([abs(B[i][0]*A[i,j]+B[j][0]*A[j,i]-B[i][0]*B[j][0]) \
                for i in range(s) for j in range(s)])
 
def lcoeffs(s,Rout,Rpivot):
    p=x^s*(x-1)^s
    l=p.derivative(x,s) #the polynomial.
    theRoots=[Rpivot((l.roots()[i][0]).real()) for i in range(0,s)]
    A=Matrix(Rpivot,s,s)
    B=Matrix(Rpivot,s,1)
    # generate a_ij and b_i coefficients:
    for j in range(0,s):
        lj=1
        for k in range(0,s):
            if k!=j:
                lj*=(x-theRoots[k])/(theRoots[j]-theRoots[k])
        Lj=lj.integral()
        for i in range(0,s):
            A[i,j]=Lj(theRoots[i])-Lj(0)
        B[j]=Lj(1)-Lj(0)
    print "Cooper,Rpivot = ",Cooper(s,A,B)
    AA=A.change_ring(Rout)
    BB=B.change_ring(Rout)
    #print AA.parent(),BB.parent()
    #
    print "Cooper,Rout = ",Cooper(s,AA,BB)
    la="initializer_list<double> la={"
    for i in range(0,s):
        for j in range(0,s):
            la+=str(Rout(A[i,j]))+","
        la+="\n"
    la=la[:-2]+"};\n"


    lb="initializer_list<double> lb={"
    for i in range(0,s):
        lb+=str(B[i][0])+","
    lb=lb[:-1]+"};\n"
    d={}
    d['la']=la
    d['lb']=lb
    return d
######################################################################
# "main" program starts here.
######################################################################

rngs=[(RDF,"double"),(RealField(112),"long double")]
      
cc="copy(la.begin(),la.end(),a);\ncopy(lb.begin(),lb.end(),b);\n}\n"

fout=open("generated_coeffs.hpp","w")
Rpivot=RealField(1024)
for rng in rngs:
    Rout=rng[0]
    nm=rng[1]
    for i in range(2,17):
        a="template<> void icoeffs<"+str(i)+","+nm+">("+nm+"* a," \
            +nm+"* b)\n{\n"
        print "\n",nm,i
        ll=lcoeffs(i,Rout,Rpivot)
        a+=ll["la"]+"\n"+ll["lb"]+"\n"
        #
        fout.write(a)
        fout.write(cc)
        fout.write("\n")
fout.close()
    
