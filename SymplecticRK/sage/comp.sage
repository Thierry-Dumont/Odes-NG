
#-----------------------------------------------------------------------------
# define s (number of stages) and RR (precision) here.
#-----------------------------------------------------------------------------
s=16 #number of steps of the formula. Adapt to your needs...
RR=RealField(53) #classical "double" of C++.
#RR=RealField(112) #this seems the precision (112 bits) of "long double" in g++
#-----------------------------------------------------------------------------
# DO NOT CHANGE ANYTHING BELOW THIS LINE.
#-----------------------------------------------------------------------------
R.<x> = RR[] #the polynomial ring.
p=x^s*(x-1)^s 
l=p.derivative(x,s) #the polynomial.
theRoots=[l.roots()[i][0] for i in range(0,s)]
A=Matrix(RR,s,s)
B=Matrix(RR,s,1)
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
#
#  Verify; theorem of Cooper (see Hairer et al.).
#  All the following output must be (about) equal zero.
#
print "Cooper test:"
for i in range(0,s):
    print [B[i][0]*A[i,j]+B[j][0]*A[j,i]-B[i][0]*B[j][0] for j in range(0,s)]
print "\n"
#
# Compute the coefficients for the continuous extrapolation
# (see Hairer et al. 2nd edition, page 326).
#
# the Vandermonde matrix: 
V=Matrix(RR,s,s)
for j in range(0,s):
    for i in range(0,s):
        V[i,j]=theRoots[j]^i
# the RHS:
W=Matrix(RR,s,s)
for j in range(0,s):
    for i in range(0,s):
        W[i,j]=((1+theRoots[j])^(i+1))/(i+1)
# the continuous extrapolation coefficients:
Q=V\W 
#
#---------------------------------------------------------------------------
#now generate the results in a C++ compatible syntax:
#
print "Coefficients (A):\n"
k=0
for i in range(0,s):
    ch=""
    for j in range(0,s):
        ch+="a["+str(k)+"]="+str(A[i,j])+"; "
        k+=1
        if k%2==0:
            ch+="\n"
    print ch
print "Coefficients (B):\n"
ch=""
for i in range(0,s):
    ch+="b["+str(i)+"]="+str(B[i][0])+";"
    if i%2:
        ch+="\n"
print ch
print "Coefficients (C):\n"
ch=""
for i in range(0,s):
    ch+="c["+str(i)+"]="+str(theRoots[i])+";"
    if i%2:
        ch+="\n"
print ch
print "Coefficients (Q (extrapolation)):\n"
k=0
for i in range(0,s):
    ch=""
    for j in range(0,s):
        ch+="q["+str(k)+"]="+str(Q[i,j])+"; "
        k+=1
        if k%2==0:
            ch+="\n"
    print ch
