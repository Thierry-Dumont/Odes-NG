def outa(n,A):
    for i in range(0,n):
        for j in range(0,i+1):
            s="A["+str(i)+"]["+str(j)+"]="+str(A[i,j])+";"
            print(s)
def outb(n,b):
    for i in range(0,n):
        s="b["+str(i)+"]="+str(b[i])+";"
        print(s)
def sumline(n,A):
    for i in range(0,n):
        su=sum(A[i,j] for j in range(0,n))
        s="S["+str(i)+"]="+str(su.n())+";"
        print(s)
print("SDIRK,ordre 4")
n=5
l1=[1/4,0,0,0,0]
l2=[1/2,1/4,0,0,0]
l3=[17/50,-1/25,1/4,0,0]
l4=[371/1360,-137/2720,15/544,1/4,0]
l5=[25/24,-49/48,125/16,-85/12,1/4]
A=matrix([l1,l2,l3,l4,l5])
print(A)
outa(n,A.n())
sumline(n,A)
Ai=A.inverse()
#print(Ai)
b=vector(l5)
#print(b)
d=b*Ai;
print(b*Ai)
outb(5,d)
#-----
print("Crouzeix:")
n=3
gamma=cos(pi/18)/sqrt(3)+1/2
print("gamma= ",gamma)
l1=[gamma,0,0]
l2=[1/2-gamma,gamma,0]
l3=[2*gamma,1-4*gamma,gamma]
A=matrix([l1,l2,l3])
print(A)
outa(3,A.n())
sumline(n,A)

Ai=A.inverse()
delta=1/6/(2*gamma-1)^2
b=vector([delta,1-2*delta,delta])
print(b)
d=b*Ai
print(d.n())
outb(3,d.n())
#--------
print("Ropp et Shadid")
n=2
gamma=1+1/sqrt(2)
print("gamma= ",gamma)
l1=[gamma,0]
l2=[1-2*gamma,gamma]
A=matrix([l1,l2])
outa(n,A.n())
sumline(n,A)

Ai=A.inverse()
b=vector([1/2,1/2])
print(b)
d=b*Ai
print(d.n())
outb(n,d.n())
