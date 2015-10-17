n=4
l1=[0,0,0,0]
l2=[1/2,0,0,0]
l3=[0,1/2,0,0]
l4=[0,0,1,0]
A=matrix(QQ,[l1,l2,l3,l4])
b=vector(QQ,[1/6,2/6,2/6,1/6])
R.<z>=PolynomialRing(QQ)
p=R(1)
k=[p]
for i in range(1,n):
    q=R(1)
    for j in range(0,i):
        q+=z*A[i,j]*k[j]
    k.append(q)
print(k)
p=R(1)
for i in range(0,n):
    p+=b[i]*z*k[i]
print(p)
