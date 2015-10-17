var("x y")
var("A a21 a22 a31 a32 a41 a42 a43 b1 b2 b3 b4")
k1=A*x
k2=A*(x+a21*k1)
k3=A*(x+a31*k1+a32*k2)
k4=A*(x+a41*k1+a42*k2+a43*k3)
y=x+b1*k1+b2*k2+b3*k3+b4*k4
y1=y.factor(x)
p=y1.factor_list()[0][0]
q=p.collect(A)
s=q.coefficients()
print(s)
t=[s[i][0] for i in [0..4]]
h1=t[0]+A*(t[1]+A*(t[2]+A*(t[3]+A*t[4])))
print(h1)
var("c1 c2 c3 c4")
r=q.substitute({b1:c1,b2:c2,b3:c3,b4:c4})
