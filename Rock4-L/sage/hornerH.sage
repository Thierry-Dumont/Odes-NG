# compute polynomial associated to the 3 steps recurrence formula.
var("x y h")
var("A a21 a22 a31 a32 a41 a42 a43 b1 b2 b3 b4")
k1=A*x
k2=A*(x+h*a21*k1)
k3=A*(x+h*(a31*k1+a32*k2))
k4=A*(x+h*(a41*k1+a42*k2+a43*k3))
y=x+h*(b1*k1+b2*k2+b3*k3+b4*k4)
y1=y.factor(x)
p=y1.factor_list()[0][0]
q=p.collect(A)
s=q.coefficients()
#print(s)
t=[s[i][0]/h for i in [0..4]]
q=[g.simplify_full() for g in t]
print(q)

