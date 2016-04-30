


origA = [2 1 0 0 0 ; 
     1 2 1 0 0 ;
     0 1 2 1 0 ;
     0 0 1 2 1 ;
     0 0 0 1 2]


q1=[1,1,1,1,1;
    1,1,1,1,1]'
x=[0,0,0,0,0]'
rhs = [1,1,1,1,1]';

Proj = eye(5) - q1*q1'

A = Proj*origA

r = rhs
p = r



k=1
res = norm(r)

pAp = p'*A*p

alpha = norm(r)/pAp

x=x+alpha*p

r=r-alpha*A*p

newrnorm = norm(r)

beta = newrnorm/res

p = r+beta*p








k=2
res = norm(r)
pAp = p'*A*p

alpha = norm(r)/pAp

x=x+alpha*p

r=r-alpha*A*p

newrnorm = norm(r)

beta = newrnorm/res

p = r+beta*p





k=3
res = norm(r)
pAp = p'*A*p

alpha = norm(r)/pAp

x=x+alpha*p

r=r-alpha*A*p

newrnorm = norm(r)

beta = newrnorm/res

p = r+beta*p




k=4
res = norm(r)
pAp = p'*A*p

alpha = norm(r)/pAp

x=x+alpha*p

r=r-alpha*A*p

newrnorm = norm(r)

beta = newrnorm/res

p = r+beta*p




k=5
res = norm(r)
pAp = p'*A*p

alpha = norm(r)/pAp

x=x+alpha*p

r=r-alpha*A*p

newrnorm = norm(r)

beta = newrnorm/res

p = r+beta*p




