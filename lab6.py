from math import exp as e
from scipy import optimize 
def f(x,y):
    return x*e((7+x)/2.)+pow(y,2)*e((7+x)/2.)-8.*y*e((7+x)/2.)+23.*e((7+x)/2.)


def fun(x):
   return f(*x)
x0 = (-1, 5)
fprime = lambda x0: optimize.approx_fprime(x0, fun, 0.001)
res = optimize.minimize(fun, x0, method = 'Newton-CG', jac = fprime)
print(res)