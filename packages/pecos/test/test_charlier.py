#!/usr/bin/python
import math

a = 0.6

def fval(x, order):
  global a

  if order == 0:
    return 1.0

  elif order == 1:
    return (a-x)/a

  elif order == 2:
    return (a*a-x-2.0*a*x+x*x)/a/a

  elif order == 3:
    fm2 = fval(x, order-2)
    fm1 = fval(x, order-1)
    fm = ((order-1+a-x)*fm1 - (order-1)*fm2)/a
    #print 'x = ',x,'\tNumerical = ', str(fm), '\tExact = ', exact3(x)
    return fm

  elif order > 3:
    fm2 = fval(x, order-2)
    fm1 = fval(x, order-1)
    return ((order-1+a-x)*fm1 - (order-1)*fm2)/a

  else:
    raise("Bad order")


def exact3(x):
  global a
  return -(pow(a,3)-2.0*x-3.0*a*x-3.0*a*a*x+3.0*x*x+3.0*a*x*x-pow(x,3))/pow(a,3)


def exact(order):
  global a
  return math.exp(a)/math.pow(a,order)*math.factorial(order)


def doSum(N, o1, o2):
  global a
  sum = 0.0
  for x in range(0,N):
    sum = sum + a**x/math.factorial(x)*fval(x,o1)*fval(x,o2)
  return sum



def main():
  order = 8
  for i in range(0,30):
    sum = doSum(i,order,order)
    print i, sum, math.fabs(sum - exact(order))
  print 'Exact = '+str(exact(order))


main()
