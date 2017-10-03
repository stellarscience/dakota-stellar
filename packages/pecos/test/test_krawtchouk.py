#!/projects/vuq/anaconda/bin/python
import math

p = 0.10
N = 20

def fval(x, order):
  global p, N

  if order == 0:
    return 1.0

  elif order == 1:
    return (p*N-x)/(p*N)

  elif order == 2:
    return (p*p*N*(1.0-N)+(1.0-2.0*p+2.0*p*N)*x-x*x)/(p*p*N*(1.0-N))

  elif order > 2:
    om1 = order-1
    fm2 = fval(x, order-2)
    fm1 = fval(x, order-1)
    return ((p*(N-om1)+om1*(1.0-p)-x)*fm1 - om1*(1.0-p)*fm2)/(p*(N-om1))

  else:
    raise("Bad order")


def nCrBad(n,r):
  f = math.factorial
  return f(n) / f(r) / f(n-r)


def nCr(n,k):
  combination = 1.
  nmk = n - k;
  if k <= nmk:
    for i in range(k):
      combination *= float(n-i)/float(k-i);
  else:
    for i in range(nmk):
      combination *= float(n-i)/float(nmk-i);

  return combination


def poch(r,n):
  if n == 0:
    return 1.0
  else:
    poch = r
    for i in range(1,n):
      poch *= r+i;
    return poch


def exact(order):
  global p, N
  return math.pow(-1,order)*math.factorial(order)/poch(-N,order)*math.pow((1.0-p)/p,order)


def doSum(n, o1, o2):
  global p, N
  sum = 0.0
  for x in range(0,n):
    sum = sum + nCr(N,x)*p**x*(1.0-p)**(N-x)*fval(x,o1)*fval(x,o2)
  return sum



def main():
  global N
  print fval(1.23, 5)
  #print poch(-2.0,1)
  #print nCr(6,3)
  order = 3
  for i in range(0,N+1):
    sum = doSum(i,order,order)
    print i, sum, (sum - exact(order))/exact(order)*100.0
  print 'Exact = '+str(exact(order))


main()
