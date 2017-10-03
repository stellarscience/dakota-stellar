#!/projects/vuq/anaconda/bin/python
import math

a = 4
b = 6
N = 10

def fval(x, order):
  global a, b, N

  if order == 0:
    return 1.0

  elif order == 1:
    return 1.0 + (2.0+a+b)/((-N)*(a+1.0))*x

  elif order == 2:
    return 1.0 - 2.0*(3.0+a+b)*x/(N*(a+1.0)) + (3.0+a+b)*(4.0+a+b)/((a+1.0)*(a+2.0)*N*(N-1.0))*x*(x-1.0)

  elif order > 2:
    om1 = order-1
    fm2 = fval(x, order-2)
    fm1 = fval(x, order-1)
    A = (om1+a+b+1.0)*(om1+a+1.0)*(N-om1)/((2.0*om1+a+b+1.0)*(2.0*om1+a+b+2.0))
    C = om1*(om1+a+b+N+1.0)*(om1+b)/((2.0*om1+a+b)*(2.0*om1+a+b+1.0))
    return ((A+C-x)*fm1 - C*fm2)/A

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
  global a, b, N
  return math.pow(-1,order)*poch((order+a+b+1.0),(N+1))*poch((b+1.0),order)*math.factorial(order)/((2.0*order+a+b+1.0)*poch((a+1.0),order)*poch(-N,order)*math.factorial(N))


def doSum(n, o1, o2):
  global a, b, N
  sum = 0.0
  for x in range(0,n+1):
    sum = sum + nCr(a+x,x)*nCr(b+N-x,N-x)*fval(x,o1)*fval(x,o2)
  return sum



def main():
  global N
  order = 7
  for i in range(0,N+1):
    sum = doSum(i,order,order)
    print i, sum, (sum - exact(order))/exact(order)*100.0
  print 'Exact = '+str(exact(order))


#print fval(5.0,5)
#print nCr(a+5,5)*nCr(b+N-5,N-5)*fval(5.0,0)*fval(5.0,0)
main()
