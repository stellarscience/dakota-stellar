#!/projects/vuq/anaconda/bin/python
import math

c = 0.10
B = 1.5

def fval(x, order):
  global c, B

  if order == 0:
    return 1.0

  elif order == 1:
    return (c*B+(c-1.0)*x)/(c*B)

  elif order == 2:
    return (c*c*B*(B+1.0) + (2.0*c*(B+1.0)-c+1.0)*(c-1.0)*x + (c-1)*(c-1)*x*x)/(c*c*B*(B+1.0))

  elif order > 2:
    om1 = order-1
    fm2 = fval(x, order-2)
    fm1 = fval(x, order-1)
    return ((om1+(om1+B)*c+(c-1.0)*x)*fm1 - om1*fm2)/(c*(om1+B))

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
  global c, B
  return math.factorial(order)/(poch(B,order)*math.pow(c,order)*math.pow(1.0-c,B))


def doSum(n, o1, o2):
  global c, B
  sum = 0.0
  for x in range(0,n):
    sum = sum + poch(B,x)/math.factorial(x)*math.pow(c,x)*fval(x,o1)*fval(x,o2)
  return sum



def main():

  order = 10
  for i in range(0,50):
    sum = doSum(i,order,order)
    print i, sum, (sum - exact(order))/exact(order)*100.0
  print 'Exact = '+str(exact(order))
  print 'Test value = '+str(fval(4.0,1))


main()
