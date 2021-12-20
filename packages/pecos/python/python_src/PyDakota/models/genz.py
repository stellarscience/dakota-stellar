#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

import numpy
from PyDakota.approximation import Function, PyFunction
class GenzFunction(Function):
    def __init__( self, func_type, num_vars, c=None, w=None, name=None ):
        """
        having c and w as arguments are required to enable pickling
        name allows name to be consistent when pickling
        because name set here is appended to when set_coefficients is called
        and that function is not called when pickling
        """
        Function.__init__( self ) # must call c++ object wrapper init function

        self.func_type = func_type
        self.num_vars = num_vars
        self.num_qoi = 1
        if name is None:
            self.name = 'genz-' + self.func_type + '-%d' %self.num_vars
        else:
            self.name=name
        self.c=c
        self.w=w

        self.evaluator=PyFunction(self.value_)

    def value(self, x, opts=dict()):
        eval_type=opts.get('eval_type','value')
        if ( ('grad' in eval_type) and
            (self.func_type == "discontinuous" or self.func_type == "continuous") ):
            msg =  "gradients cannot be computed for %s Genz function"%self.func_type
            raise Exception, msg
        assert x.min()>=0 and x.max()<=1.
        vals = self.evaluator.value(x)
        if eval_type=='value':
            return vals[:,:1]
        if eval_type=='value-grad':
            return vals
        if eval_type=='grad':
            return vals[:,1:]

    def value_( self, x ):
        assert x.ndim == 1
        if (self.func_type == "discontinuous" or self.func_type == "continuous"):
            result = numpy.empty( ( 1 ) , numpy.double )
        else:
            result = numpy.empty( ( self.num_vars+1 ) , numpy.double );
        self.num_vars = self.c.shape[0]
        if ( self.func_type == "oscillatory" ):
            result[0] = 2.0 * numpy.pi * self.w[0]
            for d in range( self.num_vars ):
                result[0] += self.c[d] * x[d];
            for d in range( self.num_vars ):
                result[1+d] = -self.c[d] * numpy.sin( result[0] );
            result[0] = numpy.cos( result[0] );
        elif ( self.func_type == "product-peak" ):
            result[0] = 1.0
            for d in range( self.num_vars ):
                result[0] *= ( 1.0 / ( self.c[d] * self.c[d] ) +
                               ( x[d] - self.w[d] ) *  ( x[d] - self.w[d] ) )
            for d in range( self.num_vars ):
                result[1+d] = 2. * ( x[d] - self.w[d] ) / result[0]
            result[0] = 1.0 / result[0]
        elif ( self.func_type == "corner-peak" ):
            result[0] = 1.0;
            for d in range( self.num_vars ):
                result[0] += self.c[d] * x[d];
            for d in range( self.num_vars ):
                result[1+d]  = -self.c[d] * (self.num_vars+1) / \
                    result[0]**( self.num_vars+2 )
            result[0]  = 1.0 / result[0]**( self.num_vars+1 )
        elif ( self.func_type == "gaussian-peak" ):
            result[0] = 0.0;
            for d in range( self.num_vars ):
                result[0] += self.c[d] * self.c[d] * ( x[d] - self.w[d] ) * \
                    ( x[d] - self.w[d] )
            for d in range( self.num_vars ):
                result[1+d] = 2. * self.c[d]**2 * (self.w[d]-x[d]) * \
                    numpy.exp( -result[0] )
            result[0] = numpy.exp( -result[0] );
        elif ( self.func_type == "continuous" ):
            result[0] = 0.0;
            for d in range( self.num_vars ):
                result[0] += self.c[d] * numpy.abs( x[d] - self.w[d] );
            result[0] = numpy.exp( -result[0] );
            result[1:] = 0.
        elif ( self.func_type == "discontinuous" ):
            result[0] = 0.0;
            if ( self.num_vars == 1 ):
                if ( x[0] <= self.w[0] ):
                    for d in range( self.num_vars ):
                            result[0] += self.c[d] * x[d];
                    result[0] = numpy.exp( result[0] );
            else:
                if ( ( x[0] <= self.w[0] ) and ( x[1] <= self.w[1] ) ):
                    for d in range( self.num_vars ):
                        result[0] += self.c[d] * x[d];
                    result[0] = numpy.exp( result[0] );
        elif (self.func_type=='corner-peak-second-order'):
            result[0]=0.
            for d in xrange(self.num_vars-1):
                result[0] += (1.+self.c[d]*x[d] +self.c[d+1]*x[d+1])**(-3)
        else:
            msg = "ensure func_num in [\"oscillatory\",\"product-peak\","
            msg += "\"corner-peak\",\"gaussian-peak\","
            msg += "\"continuous\",\"discontinuous\"]"
            raise Exception, msg

        return result

    def set_coefficients( self, c_factor, coef_type, w_factor = 0., seed=0 ):

        self.c = numpy.empty( ( self.num_vars ), numpy.double );
        self.w = numpy.empty( ( self.num_vars ), numpy.double );
        self.name += '-' + coef_type

        if ( coef_type == "no-decay" ):
            csum = 0.0
            for d in range( self.num_vars ):
                self.w[d] = w_factor
                self.c[d] = ( d + 0.5 ) / self.num_vars
                csum += self.c[d]

            self.c *= ( c_factor / csum )
        elif ( coef_type == "quadratic-decay" ):
            csum = 0.0
            for d in range( self.num_vars ):
                self.w[d] = w_factor
                self.c[d] = 1.0 / ( d + 1. )**2
                csum += self.c[d]
            for d in range( self.num_vars ):
                self.c[d] *= ( c_factor / csum )
        elif ( coef_type == "quartic-decay" ):
            csum = 0.0
            for d in range( self.num_vars ):
                self.w[d] = w_factor
                self.c[d] = 1.0 / ( d + 1. )**4
                csum += self.c[d]
            for d in range( self.num_vars ):
                self.c[d] *= ( c_factor / csum )
        elif ( coef_type == "squared-exponential-decay" ):
            csum = 0.0
            for d in range( self.num_vars ):
                self.w[d] = w_factor
                self.c[d]=numpy.exp((d+1)**2*numpy.log(1.e-15)/self.num_vars**2)
                csum += self.c[d]
            for d in range( self.num_vars ):
                self.c[d] *= ( c_factor / csum )
        elif ( coef_type == "exponential-decay" ):
            csum = 0.0
            for d in range( self.num_vars ):
                self.w[d] = w_factor
                #self.c[d] = numpy.exp( (d+1)*numpy.log( 1.e-15 )/self.num_vars )
                self.c[d] = numpy.exp( (d+1)*numpy.log( 1.e-8 )/self.num_vars )
                csum += self.c[d]
            for d in range( self.num_vars ):
                self.c[d] *= ( c_factor / csum )
        elif (coef_type=='random'):
            numpy.random.seed(seed)
            self.c = numpy.random.uniform(0.,1.,(self.num_vars))
            csum = numpy.sum(self.c)
            self.c *= c_factor/csum
            #self.w = numpy.random.uniform( 0.,1.,(self.num_vars))
            self.w = 0.
        else:
            msg = "Ensure coef_type in [\"no-decay\",\"quadratic-decay\""
            msg += ",\"quartic-decay\", \"squared-exponential-decay,\""
            msg += ",\"exponential-decay\"]"
            raise Exception, msg

    def variance( self ):
        mean = self.integrate()
        return self.recursive_uncentered_moment(self.num_vars,0.0,self.num_vars+1)-mean**2

    def recursive_uncentered_moment( self, d, integral_sum, r ):
        if ( self.func_type == "corner-peak" ):
            if ( d > 0 ):
                return ( self.recursive_uncentered_moment(
                    d-1, integral_sum, r ) -
                    self.recursive_uncentered_moment(
                        d-1, integral_sum + self.c[d-1], r ) ) / \
                        ( (d+r) * self.c[d-1] )
            else:
                return 1.0 / ( 1.0 + integral_sum )**(1+r);
        else:
            return 0.

    def integrate( self ):
        return self.recursive_integrate( self.num_vars, 0.0 )

    def recursive_integrate( self, d, integral_sum ):
        if ( self.func_type == "oscillatory" ):
            if ( d > 0 ):
                return ( self.recursive_integrate( d-1,
                                                   integral_sum+self.c[d-1] )
                         - self.recursive_integrate( d-1,
                                                     integral_sum ) ) /\
                                                     self.c[d-1];
            else:
                case = self.num_vars % 4
                if (case == 0 ):
                    return  numpy.cos( 2.0 * numpy.pi * self.w[0] + integral_sum)
                if (case == 1 ):
                    return  numpy.sin( 2.0 * numpy.pi * self.w[0] + integral_sum)
                if (case == 2 ):
                    return -numpy.cos( 2.0 * numpy.pi * self.w[0] + integral_sum)
                if (case == 3 ):
                    return -numpy.sin( 2.0 * numpy.pi * self.w[0] + integral_sum)

        elif ( self.func_type == "product-peak" ):
            prod = 1.0
            for i in xrange( self.num_vars ):
                prod = prod * self.c[i] * ( numpy.arctan( self.c[i] *
                                                       (1.0 - self.w[i] ) ) +
                                            numpy.arctan( self.c[i] * self.w[i]))
            return prod
        elif ( self.func_type == "corner-peak" ):
            if ( d > 0 ):
                return ( self.recursive_integrate( d-1, integral_sum ) -
                         self.recursive_integrate( d-1, integral_sum +
                                                   self.c[d-1] ) ) / \
                                                   ( d * self.c[d-1] )
            else:
                return 1.0 / ( 1.0 + integral_sum );

        elif ( self.func_type == "gaussian-peak" ):
            msg = "GenzModel::recursive_integrate() integration of "
            msg += "gaussian_peak function  has not been implmemented.";
            raise Exception, msg

        elif ( self.func_type == "continuous" ):
            prod = 1.0;
            for i in xrange( self.num_vars ):
                prod /= (  self.c[i] * ( 2.0 - numpy.exp( -self.c[i]*self.w[i])-
                                         numpy.exp( self.c[i]*(self.w[i]-1.0))));
            return prod;
        elif ( self.func_type == "discontinuous" ):
            prod = 1.0;
            if ( self.num_vars < 2 ):
                for i in xrange( self.num_vars ):
                    prod *= ( numpy.exp( self.c[i] * self.w[i] )-1.0 )/self.c[i];
            else:
                for i in xrange( 2 ):
                    prod *= ( numpy.exp( self.c[i] * self.w[i] ) - 1.0)/self.c[i]
                for i in xrange( 2, self.num_vars ):
                    prod *= ( numpy.exp( self.c[i] ) - 1.0 ) / self.c[i];
            return prod;
        elif (self.func_type=='corner-peak-second-order'):
            self.func_type='corner-peak'
            c = self.c.copy()
            num_vars = self.num_vars
            result = 0
            for d in xrange(num_vars-1):
                self.num_vars=2; self.c = c[d:d+2]
                result += self.recursive_integrate(2,0.)
                self.func_type='corner-peak-second-order'
            self.c = c.copy()
            self.num_vars = num_vars
            return result
        else:
            msg = "GenzModel::recursive_integrate() incorrect f_type_ ";
            msg += "was provided";
            raise Exception, msg

    def __reduce__(self):
        return (type(self),(self.func_type, self.num_vars, self.c, self.w, self.name))

# add unit test like to test pickling
#python -c "from PyDakota.models.genz import GenzFunction; g = GenzFunction('oscillatory',2); import pickle; pickle.dump(g,open('function.pkl','wb')); g1 = pickle.load(open('function.pkl','rb'))"
