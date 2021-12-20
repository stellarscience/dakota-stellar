#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

# System imports
import numpy
import sys
import unittest

from PyDakota import options_list
from PyDakota.options_list import *

class OptionsListTestCase(unittest.TestCase):
    "TestCase class for Teuchos.ParameterList"

    def test_set(self):
        opts = OptionsList()
        opts.set("key1",2.)
        assert opts.get("key1")==2.
        self.assertRaises(KeyError, opts.get,"key2")

        opts.set("key1",2)
        assert opts.get("key1")==2

        opts.set("key1","a")
        assert opts.get("key1")=="a"

        opts.set("key1",{})
        assert opts.get("key1")=={}

        array = numpy.array([1,2,3],dtype=numpy.int64)
        opts.set("key3",array)
        assert numpy.allclose(opts.get("key3"),array)

        array = numpy.array([1,2,3],dtype=numpy.int32)
        opts.set("key3",array)
        assert numpy.allclose(opts.get("key3"),array)

        array = numpy.array([1,2,3],dtype=numpy.double)
        opts.set("key3",array)
        assert numpy.allclose(opts.get("key3"),array)

        array = [1,2,3]
        opts.set("key3",array)
        assert numpy.allclose(opts.get("key3"),array)

        array = [1.,2.,3.]
        opts.set("key3",array)
        assert numpy.allclose(opts.get("key3"),array)

        array = numpy.arange(10,dtype=numpy.double).reshape(2,5)
        opts.set("key4",array)
        assert numpy.allclose(opts.get("key4"),array)

        array = [[1,2,3,4],[5,6,7,8]]
        opts.set("key4",array)
        assert numpy.allclose(opts.get("key4"),array)

        list_of_opts = []
        opts1 = {"key1":2.}
        list_of_opts.append(opts1)
        opts2 = {"key2":3.}
        list_of_opts.append(opts2)
        opts.set("key5",list_of_opts)

    def test_len(self):
        opts = OptionsList()
        assert len(opts)==0

        opts = OptionsList({'key1':'a','key2':2})
        assert len(opts)==2

    def test_pydict_to_options_list(self):
        opts = OptionsList()
        opts.set("key1",{})
        assert len(opts)==1

    def test__eq__(self):
        opts1 = OptionsList()
        opts1.set("key1",{})
        array = [1.,2.,3.]
        opts1.set("key2",array)

        opts2 = OptionsList()
        opts2.set("key1",{})
        opts2.set("key2",array)

        assert opts1==opts2

        assert opts1==opts1

        opts2.set("key2","a")
        assert opts1 != opts2

        items = [1,2.,'a',OptionsList()]
        names = ['key%s'%(i+1) for i in range(len(items))]
        opts = OptionsList()
        list_of_opts = []
        opts1 = OptionsList(); opts1.set(names[0],items[0])
        list_of_opts.append(opts1)
        opts2 = OptionsList(); opts2.set(names[1],items[1])
        list_of_opts.append(opts2)
        opts.set(names[3],list_of_opts)
        assert opts=={names[3]:[{names[0]:items[0]},{names[1]:items[1]}]}

    def test__get_item__(self):
        opts = OptionsList()
        opts.set("key1",2.)
        assert opts["key1"]==2.

        items = [1,2.,'a',OptionsList()]
        names = ['key%s'%(i+1) for i in range(len(items))]
        opts = OptionsList()
        list_of_opts = []
        opts1 = OptionsList(); opts1.set(names[0],items[0])
        list_of_opts.append(opts1)
        opts2 = OptionsList(); opts2.set(names[1],items[1])
        list_of_opts.append(opts2)
        opts.set(names[3],list_of_opts)
        assert opts[names[3]]==[{names[0]:items[0]},{names[1]:items[1]}]

    def test__contains__(self):
        opts = OptionsList()
        opts.set("key1",2.)
        assert "key1" in opts

    def test__set_item__(self):
        opts = OptionsList()
        opts["key1"]=2.
        opts["key1"]==2.

    def test__str__(self):
        opts = OptionsList()
        opts.set("key1",[1,2,3])
        #print opts
        #print opts.__repr__

        items = [1,2.,'a',OptionsList()]
        names = ['key%s'%(i+1) for i in range(len(items))]
        opts = OptionsList()
        list_of_opts = []
        opts1 = OptionsList(); opts1.set(names[0],items[0])
        list_of_opts.append(opts1)
        opts2 = OptionsList(); opts2.set(names[1],items[1])
        list_of_opts.append(opts2)
        opts.set(names[3],list_of_opts)
        #print opts


if __name__ == "__main__":

    # Create the test suite object
    suite = unittest.TestSuite()

    # Add the test cases to the test suite
    suite.addTest(unittest.makeSuite(OptionsListTestCase))

    # Run the test suite
    result = unittest.TextTestRunner(verbosity=1).run(suite)

    
