#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

from PyDakota.options_list import OptionsList
from PyDakota.swig_examples.options_list_interface import *
def test_options_list_typemap_in():
    """
    Test that a wrapped function that takes an OptionsList as input
    also accepts a PythonDictionary
    """
    opts = OptionsList()
    str_types = ['int','double','string','optionslist']
    items = [1,2.,'a',OptionsList()]
    names = ['key%s'%(i+1) for i in range(len(items))]
    for i,type_str in enumerate(str_types):
        item = items[i]
        name = names[i]
        set_entry = PyDakota.swig_examples.options_list_interface.__dict__[type_str + "set_entry"]
        tmp_opts = set_entry({},name,item)
        opts = set_entry(opts,name,item)
        assert tmp_opts=={names[i]:items[i]}

    pydict = dict((names[i],items[i]) for i in range(len(items)))
    assert opts == pydict

    opts = OptionsList()
    list_of_opts = []
    opts1 = OptionsList(); opts1.set(names[0],items[0])
    list_of_opts.append(opts1)
    opts2 = OptionsList(); opts2.set(names[1],items[1])
    list_of_opts.append(opts2)
    opts.set(names[3],list_of_opts)
    #tmp_opts = intset_entry(opts,names[2],1)
    #assert tmp_opts=={names[2]:1,names[3]:[{names[0]:items[0]},{names[1]:items[1]}]}
    print opts,{names[3]:[{names[0]:items[0]},{names[1]:items[1]}]}
    assert opts=={names[3]:[{names[0]:items[0]},{names[1]:items[1]}]}

if __name__ == "__main__":
    test_options_list_typemap_in()
