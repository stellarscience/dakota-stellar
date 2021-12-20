#!/usr/bin/env python
from __future__ import print_function
import unittest
import os

__author__ = 'J. Adam Stephens'
__copyright__ = 'Copyright 2014-2020 National Technology & Engineering Solutions of Sandia, LLC (NTESS)'
__license__ = 'GNU Lesser General Public License'

try: # Python 2/3 compatible import of StringIO
    import StringIO
except ImportError:
    import io as StringIO
import dakota.interfacing as di
from dakota.interfacing import parallel

apreproParams = """                    { DAKOTA_VARS     =                      3 }
                    { x1              =  7.488318331306800e-01 }
                    { x2              =  2.188638686202466e-01 }
                    { dussv_1         =                      "foo bar" }
                    { DAKOTA_FNS      =                      1 }
                    { ASV_1:response_fn_1 =                      %d }
                    { DAKOTA_DER_VARS =                      2 }
                    { DVV_1:x1        =                      1 }
                    { DVV_2:x2        =                      2 }
                    { DAKOTA_AN_COMPS =                      2 }
                    { AC_1:first.sh   =                      "a b" }
                    { AC_2:first.sh   =                      "b" }
                    { DAKOTA_EVAL_ID  =                      1 }
"""

dakotaParams = """                                          3 variables
                      7.488318331306800e-01 x1
                      2.188638686202466e-01 x2
                                    foo bar dussv_1
                                          1 functions
                                          %d ASV_1:response_fn_1
                                          2 derivative_variables
                                          1 DVV_1:x1
                                          2 DVV_2:x2
                                          2 analysis_components
                                        a b AC_1:first.sh
                                          b AC_2:first.sh
                                          1 eval_id
"""

dakotaSingleEvalBatchParams = """                                          3 variables
                      7.488318331306800e-01 x1
                      2.188638686202466e-01 x2
                                    foo bar dussv_1
                                          1 functions
                                          1 ASV_1:response_fn_1
                                          2 derivative_variables
                                          1 DVV_1:x1
                                          2 DVV_2:x2
                                          0 analysis_components
                                        2:1 eval_id
"""


dakotaBatchParams = """                                          2 variables
                      2.109465794009156e-01 x1
                     -9.675715913879684e-01 x2
                                          1 functions
                                          1 ASV_1:obj_fn
                                          2 derivative_variables
                                          1 DVV_1:x1
                                          2 DVV_2:x2
                                          0 analysis_components
                                        1:1 eval_id
                                          2 variables
                      3.222614371054804e-01 x1
                      4.946014218730854e-02 x2
                                          1 functions
                                          1 ASV_1:obj_fn
                                          2 derivative_variables
                                          1 DVV_1:x1
                                          2 DVV_2:x2
                                          0 analysis_components
                                        1:2 eval_id
"""

# Helper functions needed for Python < 2.7
def set_function(r):
    r["response_fn_1"].function = 5.0

def set_gradient(r):
    r["response_fn_1"].gradient = [1.0189673084127668E-266, -6.3508646783183408E-264]

def set_hessian(r):
    r["response_fn_1"].hessian = [[1.0, 2.0],[2.0,3.0]]


class dakotaInterfacingTestCase(unittest.TestCase):

    def test_aprepro_format(self):
        """Confirm that aprepro format Parameters files are parsed correctly."""
        pio = StringIO.StringIO(apreproParams % 1)
        p, r = di.interfacing._read_parameters_stream(stream=pio, results_file="results.out")
        self.assertEqual(p.num_variables, 3)
        self.assertEqual(p.descriptors, ["x1","x2","dussv_1"])
        self.assertAlmostEqual(p["x1"], 7.488318331306800e-01)
        self.assertAlmostEqual(p["x2"], 2.188638686202466e-01)
        self.assertEqual(p["dussv_1"], "foo bar")
        self.assertEqual(p.num_an_comps, 2)
        self.assertEqual(p.an_comps[0], "a b")
        self.assertEqual(p.an_comps[1], "b")
        self.assertEqual(p.eval_id, "1")
        self.assertTrue(p.aprepro_format)
        
        self.assertEqual(r.num_responses, 1)
        self.assertEqual(r.descriptors, ["response_fn_1"])
        self.assertEqual(r["response_fn_1"].asv, (True, False, False))
        self.assertEqual(r.num_deriv_vars, 2)
        self.assertEqual(r.deriv_vars, ["x1","x2"])
        self.assertTrue(r.aprepro_format)
        self.assertEqual(r.results_file, "results.out")

    def test_dakota_format(self):
        """Confirm that Dakota format Parameters files are parsed correctly."""
        pio = StringIO.StringIO(dakotaParams % 1)
        p, r = di.interfacing._read_parameters_stream(stream=pio, results_file="results.out")
        self.assertEqual(p.num_variables, 3)
        self.assertEqual(p.descriptors, ["x1","x2","dussv_1"])
        self.assertAlmostEqual(p["x1"], 7.488318331306800e-01)
        self.assertAlmostEqual(p["x2"], 2.188638686202466e-01)
        self.assertEqual(p["dussv_1"], "foo bar")
        self.assertEqual(p.num_an_comps, 2)
        self.assertEqual(p.an_comps[0], "a b")
        self.assertEqual(p.an_comps[1], "b")
        self.assertEqual(p.eval_id, "1")
        self.assertFalse(p.aprepro_format)
        
        self.assertEqual(r.num_responses, 1)
        self.assertEqual(r.descriptors, ["response_fn_1"])
        self.assertEqual(r["response_fn_1"].asv, (True, False, False))
        self.assertEqual(r.num_deriv_vars, 2)
        self.assertEqual(r.deriv_vars, ["x1","x2"])
        self.assertEqual(r.eval_id, "1")
        self.assertFalse(r.aprepro_format)
        self.assertEqual(r.results_file, "results.out")

    def test_asv(self):
        """Results behaves according to the ASV when response data is set."""
        # Function only
        pio = StringIO.StringIO(dakotaParams % 1)
        p, r = di.interfacing._read_parameters_stream(stream=pio)
        set_function(r)
        self.assertRaises(di.interfacing.ResponseError, set_gradient, r)
        self.assertRaises(di.interfacing.ResponseError, set_hessian, r)
        r.write(StringIO.StringIO())
        # Gradient only
        pio = StringIO.StringIO(dakotaParams % 2)
        p, r = di.interfacing._read_parameters_stream(stream=pio)
        self.assertRaises(di.interfacing.ResponseError, set_function, r)
        set_gradient(r)
        self.assertRaises(di.interfacing.ResponseError, set_hessian, r)
        r.write(StringIO.StringIO())
        # Hessian only
        pio = StringIO.StringIO(dakotaParams % 4)
        p, r = di.interfacing._read_parameters_stream(stream=pio)
        self.assertRaises(di.interfacing.ResponseError, set_function, r)
        self.assertRaises(di.interfacing.ResponseError, set_gradient, r)
        set_hessian(r)
        r.write(StringIO.StringIO())
       
    def test_ignore_asv(self):
        """Confirm that ASV is ignored."""
        # Test exceptions
        for i in range(1,8):
            sio = StringIO.StringIO(dakotaParams % i)
            p, r = di.interfacing._read_parameters_stream(stream=sio,ignore_asv=True)
            set_function(r) 
            set_gradient(r) 
            set_hessian(r)
            r.write(StringIO.StringIO())
        # Test write-time ignoring
        sio = StringIO.StringIO(dakotaParams % 3)
        p, r = di.interfacing._read_parameters_stream(stream=sio,ignore_asv=False)
        set_function(r)
        rio = StringIO.StringIO()
        r.write(stream=rio,ignore_asv=True)
        self.assertEqual(rio.getvalue(), "  5.0000000000000000E+00 response_fn_1\n")
    
    def test_results_write(self):
        """Verify Written test results"""
        sio = StringIO.StringIO(dakotaParams % 7)
        p, r = di.interfacing._read_parameters_stream(stream=sio)
        set_function(r) 
        set_gradient(r) 
        set_hessian(r)
        rio = StringIO.StringIO()
        r.write(stream=rio)
        expected = """  5.0000000000000000E+00 response_fn_1
[   1.0189673084127668E-266 -6.3508646783183408E-264 ]
[[   1.0000000000000000E+00   2.0000000000000000E+00
     2.0000000000000000E+00   3.0000000000000000E+00 ]]
"""
        self.assertEqual(rio.getvalue(), expected)
        # Test simulation failure flag
        r.fail()
        rio = StringIO.StringIO()
        r.write(stream=rio)
        expected = "FAIL\n"
        self.assertEqual(rio.getvalue(), expected)

    def test_batch_support(self):
        """Verify that batch-format params and results files are handled correctly"""
       
        num_vars = 2
        num_fns = 1
        x1_vals = [2.109465794009156e-01, 3.222614371054804e-01]
        batch_id = 1
        eval_ids = ["1:1", "1:2"]
        sio = StringIO.StringIO(dakotaBatchParams)
        p, r = di.interfacing._read_parameters_stream(sio, False, True, None)
        # correct number of parameter sets and response sets
        self.assertEqual(len(p),2)
        self.assertEqual(len(r),2)
        # batch ids
        self.assertEqual(p.batch_id, "1") 
        self.assertEqual(r.batch_id, "1") 
        # Indexing
        for i in range(len(p)):
            self.assertAlmostEqual(p[i]["x1"], x1_vals[i])
            self.assertEqual(r[i].eval_id, eval_ids[i])
        with self.assertRaises(IndexError):
            p[2]
        with self.assertRaises(IndexError):
            r[2]
        
        # Test iteration
        for i, pi in enumerate(p):
             self.assertEqual(pi.eval_id, eval_ids[i])
        # Test write locking
        with self.assertRaises(di.interfacing.BatchWriteError):
            r[0].write("foo")
        # Verify writes
        r[0][0].function = 1.0
        r[1][0].function = 2.0
        results_file = StringIO.StringIO()
        r.write(results_file)
        results_strings = results_file.getvalue().split("\n")
        self.assertEqual(results_strings[0], "#")
        val, label = results_strings[1].split()
        self.assertAlmostEqual(float(val),1.0)
        self.assertEqual(results_strings[2], "#")
        val, label = results_strings[3].split()
        self.assertAlmostEqual(float(val),2.0)
    
    def test_single_eval_batch(self):
        """Verify that batch objects are returned for batch = True"""
        pio = StringIO.StringIO(dakotaSingleEvalBatchParams)
        p, r = di.interfacing._read_parameters_stream(stream=pio, 
                           batch=True, results_file="results.out")
        self.assertIsInstance(p, di.interfacing.BatchParameters)
        self.assertEqual(len(p), 1)
        self.assertIsInstance(r, di.interfacing.BatchResults)
        self.assertEqual(len(r), 1)

    def test_batch_setting_exception(self):
        """Ensure BatchSettingError is raised when batch=True"""
        pio = StringIO.StringIO(dakotaParams % 1)
        with self.assertRaises(di.interfacing.BatchSettingError):
            p, r = di.interfacing._read_parameters_stream(stream=pio, 
                           batch=True, results_file="results.out")

    def test_dprepro(self):
        """Verify that templates are substituted correctly"""
 
        sio = StringIO.StringIO(dakotaParams % 3)
        p, r = di.interfacing._read_parameters_stream(stream=sio)

        # Test insertion of DakotaParams and DakotaResults
        tpl = """{"%3.1f" % DakotaParams["x1"]}
{"%3.1f" % x1}
{DakotaResults[0].asv.function}
"""
        result = di.dprepro(tpl, parameters=p, results=r)
        solution = "0.7\n0.7\nTrue\n"
        self.assertEqual(result,solution) 
        # Test insertion of extra "include" parameters
        extra = {"foo":5}
        tpl = """{"%3.1f" % DakotaParams["x1"]}
{"%3.1f" % x1}
{DakotaResults[0].asv.function}
{foo}
"""  
        result = di.dprepro(tpl, parameters=p, results=r, include=extra)
        solution = "0.7\n0.7\nTrue\n5\n"
        self.assertEqual(result,solution)
        # Test without DakotaParams or DakotaResults
        tpl = "{foo}"
        result = di.dprepro(tpl,include=extra)
        solution = "5"
        self.assertEqual(result,solution)
 
    def test_slurm_info(self):
        """Info correctly extracted from the environment."""
        # Check env without SLURM_JOBID set
        try:
            del os.environ["SLURM_JOBID"]
        except KeyError:
            pass
        self.assertRaises(parallel.MgrEnvError,parallel._get_job_info)

        # Set JOBID and set TASKS_PER_NODE
        os.environ["SLURM_JOBID"] = "1337"
        os.environ["SLURM_TASKS_PER_NODE"] = "16(x2)"
        nnodes, tpn, job_id = parallel._get_job_info()
        self.assertEqual(nnodes, 2)
        self.assertEqual(tpn, 16)
        self.assertEqual(job_id,"1337")

        # Single node
        os.environ["SLURM_TASKS_PER_NODE"] = "16"
        nnodes, tpn, job_id = parallel._get_job_info()
        self.assertEqual(nnodes, 1)
        self.assertEqual(tpn, 16)

    def test_node_list(self):
        """Relative node list is built correctly."""

        self.assertEqual( "+n0",
                parallel._get_node_list(tile=0, tile_size=16,
                                             tasks_per_node=16))
        self.assertEqual( "+n0,+n1",
                parallel._get_node_list(tile=0, tile_size=32,
                                             tasks_per_node=16))
        self.assertEqual( "+n0",
                parallel._get_node_list(tile=0, tile_size=8,
                                             tasks_per_node=16))
        self.assertEqual( "+n0",
                parallel._get_node_list(tile=1, tile_size=8,
                                             tasks_per_node=16))
        self.assertEqual( "+n1",
                parallel._get_node_list(tile=0, tile_size=16,
                    tasks_per_node=16, dedicated_master='NODE'))

        self.assertEqual( "+n1,+n2",
                parallel._get_node_list(tile=0, tile_size=32,
                    tasks_per_node=16, dedicated_master='NODE'))

        self.assertEqual( "+n0",
                parallel._get_node_list(tile=0, tile_size=8,
                    tasks_per_node=16, dedicated_master='TILE'))

        self.assertEqual( "+n1",
                parallel._get_node_list(tile=1, tile_size=8,
                    tasks_per_node=16, dedicated_master='TILE'))

        self.assertEqual( "+n1",
                parallel._get_node_list(tile=0, tile_size=16,
                    tasks_per_node=16, dedicated_master='TILE'))



    def test_static_tile(self):
        """Correct tile is returned"""
        self.assertEqual(0,
                parallel._calc_static_tile(eval_num=1, num_tiles=4))
        self.assertEqual(3,
                parallel._calc_static_tile(eval_num=4, num_tiles=4))
        self.assertEqual(0,
                parallel._calc_static_tile(eval_num=5, num_tiles=4))
        self.assertEqual(0,
                parallel._calc_static_tile(eval_num=9, num_tiles=4))
        self.assertEqual(0,
                parallel._calc_static_tile(eval_num=1, num_tiles=1))
        self.assertEqual(0,
                parallel._calc_static_tile(eval_num=100, num_tiles=1))

    def test_num_tiles(self):
        self.assertEqual(10,
                parallel._calc_num_tiles(tile_size=16, tasks_per_node=16,
                    num_nodes=10))
        self.assertEqual(1,
                parallel._calc_num_tiles(tile_size=16, tasks_per_node=16,
                    num_nodes=1))
        self.assertEqual(9,
                parallel._calc_num_tiles(tile_size=16, tasks_per_node=16,
                    num_nodes=10, dedicated_master='NODE'))
        self.assertEqual(9,
                parallel._calc_num_tiles(tile_size=16, tasks_per_node=16,
                    num_nodes=10, dedicated_master='TILE'))
        self.assertEqual(18,
                parallel._calc_num_tiles(tile_size=8, tasks_per_node=16,
                    num_nodes=10, dedicated_master='NODE'))
        self.assertEqual(19,
                parallel._calc_num_tiles(tile_size=8, tasks_per_node=16,
                    num_nodes=10, dedicated_master='TILE'))
        self.assertEqual(1,
                parallel._calc_num_tiles(tile_size=8, tasks_per_node=16,
                    num_nodes=1, dedicated_master='TILE'))
        self.assertRaises(parallel.ResourceError,parallel._calc_num_tiles,
                    tile_size=8, tasks_per_node=16, num_nodes=1, 
                    dedicated_master='NODE')
        self.assertRaises(parallel.ResourceError,parallel._calc_num_tiles,
                    tile_size=16, tasks_per_node=16, num_nodes=1, 
                    dedicated_master='TILE')
# todo: test iteration, integer access



unittest.main()
