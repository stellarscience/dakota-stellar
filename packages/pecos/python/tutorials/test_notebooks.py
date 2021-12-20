#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

import os, glob, subprocess, tempfile, nbformat

def run_notebook(path):
    """
    Execute a notebook via nbconvert and collect the output.

    Parameters
    ----------
    path : string
        The full filename (including path) of the notebook

    Returns
    -------
    notebook : notebook object
        The parsed notebook

    errors : list
        errors found whilst executing notebook
    """
    dirname, __ = os.path.split(path)
    #print dirname == None
    #os.chdir(dirname)
    with tempfile.NamedTemporaryFile(suffix=".ipynb") as fout:
        args = ["jupyter", "nbconvert", "--to", "notebook", "--execute",
          "--ExecutePreprocessor.timeout=60",
          "--ExecutePreprocessor.kernel_name=python",
          "--ExecutePreprocessor.allow_errors=True",
          "--output", fout.name, path]
        subprocess.call(args)

        fout.seek(0)
        notebook = nbformat.read(fout,nbformat.current_nbformat)

    errors = [output for cell in notebook.cells if "outputs" in cell
                     for output in cell["outputs"]\
                     if output.output_type == "error"]

    return notebook, errors

if __name__== '__main__':

    filenames = ['Multivariate Function Approximation With Dakota.ipynb']
    for filename in filenames:
        # Test if any Exceptions or assert violations were thrown
        notebook, errors = run_notebook(filename)
        assert errors == []
