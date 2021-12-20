/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014 Sandia Corporation.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:       Demo_Opt
//- Description: Demo TPL initialize
//- Owner:       Russell Hooper
//- Checked by:  ...

/** ... description ...
 * */

#ifndef DEMO_INITIALIZE_H
#define DEMO_INITIALIZE_H

#include <string>
#include <map>
#include <memory>
#include <vector>

class Demo_Opt
{
  public:

    // API (Pure virtual) callback for objective Fn (only one for now)
    class ObjectiveFn
    {
      public:
        virtual double compute_obj(const std::vector<double> & x, bool verbose = false) = 0;
    };

    // API (Pure virtual) callback for Nonlinear Equality Constaints
    class NonlinearEqFn
    {
      public:
        virtual int get_num_nln_eq(bool verbose = false) = 0;
        virtual void compute_nln_eq(std::vector<double> &c, const std::vector<double> &x, bool verbose = false) = 0;
    };

    // API (Pure virtual) callback for Nonlinear Inequality Constaints
    class NonlinearIneqFn
    {
      public:
        virtual int get_num_nln_ineq(bool verbose = false) = 0;
        virtual void compute_nln_ineq(std::vector<double> &c, const std::vector<double> &x, bool verbose = false) = 0;
    };

    // Default ctor
    Demo_Opt();

    // Allow specification of options filename
    bool set_solver_options(const std::string & filename, bool verbose = false);

    // A simple initialization method
    bool initialize(bool verbose = false);

    // A simple initialization method
    bool execute(bool verbose = false);

    // Parameter settings
    template< typename T >
      void set_param(const std::string & param, const T & val)
        { set_parameter_value(param, val); }

    // Register an objective fn callback interface
    void register_obj_fn(ObjectiveFn* obj_fn)
      { obj_fn_callback_ = obj_fn; }

    // Register an nonlinear equality fn callback interface
    void register_nln_eq_fn(NonlinearEqFn* obj)
      { nln_eq_fn_callback_ = obj; }

    // Register an nonlinear inequality fn callback interface
    void register_nln_ineq_fn(NonlinearIneqFn* obj)
      { nln_ineq_fn_callback_ = obj; }

    // Specify problem data
    void set_problem_data(const std::vector<double> &,   //  "Initial Guess"
                          const std::vector<double> &,   //  "Lower Bounds"
                          const std::vector<double> & ); //  "Upper Bounds"

    // Get best current state
    const std::vector<double> & get_best_x() const
      { return best_x_;}

    // Get best current objective values
    double get_best_f() const
      { return best_f_;}

    // Get best nonlinear equality constraint values
    const std::vector<double> & get_best_nln_eqs() const
      { return best_nln_eqs_;}

    // Get best nonlinear equality constraint values
    const std::vector<double> & get_best_nln_ineqs() const
      { return best_nln_ineqs_;}


  private:

    void set_parameter_value(const std::string & param, const int & val)
      { int_params_[param] = val; }
    void set_parameter_value(const std::string & param, const double & val)
      { dbl_params_[param] = val; }

    std::string options_file_;

    std::map<std::string, int>    int_params_;
    std::map<std::string, double> dbl_params_;

    ObjectiveFn     * obj_fn_callback_;
    NonlinearEqFn   * nln_eq_fn_callback_;
    NonlinearIneqFn * nln_ineq_fn_callback_;

    std::vector<double> init_vals_;
    std::vector<double> lower_bnds_;
    std::vector<double> upper_bnds_;

    std::vector<double> best_x_;
    double best_f_;
    std::vector<double> best_nln_eqs_;
    std::vector<double> best_nln_ineqs_;
};

#endif
