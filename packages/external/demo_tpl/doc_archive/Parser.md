
## Passing Options

 Dakota maintains a master list of hierarchical options in its
 $DAKTOA_SRC/src/dakota.xml file.  Several common options associated
 with optimizers are already supported and then only need to be
 exposed within the correct hierarchy (scope).  Initially, however, it
 is often best to simply pass a native options file directly to a TPL.
 Dakaota supprts this approach via the `options_file` specification as
 shown in the test input file.  To both expose the `demo_tpl`
 optimizer and associate the `options_file` specification with it, the
 dakota.xml file would be modified as follows:

 ```
   # File $DAKTOA_SRC/src/dakota.xml

   ... snip ...
    <!-- **** TOPLEVEL *** -->
    <keyword id="method" name="method" minOccurs="1" maxOccurs="unbounded" code="{N_mdm3(start,0,stop)}" label="Method" >

     ... snip ...
      <!-- Primary method selection alternation -->
      <oneOf label="Method (Iterative Algorithm)">

       ... snip ...

        <keyword  id="demo_tpl" name="demo_tpl" code="{N_mdm(utype,methodName_DEMO_TPL)}" label="demo_tpl" help="" minOccurs="1" group="Optimization: Local" >
          <keyword  id="options_file" name="options_file" code="{N_mdm(str,advancedOptionsFilename)}" label="Advanced Options File"  minOccurs="0" default="no advanced options file" >
            <param type="INPUT_FILE" />
          </keyword>
        </keyword>

       ... end snip ...
     ... end snip ...
   ... end snip ...
 ```

 Dakota's current parser system next needs to connect this change to
 it's internal options database and to its list of methods.  This is
 accomplished by modifying a few files as follows, eg

 ```
   # File $DAKTOA_SRC/src/NIDRProblemDescDB.cpp
   <... snip ...>
     MP2s(methodName,ROL),      // existing method
     MP2s(methodName,DEMO_TPL), // -----  our new demo_tpl method -----
     MP2s(methodName,NL2SOL),   // existing method
   <... end snip ...>
 ```

 ```
   # File $DAKTOA_SRC/src/DataMethod.hpp
   <... snip ...>
       GENIE_OPT_DARTS, GENIE_DIRECT,
       // Place Demo Opt TPL here based on current state of non-gradient flavor
       DEMO_TPL,                // -----  our new demo_tpl method -----
       // Gradient-based Optimizers / Minimizers:
       NONLINEAR_CG, OPTPP_CG, OPTPP_Q_NEWTON, OPTPP_FD_NEWTON, OPTPP_NEWTON,
   <... end snip ...>
 ```


 ```
   # File $DAKTOA_SRC/src/DataIterator.cpp
   <... snip ...>
        #ifdef HAVE_DEMO_TPL
        #include "DemoOptimizer.hpp"
        #endif
   <... end snip ...>

        Iterator* Iterator::get_iterator(ProblemDescDB& problem_db, Model& model)
        {
          unsigned short method_name = problem_db.get_ushort("method.algorithm");
   <... snip ...>
        #ifdef HAVE_DEMO_TPL
            case DEMO_TPL:      // -----  our new demo_tpl method -----
              return new DemoTPLOptimizer(problem_db, model); break;
        #endif
            default:
              switch (method_name) {
   <... end snip ...>


        /// bimap between method enums and strings; only used in this
        /// compilation unit
        static UShortStrBimap method_map =
          boost::assign::list_of<UShortStrBimap::relation>
          (HYBRID,                          "hybrid")
   <... snip ...>
          (DEMO_TPL,                        "demo_tpl")
          ;
   <... end snip ...>
 ```

 The next time Dakota is configured with the option `-D
 ENABLE_SPEC_MAINT:BOOL=ON` defined Dakota will automatically generate
 a file, $DAKTOA_SRC/src/dakota.input.nspec, based on the dakota.xml
 file.

 Once Dakota has been compiled with these changes, the simple test
 input file should parse and attempt to call the
 DemoTPLOptimizer::core_run() method to perform the optimization of
 the Dakota "text_book" example problem.
