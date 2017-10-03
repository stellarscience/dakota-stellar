/*-----------------------------------------------------*/
/*  how to use the NOMAD library with a user function  */
/*-----------------------------------------------------*/
#include "nomad.hpp"
#include <pthread.h>
#include <semaphore.h>

using namespace std;
// using namespace NOMAD; avoids putting NOMAD:: everywhere

// Number of threads to be used in parallel 
#define NUM_THREADS    4 

// A semaphore to manage parallel evaluations with NUM_THREADS
sem_t *mySemaphore;

// A structure to pass arguments to the evaluation wrapper function
class My_Evaluator;
typedef struct Arg_Eval_tag {
	NOMAD::Eval_Point	*	x	;
	NOMAD::Double			h_max;
	bool                *	count_eval;
	const My_Evaluator	*	ev;     
} Arg_Eval_t;



class My_Evaluator : public NOMAD::Evaluator {
public:
	My_Evaluator  ( const NOMAD::Parameters & p ) :
    NOMAD::Evaluator ( p ) {}
	
	~My_Evaluator ( void ) {}
	
	/*----------------------------------------*/
	/*     The problem evaluatiom             */
	/*----------------------------------------*/	
	// Provide the objective and constraints here		
	bool eval_x ( NOMAD::Eval_Point   & x          ,
				 const NOMAD::Double & h_max      ,
				 bool                & count_eval   ) const
	{
		
		NOMAD::Double c1 = 0.0 , c2 = 0.0;
		for ( int i = 0 ; i < 5 ; i++ ) 
		{
			c1 += (x[i]-1).pow2();
			c2 += (x[i]+1).pow2();
		}
		x.set_bb_output  ( 0 , x[4]  ); // objective value
		x.set_bb_output  ( 1 , c1-25 ); // constraint 1
		x.set_bb_output  ( 2 , 25-c2 ); // constraint 2
		
		count_eval = true; // count a black-box evaluation
		
		return true;       // the evaluation succeeded
	}
	
	// Wrapper of eval_x used for parallel evaluation (pthreads). 	
	static void * wrapper_eval_x ( void * dummy_eval_arg   )
	{
		
		Arg_Eval_t * eval_arg = static_cast<Arg_Eval_t *>(dummy_eval_arg);
		
		NOMAD::Eval_Point & x=*(eval_arg->x);
		bool & count_eval = *(eval_arg->count_eval);
		const NOMAD::Double & h_max = (eval_arg->h_max);		
		const My_Evaluator *ev=(eval_arg->ev);
		
		ev->eval_x(x,h_max,count_eval);
		
		pthread_exit(NULL);
		
	}
	
	// Implementation
	bool eval_x ( std::list<NOMAD::Eval_Point *>	&list_x  ,
				 const NOMAD::Double				& h_max,
				 std::list<bool>					& list_count_eval ) const
	{
		int rc;
		

		
		list_count_eval.assign(list_count_eval.size(),false);  // Evaluations are not counted until eval_x is called and set count_eval
		
		// All threads are created
		pthread_t threads[list_x.size()];
		
		// The semaphore will allow to run NUM_THREADS in parallel
		mySemaphore=sem_open("sem",0,NUM_THREADS);
		
		// Arguments passed to the evaluation wrapper
		Arg_Eval_t * eval_arg=new Arg_Eval_t[list_x.size()];		

		std::list<NOMAD::Eval_Point *>::iterator it_x=list_x.begin(),end_x=list_x.end();
		std::list<bool>::iterator it_count=list_count_eval.begin();		
		int i=0;
		// The list of points is evaluated under the control of a semaphore
		for (it_x=list_x.begin(); it_x!=end_x; ++it_x,++it_count,++i)
		{
			eval_arg[i].x=(*it_x);
			eval_arg[i].h_max=h_max.value();
			eval_arg[i].count_eval=&(*it_count);
			eval_arg[i].ev=this;
			rc=pthread_create(&threads[i], NULL, wrapper_eval_x,&eval_arg[i]);
			if (rc)
			{
				cout << "Error:unable to create thread," << rc << endl;
				return false;
			}
			// wait until value of semaphore is greater than 0
			// decrement the value of semaphore by 1
			sem_wait(mySemaphore);
			
		}
		
		delete []eval_arg;
		int ret;
		for (i=0; i<list_x.size(); ++i)
		{
			// Wait for all the threads to finish
			ret=pthread_join(threads[i],0);
			if (ret!=0)
			{
				perror("pthread join has failed");
				return false;
			}
		}
		
		sem_unlink("sem");
		
		return true;       // the evaluation succeeded
	}
	
	
	
};

/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main ( int argc , char ** argv ) {

  // display:
  NOMAD::Display out ( std::cout );
  out.precision ( NOMAD::DISPLAY_PRECISION_STD );

  try {

    // NOMAD initializations:
    NOMAD::begin ( argc , argv );

    // parameters creation:
    NOMAD::Parameters p ( out );

    p.set_DIMENSION (5);             // number of variables

    vector<NOMAD::bb_output_type> bbot (3); // definition of
    bbot[0] = NOMAD::OBJ;                   // output types
    bbot[1] = NOMAD::PB;
    bbot[2] = NOMAD::EB;
    p.set_BB_OUTPUT_TYPE ( bbot );

//    p.set_DISPLAY_ALL_EVAL(true);   // displays all evaluations.
    p.set_DISPLAY_DEGREE(2);  
    p.set_DISPLAY_STATS ( "bbe ( sol ) obj" );

    p.set_X0 ( NOMAD::Point(5,0.0) );  // starting point

    p.set_LOWER_BOUND ( NOMAD::Point ( 5 , -6.0 ) ); // all var. >= -6
    NOMAD::Point ub ( 5 );                    // x_4 and x_5 have no bounds
    ub[0] = 5.0;                              // x_1 <= 5
    ub[1] = 6.0;                              // x_2 <= 6
    ub[2] = 7.0;                              // x_3 <= 7
    p.set_UPPER_BOUND ( ub );

    p.set_MAX_BB_EVAL (100);     // the algorithm terminates after
                                 // 100 black-box evaluations
	  
    // Max number of points to be given as a block for evaluation
	// This option is required to perform parallel evaluations 
	p.set_BB_MAX_BLOCK_SIZE(NUM_THREADS);
  
    // parameters validation:
    p.check();

    // custom evaluator creation:
    My_Evaluator ev   ( p );

    // algorithm creation and execution:
    NOMAD::Mads mads ( p , &ev );
    mads.run();
  }
  catch ( exception & e ) {
    cerr << "\nNOMAD has been interrupted (" << e.what() << ")\n\n";
  }

  NOMAD::Slave::stop_slaves ( out );
  NOMAD::end();

  return EXIT_SUCCESS;
}
