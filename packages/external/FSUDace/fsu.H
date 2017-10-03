char ch_cap ( char c );
bool ch_eqi ( char c1, char c2 );
int ch_to_digit ( char c );
double chi_measure ( int ndim, int n, double z[], int ns, int seed_init );
void cvt_iterate ( int ndim, int n, int batch, int sample, bool reset, 
  int sample_num, int *seed, double r[], double *energy );
void cvt_sample ( int ndim, int n, int n_now, int sample, bool reset, 
  int *seed, double r[] );
void cvt_write ( int ndim, int n, int batch, int seed_init, int seed, 
  char *init_string, int it_max, int it_num, char *sample_string, 
  int sample_num, double r[], char *file_out_name );
double d_epsilon ( void );
double d_huge ( void );
double d_max ( double x, double y );
double d_measure ( int ndim, int n, double z[], int ns, int seed_init );
double d_min ( double x, double y );
double dge_det ( int n, double a[] );
char digit_to_ch ( int i );
bool dmat_in_01 ( int m, int n, double a[] );
void dmat_transpose_print ( int m, int n, double a[], const char *title );
void dmat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, const char *title );
void dmat_uniform_01 ( int m, int n, int *seed, double r[] );
double *dtable_data_read ( const char *input_filename, int m, int n );
void dtable_data_write ( int m, int n, double table[], std::ofstream &output );
void dtable_header_read ( const char *input_filename, int *m, int *n );
void dtable_header_write ( int m, int n, char *output_filename, std::ofstream &output );
int *dvec_sort_heap_index_a ( int n, double a[] );
void dvec_uniform_01 ( int n, int *seed, double r[] );
int file_column_count ( const char *input_filename );
void file_name_ext_get ( const char *file_name, int *i, int *j );
char *file_name_ext_swap ( const char *file_name, const char *ext );
int file_row_count ( const char *input_filename );
void find_closest ( int ndim, int n, int sample_num, double s[], double r[],
  int nearest[] );
void fsu_cvt ( int ndim, int n, int batch, int init, int sample, int sample_num, 
  int it_max, int *seed, double r[], int *it_num );
void fsu_halton ( int ndim, int n, int step, int seed[], int leap[],
  int base[], double r[] );
void fsu_hammersley  ( int ndim, int n, int step, int seed[], int leap[],
  int base[], double r[] );
void fsu_latinize ( int m, int n, double table[] );
int get_seed ( void );
double h_measure ( int ndim, int n, double z[], int ns, int seed_init );
bool halham_leap_check ( int ndim, int leap[] );
bool halham_n_check ( int n );
bool halham_ndim_check ( int ndim );
bool halham_seed_check ( int ndim, int seed[] );
bool halham_step_check ( int step );
void halham_write ( int ndim, int n, int step, int seed[], int leap[], int base[], 
  double r[], char *file_out_name );
bool halton_base_check ( int ndim, int base[] );
bool hammersley_base_check ( int ndim, int base[] );
int i_log_10 ( int i );
int i_max ( int i1, int i2 );
int i_min ( int i1, int i2 );
char *i_to_s ( int i );
void ivec_transpose_print ( int n, int a[], const char *title );
double *pointset_spacing ( int ndim, int n, double z[] );
int prime ( int n );
unsigned long random_initialize ( int seed );
bool s_eqi ( const char *s1, const char *s2 );
int s_index_last_c ( const char *s, char c );
int s_len_trim ( const char *s );
double s_to_d ( char *s, int *lchar, bool *error );
bool s_to_dvec ( char *s, int n, double rvec[] );
double tau_measure ( int ndim, int n, double z[], int ns, int seed_init );
void timestamp ( void );
char *timestring ( void );
void tuple_next_fast ( int m, int n, int rank, int x[] );
int word_count ( char *s );
