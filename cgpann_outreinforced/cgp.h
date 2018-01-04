/* cgp.h
Julian F. Miller (c), 2009
version 1.1 of first public release on 20-July-2009
Dept. of Electronics, University of York, UK
*/

#include <stdio.h>

#ifndef __CGP_H__

#define __CGP_H__


#define DATA_IS_DOUBLE



/* this if statement defines the type of data that
is processed by cgp genotypes
if you have a user defined data type
that is not on this list then you need to add a case
to the if statement
Note this will mean that you will have to define a new function set
this is defined in the c program file node_function.c
*/

/* this defines how the maximum constant inputs are input to cgp */

#ifdef DATA_IS_DOUBLE
typedef double			data_type;
#endif

/* this global defines excatly how many constant inputs are used
in cgp. This varies for the data_type chosen. For unsigned int
data type num_constant_inputs is set to zero */

#define PI							3.1415926
#define MAX_NUM_ROWS				10
#define MAX_NUM_COLS				500
#define MAX_NUM_INPUTS				10
#define MAX_NUM_OUTPUTS				4
#define MAX_NUM_DATA				300
#define MAX_ARITY					MAX_NUM_INPUTS
#define MAX_NUM_NODES				MAX_NUM_ROWS*MAX_NUM_COLS
#define MAX_NUM_NODES_PLUS_INPUTS	MAX_NUM_NODES+MAX_NUM_INPUTS
#define MAX_NUM_GENES_PER_NODE		((MAX_ARITY * 2) + 1)
#define MAX_NUM_GENES				MAX_NUM_GENES_PER_NODE*MAX_NUM_NODES
#define MAX_OUTPUT_SIZE				MAX_NUM_INPUTS + MAX_NUM_NODES + MAX_NUM_OUTPUTS
#define MAX_NUM_RUNS				200
#define ERROR_THRESHOLD				0.0001
#define NUM_CLASSES					3
#define outs						4
//
//#define HITS_BASED_FITNESS


#define MAX_NUM_CHROMOSOMES			100  /* max population size */
#define MAX_NUM_LETTERS				100
#define MAX_NUM_FUNCTIONS			2
#define MAXNUM						4294967295
#define MAX_NUM_GENERATIONS			100000

/* these are all read from the .par file and never change after that */
int				population_size;
double			per_cent_mutate;
int				num_generations;
int				num_runs_total;
int				num_rows;
int				num_cols;
int				levels_back;
int				progress_report;
int				report_interval;
unsigned		global_seed;
int				save_best_chrom;
int				run_from_chrom;
int				shrink_phenotype;
int				arity;

/*  global constants calculated in get_parameters */
int				num_functions;
int				num_genes;
int				num_nodes, num_genes_per_node;
int				number[MAX_NUM_FUNCTIONS];
char			node_types[MAX_NUM_FUNCTIONS][20];



/* this stores the node function address in
allowed_functions[] */
int             allowed_functions[MAX_NUM_FUNCTIONS];

double			fitarray[MAX_NUM_GENERATIONS];// Gulrukh

/* data defining the computational problem read from .dat file */
int             num_inputs, num_outputs, num_tests;
data_type	    data_inputs[MAX_NUM_DATA][MAX_NUM_INPUTS];
data_type	    data_outputs[MAX_NUM_DATA][MAX_NUM_INPUTS];
data_type	    cgpout[MAX_NUM_DATA][MAX_NUM_OUTPUTS];
data_type	    cgpout2[MAX_NUM_DATA][MAX_NUM_OUTPUTS];//Gulrukh
/* calculated global constants  */
int				perfect, bit_width;


/* macros */
#define pow2(x) (1<<x)                                   /* returns 2 to power x (x>=0) */
#define getbit(decimal,nthbit) ((decimal>>nthbit) & 01)  /* gets nth bit */



/* function prototypes */

void validate_command_line(int argc, char* argv[], char parfile[], char datafile[]);

void get_parameters(char parfile[MAX_NUM_LETTERS], char datafile[MAX_NUM_LETTERS]);

void write_cgp_info(char command[], char datafile[MAX_NUM_LETTERS]);

int newrand(int range);

data_type myfscanf(FILE* fp);

void read_data(char datafile[MAX_NUM_LETTERS]);

void define_perfect(void);

void fprint_a_chromosome(int* chromosome, char name[], int append);

void print_a_chromosome(int* chromosome);

void fprint_a_raw_chromosome(int* chromosome, char name[], int append);

void fprint_node_used(int size, int array[MAX_NUM_NODES_PLUS_INPUTS], char name[], int append);

int get_nodes_to_process(int* chromosome, int nodes_to_process[MAX_NUM_NODES]);

void fprint_active_genes(int* chromosome, char name[30]);

void read_from_chrom(int** chromosomes);

data_type  node_type(data_type in[MAX_NUM_GENES_PER_NODE],
	int function_gene, int arity);

void decode_cgp(int* chromosome,
	data_type data_inputs[MAX_NUM_DATA][MAX_NUM_INPUTS],
	data_type cgp_outputs[MAX_NUM_OUTPUTS],
	int num_nodes_to_process,
	int nodes_to_process[MAX_NUM_NODES],
	int fitness_test);


double correctness_test_boolean(data_type data_output, data_type cgp_output);

double correctness_test2(data_type data_output, data_type cgp_output, data_type change);

double evaluate_cgp_outputs(data_type cgp_outputs[MAX_NUM_OUTPUTS],
	int test);

double fitness(int* chromosome);

void get_gene_limits(int column, int* limit_min, int* limit);

int get_connection_gene(int limit_min, int limit);

int get_function_gene(void);

int get_output_gene(void);

int get_weight(void);

int is_not_output_gene(int gene);

int is_function_gene(int gene, int locus);

void generate_a_random_chromosome(int* chromosome);

void initialise(int** chromosomes);

double  get_best_chromosome(int** chromosomes,
	int*  best_chromosome,
	double previous_best_fitness,
	int    gen);

int get_num_mutant(int num_genes, double per_cent_mutate);


void mutate_a_gene(int*  chromosome);

void mutate_a_chromosome(int *chromosome, int num_mutations);

void generate_new_pop_es(int** chromosomes,
	int*  best_chromosome);

int** create_2d_space(int num_horizontals, int num_verticals);

int* create_1d_space(int size);

void free_array2d(int num_verticals, int** array2d);

void write_generation_to_screen(int gen_index);

void write_progress_info_to_screen(int generation, double fit);

void write_result_of_EA_to_file(int run, int bgen, double best_fit,
	int* best_chromosome);

void check_if_improvement(double best_fit, double* previous_best_fit, int* best_gen, int gen,
	char prog[MAX_NUM_LETTERS],
	int* best_chromosome);

double  EA(int *gen_of_best, int* num_nodes_active, int run,
	char prog[MAX_NUM_LETTERS]);

void setup_report_files(int run, char prog[MAX_NUM_LETTERS]);

void report_final_results(double av_fitness, double st_dev,
	double av_num_nodes_active,
	double best_of_best_fit,
	double worst_of_best_fit,
	double av_best_gen,
	int    num_perfect);

double average(int num_items, double items[MAX_NUM_RUNS]);

double get_standard_deviation(int num_items, double average,
	double items[MAX_NUM_RUNS]);
void get_cgp_outputs(data_type cgpout[MAX_NUM_DATA][MAX_NUM_OUTPUTS]);//Gulrukh

void run_EA(int num_runs_total);

void run_test(char testfile[MAX_NUM_LETTERS], char outfile[MAX_NUM_LETTERS]);

void get_cgp_outputs(int* best_chromosome, char file[MAX_NUM_LETTERS]);

#endif
