/* cgp-functions.c
Julian F. Miller (c) 2009
version 1.1 of first public release on 20-July-2009
Dept. of Electronics, University of York, UK

IMPORTANT: For Boolean unsigned problems
program outputs are arranged most significant on the left
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cgp.h"

/* validate cgp program command line and create file strings for .par and data files */
void validate_command_line(int argc, char* argv[], char parfile[], char datafile[])
{
	puts("");
	puts("*********    WELCOME TO CARTESIAN GENETIC PROGRAMMING ANN      *********");
	puts("********* Validating command line arguments to cgp program *********");

	if (argc != 3)
	{
		puts("INCORRECT NUMBER OF ARGUMENTS");
		puts("Type cgp <file.par> <file> then return");
		puts("where file is a formatted file of test input/output data");
		exit(1);
	}
	if (strlen(argv[1])>(MAX_NUM_LETTERS - 1))
	{
		puts("filename for parameter file is too long");
		printf("It should be less than %d characters\n", MAX_NUM_LETTERS);
		exit(2);
	}
	strcpy(parfile, argv[1]);
	strcpy(datafile, argv[2]);

	if (strlen(argv[2])>(MAX_NUM_LETTERS - 1))
	{
		puts("filename for data file is too long");
		printf("It should be less than %d characters\n", MAX_NUM_LETTERS);
		exit(3);
	}

}

/* read from the parameter file all the global parameters */
void get_parameters(char parfile[MAX_NUM_LETTERS], char datafile[MAX_NUM_LETTERS])
{
	int		i;
	int     max_arity;
	char	dummy[50];
	FILE*	fp;

	printf("\n********* Reading parameters defined in %s *********\n", parfile);
	fp = fopen(parfile, "r");
	if (!fp)
	{
		printf("Missing file: %s\n", parfile);
		exit(1);
	}
	fscanf(fp, "%d %s", &population_size, dummy);
	fscanf(fp, "%lf %s", &per_cent_mutate, dummy);
	fscanf(fp, "%d %s", &num_generations, dummy);
	fscanf(fp, "%d %s", &num_runs_total, dummy);
	fscanf(fp, "%d %s", &num_rows, dummy);
	fscanf(fp, "%d %s", &num_cols, dummy);
	fscanf(fp, "%d %s", &levels_back, dummy);
	fscanf(fp, "%d %s", &progress_report, dummy);
	fscanf(fp, "%d %s", &report_interval, dummy);
	fscanf(fp, "%u %s", &global_seed, dummy);
	fscanf(fp, "%d %s", &save_best_chrom, dummy);
	fscanf(fp, "%d %s", &run_from_chrom, dummy);
	fscanf(fp, "%d %s", &shrink_phenotype, dummy);
	fscanf(fp, "%d %s", &arity, dummy);

	num_nodes = num_rows * num_cols;
	num_functions = 0;
	
	if (arity > MAX_NUM_INPUTS) max_arity = MAX_NUM_INPUTS; else max_arity = arity;
	for (i = 0; i < MAX_NUM_FUNCTIONS; i++)
	{
		/* number[] holds whether the function is used or not */
		fscanf(fp, "%d%s", &number[i], &node_types[i]);
		if (number[i])
		{
			allowed_functions[num_functions] = i;
			num_functions++;

			}
	}
	fclose(fp);

	/* each node is assigned max_arity connection genes */
	num_genes_per_node = (2 * max_arity) + 1;

	/* get input data */
	read_data(datafile);

	/* calculate the perfect score and
	for the boolean case what the bit width is

	*/
	define_perfect();


	if (population_size > MAX_NUM_CHROMOSOMES)
	{
		printf("Too large a population size (<= %d)\n", MAX_NUM_CHROMOSOMES);
		exit(0);
	}

	if (num_genes > MAX_NUM_GENES)
	{
		printf("Too many genes selected (<= %d)\n", MAX_NUM_GENES);
		exit(0);
	}

	if (num_runs_total < 1)
	{
		puts("Number of runs of EA must be at least 1");
		exit(0);
	}
	else if (num_runs_total > MAX_NUM_RUNS)
	{
		printf("Number of runs of EA must be less than %d\n", MAX_NUM_RUNS);
		exit(0);
	}

	if (num_genes < 10)
	{
		puts("Number of genes/bits must be at least 10");
		exit(0);
	}

	if ((progress_report< 0) || (progress_report > 1))
	{
		puts("Progress report parameter must be 0 or 1");
		exit(0);
	}

	srand(global_seed);

	/* assigned global constants */
	num_nodes = num_rows*num_cols;

	puts("********* Beginning execution *********");
}

/* write out parameter values in results file */
void write_cgp_info(char command[], char datafile[MAX_NUM_LETTERS])
{
	int		i;
	FILE*	fp;

	fp = fopen("cgp.txt", "w");
	fprintf(fp, "The program is        %s\n", command);
	fprintf(fp, "The data file is       %s\n", datafile);
	fprintf(fp, "population_size is    %d\n", population_size);
	fprintf(fp, "mutation rate is      %6.2lf\n", per_cent_mutate);
	fprintf(fp, "num_generations is    %d\n", num_generations);
	fprintf(fp, "num_runs is           %d\n", num_runs_total);
	fprintf(fp, "num_rows is           %d\n", num_rows);
	fprintf(fp, "num_cols is           %d\n", num_cols);
	fprintf(fp, "levels_back is        %d\n", levels_back);
	fprintf(fp, "progress report is    %d\n", progress_report);
	fprintf(fp, "report interval is    %d\n", report_interval);
	fprintf(fp, "global_seed is        %u\n", global_seed);
	fprintf(fp, "save_best_chrom is    %d\n", save_best_chrom);
	fprintf(fp, "run_from_chrom is     %d\n", run_from_chrom);
	fprintf(fp, "shrink_phenotype is     %d\n", shrink_phenotype);
	fprintf(fp, "perfect score is     %d\n", perfect);
	fprintf(fp, "arity is     %d\n", arity);

	for (i = 0; i<MAX_NUM_FUNCTIONS; i++)
	{
		fprintf(fp, "%d %s\n", number[i], node_types[i]);
	}
	fprintf(fp, "\nHere are the Results\n");
	fclose(fp);
}

/*  returns a random integer between 0 and range-1 */
int newrand(int range)
{
	int temp;

	temp = rand() % range;
	return(temp);
}

/* read an item of data from problem
specification file. Format of data
depends on defined data type
*/
data_type myfscanf(FILE* fp)
{
	data_type datum_read;


#ifdef DATA_IS_DOUBLE
	fscanf(fp, "%lf", &datum_read);
#endif

	return datum_read;
}


/* reads input and output data from a file  (e.g. compressed truth table .plu)
defines number of inputs to be num_inputs+num_constant_inputs
*/
void read_data(char datafile[MAX_NUM_LETTERS])
{
	int		i, j;
	char	dummy[MAX_NUM_LETTERS];
	FILE*	fp;

	fp = fopen(datafile, "r");
	if (!fp)
	{
		puts("ERROR. Missing input data file (e.g. .plu (compressed Boolean) .dat");
		exit(1);
	}
	else
	{
		fscanf(fp, "%s %d", dummy, &num_inputs);
		fscanf(fp, "%s %d", dummy, &num_outputs);
		fscanf(fp, "%s %d", dummy, &num_tests);
		if (num_tests >= MAX_NUM_DATA)
		{
			printf("\nERROR. Too many test cases (in datafile)\n");
			exit(0);
		}

		for (i = 0; i < num_tests; i++)
		{
			for (j = 0; j<num_inputs; j++)
				data_inputs[i][j] = myfscanf(fp);
			for (j = 0; j<num_outputs; j++)
				data_outputs[i][j] = myfscanf(fp);
		}
		fclose(fp);
	}
	num_genes = num_genes_per_node*num_nodes + num_outputs;


}

/* calculates a perfect score (if any) */
void define_perfect(void)
{
#ifdef DATA_IS_DOUBLE
	perfect = num_tests;
#endif
}


/* prints a chromosome to a file
when append is 1, the function appends the information to the file
when append is 0, the function creates a new file
*/
void fprint_a_chromosome(int* chromosome, char name[], int append)
{
	int		i, node_label;
	int      write_bracket = 1;
	FILE*	fp;

	if (append)
		fp = fopen(name, "a");
	else
		fp = fopen(name, "w");
	node_label = num_inputs - 1;
	for (i = 0; i<num_nodes*num_genes_per_node; i++)
	{
		if ((i + 1) % num_genes_per_node == 0)
		{
			node_label++;
			fprintf(fp, "[%d]:%d)\t", chromosome[i], node_label);
			write_bracket = 1;
		}
		else
		{
			if (write_bracket)
				fprintf(fp, "(");

			fprintf(fp, "%d,", chromosome[i]);
			write_bracket = 0;
		}
	}
	fprintf(fp, "\t\t");
	for (i = 0; i<num_outputs; i++)
		fprintf(fp, " %d", chromosome[num_nodes*num_genes_per_node + i]);
	fprintf(fp, "\n\n");
	fclose(fp);
}

/* prints a chromosome to the screen */
void print_a_chromosome(int* chromosome)
{
	int i;

	for (i = 0; i<num_nodes*num_genes_per_node; i++)
	{
		if ((i + 1) % num_genes_per_node == 0)
			printf(" %d\t", chromosome[i]);
		else
			printf(" %d", chromosome[i]);
	}
	printf("\t\t");
	for (i = 0; i<num_outputs; i++)
		printf(" %d", chromosome[num_nodes*num_genes_per_node + i]);
	printf("\n");
}

/* prints a chromosome to file in raw format,
so that it can be read by the program
*/
void fprint_a_raw_chromosome(int* chromosome, char name[], int append)
{
	int i;
	FILE* fp;

	if (append)
		fp = fopen(name, "a");
	else
		fp = fopen(name, "w");

	for (i = 0; i < num_nodes*num_genes_per_node; i++)
	{
		if ((i + 1) % num_genes_per_node == 0)
			fprintf(fp, " %d\t", chromosome[i]);
		else
			fprintf(fp, " %d", chromosome[i]);
	}
	fprintf(fp, "\t\t");

	for (i = 0; i<num_outputs; i++)
		fprintf(fp, " %d", chromosome[num_nodes*num_genes_per_node + i]);

	printf("\n");

	fclose(fp);
}

/* prints out array to a file */
void fprint_node_used(int size, int array[MAX_NUM_NODES_PLUS_INPUTS], char name[], int append)
{
	int i;
	FILE* fp;

	if (append)
		fp = fopen(name, "a");
	else
		fp = fopen(name, "w");

	fprintf(fp, "\nnode_used is now\n");
	fprintf(fp, "  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20\n");
	for (i = 0; i < size; i++)
		fprintf(fp, "%3d", array[i]);

	fclose(fp);
}


/* calculates the addresses of nodes used
and stores them in node_to_process
returns the number of nodes used
*/
int get_nodes_to_process(int* chromosome, int nodes_to_process[MAX_NUM_NODES])
{
	int		i, j, index;
	int		num_nodes_to_process;
	int		node_genes[MAX_NUM_GENES_PER_NODE];
	int     node_used[MAX_NUM_NODES_PLUS_INPUTS];
	int		max_size_node_used;

	max_size_node_used = num_nodes + num_inputs;

	/* say all nodes not used */
	for (i = 0; i < max_size_node_used; i++)
		node_used[i] = 0;

	/* all the nodes whose output is given by the output genes are active */
	/* last num_outputs genes are all output genes */
	for (i = num_genes - num_outputs; i < num_genes; i++)
		node_used[chromosome[i]] = 1;

	for (i = max_size_node_used - 1; i >= num_inputs; i--)
	{
		if (node_used[i])
		{
			/* get input addresses and type of this gate */
			index = num_genes_per_node*(i - num_inputs);

			/* write genes for node into array node_genes */
			for (j = 0; j < num_genes_per_node; j++)
				node_genes[j] = chromosome[index + j];

			/* each function has an arity stored in
			allowed_functions[][1].
			Find the nodes whose data is used
			*/


			for (j = 0; j < (2 * arity); j+=2)
				node_used[node_genes[j]] = 1;

		}
	}
	/* find number of used nodes */
	num_nodes_to_process = 0;
	for (i = num_inputs; i < max_size_node_used; i++)
	{
		if (node_used[i])
		{
			nodes_to_process[num_nodes_to_process] = i;
			num_nodes_to_process++;
		}
	}
	return num_nodes_to_process;
}



/**** print chromosome to file and indicate inactive genes with -1 */
void fprint_active_genes(int* chromosome, char name[30])
{
	int		i, j, index;
	int		write_bracket;
	int		node_label;
	int		num_unused_nodes, num_nodes_active;
	int		node_genes[MAX_NUM_GENES_PER_NODE];
	int     node_used[MAX_NUM_NODES_PLUS_INPUTS];
	int		max_size_node_used;

	int*	active_chromosome = NULL;
	FILE*	fp;


	active_chromosome = create_1d_space(num_genes);

	max_size_node_used = num_nodes + num_inputs;

	fp = fopen(name, "a");
	for (i = 0; i<num_genes; i++)
		active_chromosome[i] = -1;
	for (i = num_genes - num_outputs; i<num_genes; i++)
		active_chromosome[i] = chromosome[i];

	/* say all nodes not used */
	for (i = 0; i < max_size_node_used; i++)
		node_used[i] = 0;

	/* all the nodes whose output is given by the output genes are active */
	/* last num_outputs genes are all output genes */
	for (i = num_genes - num_outputs; i < num_genes; i++)
		node_used[chromosome[i]] = 1;

	for (i = max_size_node_used - 1; i >= num_inputs; i--)
	{
		if (node_used[i])
		{
			/* get input addresses and type of this gate */
			index = num_genes_per_node*(i - num_inputs);

			/* write genes for node into array node_genes */
			for (j = 0; j < num_genes_per_node; j++)
			{
				node_genes[j] = chromosome[index + j];
				active_chromosome[index + j] = node_genes[j];
			}

			/* each function has an arity stored in
			allowed_functions[][1].
			Find the nodes whose data is used
			*/
			for (j = 0; j < (2 * arity); j+=2)
				node_used[node_genes[j]] = 1;

		}
	}

	node_label = num_inputs - 1;
	write_bracket = 1;

	for (i = 0; i < num_inputs; i++)
	{
		if (node_used[i] > 0) fprintf(fp, "(%d)\t", i);
	}

	for (i = 0; i < num_nodes*num_genes_per_node; i++)
	{
		if ((i + 1) % num_genes_per_node == 0)
		{
			node_label++;
			if (active_chromosome[i]<0)
				fprintf(fp, "[*]:%d)\t", node_label);
			else
				fprintf(fp, "[%d]:%d)\t", active_chromosome[i], node_label);
			write_bracket = 1;
		}
		else
		{
			if (write_bracket == 1)
				fprintf(fp, "(");

			if (active_chromosome[i] == -1)
				fprintf(fp, "*,");
			else
				fprintf(fp, "%d,", active_chromosome[i]);

			write_bracket = 0;
		}
	}
	fprintf(fp, "\t\t");
	for (i = 0; i < num_outputs; i++)
		fprintf(fp, " %d", active_chromosome[num_nodes*num_genes_per_node + i]);

	num_unused_nodes = 0;
	for (i = num_inputs; i < num_inputs + num_nodes; i++)
	if (!node_used[i])
		num_unused_nodes++;

	num_nodes_active = num_nodes - num_unused_nodes;

	fprintf(fp, "\nnumber of active gates is %d\n\n", num_nodes_active);
	fclose(fp);

	free(active_chromosome);
}

/**** print chromosome to file and indicate inactive genes with -1 */
void fprint_active_genes2(int* chromosome, char name[30])
{
	int		i, j, index;
	int		write_bracket;
	int		node_label;
	int		num_unused_nodes, num_nodes_active;
	int		node_genes[MAX_NUM_GENES_PER_NODE];
	int     node_used[MAX_NUM_NODES_PLUS_INPUTS];
	int		max_size_node_used;

	int*	active_chromosome = NULL;
	FILE*	fp;


	active_chromosome = create_1d_space(num_genes);

	max_size_node_used = num_nodes + num_inputs;

	fp = fopen(name, "w");
	for (i = 0; i<num_genes; i++)
		active_chromosome[i] = -1;
	for (i = num_genes - num_outputs; i<num_genes; i++)
		active_chromosome[i] = chromosome[i];

	/* say all nodes not used */
	for (i = 0; i < max_size_node_used; i++)
		node_used[i] = 0;

	/* all the nodes whose output is given by the output genes are active */
	/* last num_outputs genes are all output genes */
	for (i = num_genes - num_outputs; i < num_genes; i++)
		node_used[chromosome[i]] = 1;

	for (i = max_size_node_used - 1; i >= num_inputs; i--)
	{
		if (node_used[i])
		{
			/* get input addresses and type of this gate */
			index = num_genes_per_node*(i - num_inputs);

			/* write genes for node into array node_genes */
			for (j = 0; j < num_genes_per_node; j++)
			{
				node_genes[j] = chromosome[index + j];
				active_chromosome[index + j] = node_genes[j];
			}

			/* each function has an arity stored in
			allowed_functions[][1].
			Find the nodes whose data is used
			*/
			for (j = 0; j < arity; j++)
				node_used[node_genes[j]] = 1;

		}
	}

	node_label = num_inputs - 1;

	num_unused_nodes = 0;
	for (i = num_inputs; i < num_inputs + num_nodes; i++)
	if (!node_used[i])
		num_unused_nodes++;

	num_nodes_active = num_nodes - num_unused_nodes;

	fprintf(fp, "%d\n", num_inputs);
	fprintf(fp, "%d\n", num_outputs);
	fprintf(fp, "%d\n", num_nodes);
	fprintf(fp, "%d\n", num_genes_per_node);
	fprintf(fp, "%d\n", num_nodes_active);

	for (i = 0; i < num_nodes*num_genes_per_node; i++)
	{
		if ((i + 1) % num_genes_per_node == 0)
		{
			node_label++;
			if (active_chromosome[i]> -1)
			{
				fprintf(fp, "%d %d\t", active_chromosome[i], node_label);
			}

		}
		else
		if (active_chromosome[i]> -1)
			fprintf(fp, "%d ", active_chromosome[i]);
	}

	fprintf(fp, "\n");
	for (i = 0; i < num_outputs; i++)
		fprintf(fp, "%d ", active_chromosome[num_nodes*num_genes_per_node + i]);

	fclose(fp);

	free(active_chromosome);
}

/* generate a starting population
from a chromosome read from a file (cgp.chr)
*/
void read_from_chrom(int** chromosomes)
{
	int	  i, j;
	FILE* fp;

	fp = fopen("cgp.chr", "r");
	if (!fp)
	{
		puts("Missing file cgp.chr (contains a chromosome)");
		exit(1);
	}
	else
	{
		/* make starting population copies of loaded chromosome */
		for (j = 0; j<population_size; j++)
		{
			if (j == 0)
			{
				i = 0;
				do
				{
					fscanf(fp, "%d", &chromosomes[j][i]);
					i++;
				} while (!feof(fp));

				if (i != num_genes)
				{
					puts("ERROR. Number of genes in cgp.chr does not match the expected number");
					printf("\nnum_genes required is %d, num_genes read is %d", num_genes, i);
					puts("Check the number of genes in the .par file");
					exit(0);
				}
			}
			else
			{
				for (i = 0; i<num_genes; i++)
					chromosomes[j][i] = chromosomes[0][i];
			}
		}
		fclose(fp);
	}
}


/* this decodes the cgp chromosome.
It is given data_inputs corresponding to a single
test case and it calculates what the cgp genotype gives
for the data outputs (cgp_outputs)
It only processes nodes that are used (nodes_to_process)
*/
void decode_cgp(int* chromosome,
	data_type data_inputs[MAX_NUM_DATA][MAX_NUM_INPUTS],
	data_type cgp_outputs[MAX_NUM_OUTPUTS],
	int num_nodes_to_process,
	int nodes_to_process[MAX_NUM_NODES],
	int fitness_test)
{
	int			i, j;
	int			first_node_gene, function_type;
	int			node_index;
	data_type	in[MAX_NUM_GENES_PER_NODE];
	data_type	output[MAX_OUTPUT_SIZE];
	data_type cgp_data;

	/* load test data into output array */
	for (j = 0; j < num_inputs; j++)
		output[j] = data_inputs[fitness_test][j ];

	/* only process nodes that are used */
	for (j = 0; j< num_nodes_to_process; j++)
	{
		/* get address of node */
		node_index = nodes_to_process[j] - num_inputs;

		/* get address of first used gene in node */
		first_node_gene = num_genes_per_node*node_index;

		for (i = 0; i < num_genes_per_node - 1; i++)						/* get input data to node */
		{
			if (i % 2 == 0) in[i] = output[chromosome[first_node_gene + i]];
			else in[i] = chromosome[first_node_gene + i];
		}
		function_type = chromosome[first_node_gene + num_genes_per_node - 1];	/* get node function */

		output[node_index + num_inputs] = node_type(in, function_type, arity);				    /* compute output of node and store */

	}

	/* process outputs */
	for (j = 0; j < num_outputs; j++)
	{
		//cgp_outputs[j] = output[chromosome[num_genes - num_outputs + j]];    /* get output from nodes referenced in output genes */
		// if threshold is used
		cgp_data = output[chromosome[num_genes - num_outputs + j]];    /* get output from nodes referenced in output genes */
		cgpout2[fitness_test][j] = cgp_data;
		/*for (i = 1; i < (NUM_CLASSES) ; i++)
		{
		if (cgp_data < (double)i / NUM_CLASSES) cgp_outputs[j] = (double)(i - 1) / (NUM_CLASSES);
		}*/
		if (cgp_data < 1.0 / NUM_CLASSES) cgp_outputs[j] = 0;
		else if (cgp_data < 2.0 / NUM_CLASSES) cgp_outputs[j] = 0.5;
		else cgp_outputs[j] = 1.0;
	}
}



#ifdef DATA_IS_DOUBLE

/* checks how close the output double cgp produces
is to the desired output integer
returns a 1 each time the absolute difference between
the desired output and the cgp output is less than a
user defined error. This is a hits based fitness measure
suitable for symbolic regression problems
*/
double correctness_test(data_type data_output,
	data_type cgp_output)
{
	double	result = 0.0;
	double x;

	x = fabs(cgp_output - data_output);

	/* hits based fitness */
#ifdef HITS_BASED_FITNESS
	if (x < ERROR_THRESHOLD)
		result = 1;
#else
	result = 1.0 / (1.0 + x); /* error based fitness */
#endif

	return result;
}

#endif


/* evaluate the fitness at a single test point */
double evaluate_cgp_outputs(data_type cgp_outputs[MAX_NUM_OUTPUTS],
	int test)
{
	int i;
	double fit = 0.0, x;
	
	for (i = 0; i < num_outputs; i++)
	{
		x = correctness_test(data_outputs[test][i], cgp_outputs[i]);
		fit = fit + x;
	}

	return fit;
}



/* this is the EA fitness function
*/
double fitness(int* chromosome)
{
	int			fitness_test;
	int			num_nodes_to_process;
	int			nodes_to_process[MAX_NUM_NODES];
	double		fit = 0.0;
	data_type	cgp_outputs[MAX_NUM_OUTPUTS];

	/* find out how many nodes there are in the phenotype */
	num_nodes_to_process = get_nodes_to_process(chromosome, nodes_to_process);

	/* apply all fitness tests */
	for (fitness_test = 0; fitness_test < num_tests; fitness_test++)
	{
		decode_cgp(chromosome, data_inputs, cgp_outputs,
			num_nodes_to_process, nodes_to_process, fitness_test);

		fit = fit + evaluate_cgp_outputs(cgp_outputs, fitness_test);

	}

	/* if a perfect solution is found and shrink_phenotype is set to 1
	this adds how many unused nodes there are to the fitness score.
	this results in phenotypes shrinking
	*/
	if ((fit == perfect) && (shrink_phenotype))
		fit = fit + num_nodes - num_nodes_to_process;

	return fit;
}

/* This calculates the limits that are used in the claculation
of allowed gene values (alleles) */
void get_gene_limits(int column, int* limit_min, int* limit)
{
	int limit_max;

	limit_max = num_inputs + column*num_rows;
	if (column<levels_back)
		*limit_min = 0;
	else
		*limit_min = num_inputs + (column - levels_back)*num_rows;

	*limit = limit_max - (*limit_min);
}

/* returns a random valid connection gene that
obeys the constraints imposed by levels_back.
Also allows program inputs to disobey levels_back */
int get_connection_gene(int limit_min, int limit)
{
	int limit_plus, rand_num;
	int gene;

	if (limit_min == 0)
		gene = newrand(limit);
	else /* allows inputs to disobey levels_back */
	{
		limit_plus = limit + num_inputs;
		rand_num = newrand(limit_plus);
		if (rand_num<limit)
			gene = rand_num + limit_min;
		else
			gene = rand_num - limit;
	}

	return gene;
}

/* returns a random valid function gene */
int get_function_gene(void)
{
	return allowed_functions[newrand(num_functions)];
}

/* returns a random valid output gene */
int get_output_gene(void)
{
	int limit_min, limit;
	int output_gene;

	limit_min = num_inputs + (num_cols - levels_back)*num_rows;
	limit = levels_back*num_rows;

	output_gene = newrand(limit) + limit_min;

	return output_gene;
}


/* checks to see if the gene is not an output gene */
int is_not_output_gene(int gene)
{
	return (gene < num_genes_per_node*num_nodes);
}

int get_weight(void)
{
	int num;
	num = rand() - (RAND_MAX / 2);
	return num;
}
/* checks to see if the gene is a function gene */
int is_function_gene(int gene, int locus)
{
	return (is_not_output_gene(gene) && (locus == (num_genes_per_node - 1)));

}

int is_weight_gene(int gene, int locus)
{
	return (is_not_output_gene(gene) && (locus % 2 == 1));

}

/* generates a random chromosome. Used by initialise */
void generate_a_random_chromosome(int* chromosome)
{
	int i, count = 0;
	int row, col, limit, limit_min;

	for (col = 0; col < num_cols; col++)
	{
		get_gene_limits(col, &limit_min, &limit);

		for (row = 0; row<num_rows; row++)
		{
			/* get random input genes */
			for (i = 0; i < num_genes_per_node - 1; i++)
			{
				switch (i % 2)
				{
				case 0:
					chromosome[count + i] = get_connection_gene(limit_min, limit);// here connection refers to connection to input
					break;
				case 1:
					chromosome[count + i] = get_weight();// connection gene can be 0 or 1
					break;
				default:
					break;
				}

			}

			/* get random function gene */
			chromosome[count + num_genes_per_node - 1] = get_function_gene();
			count = count + num_genes_per_node;
		}
	}
	/* get random function genes */
	for (i = 0; i < num_outputs; i++)
		chromosome[count + i] = get_output_gene();
}


/* creates initial population of chromosomes
either having been generated from a single
chromosome from a file or by generating
an entire random population
*/
void initialise(int** chromosomes)
{
	int  i;

	if (run_from_chrom)
		read_from_chrom(chromosomes);
	else  /* generate random population */
	for (i = 0; i < population_size; i++)
		generate_a_random_chromosome(chromosomes[i]);
}



/* calculate best population fitness and the best chromosome */
double  get_best_chromosome(int** chromosomes,
	int*  best_chromosome,
	double previous_best_fitness,
	int    gen)
{
	int		i;
	double	fitness_max, fit;
	int		best_member;

	fitness_max = -1.0;
	best_member = 0;

	for (i = 0; i < population_size; i++)
	{

		if ((i == population_size - 1) && (gen > 1))
			fit = previous_best_fitness;
		else
			fit = fitness(chromosomes[i]);

		if (fit > fitness_max)
		{
			fitness_max = fit;
			best_member = i;
		}

		/* break out of this as soon as we get a perfect score
		and shrink_phenotype is not required */
		if ((fit == perfect) && (shrink_phenotype == 0))
			break;
	}

	/* store the best chromosome */
	for (i = 0; i<num_genes; i++)
		best_chromosome[i] = chromosomes[best_member][i];

	return fitness_max;
}



/* calculates how many mutations to do per chromosome */
int get_num_mutant(int num_genes, double per_cent_mutate)
{
	return (int)(num_genes*per_cent_mutate / 100.0);
}

/* mutate one gene in a chromosome */
void mutate_a_gene(int*  chromosome)
{
	int which_gene, which_locus;
	int limit, limit_min;
	int col;

	which_gene = newrand(num_genes);
	which_locus = which_gene % num_genes_per_node;

	if (is_not_output_gene(which_gene))
	{
		if (is_function_gene(which_gene, which_locus))
		{
			if (num_functions == 1) /* redirect the mutation to a connection */
			{
				which_locus = newrand(num_genes_per_node - 1);
				which_gene = which_gene - num_genes_per_node - 1 + which_locus;
			}
			chromosome[which_gene] = get_function_gene();
		}

		else if (is_weight_gene(which_gene, which_locus)) 	chromosome[which_gene] = get_weight();

		else /* it is a connection gene */
		{
			col = which_gene / (num_genes_per_node*num_rows);

			get_gene_limits(col, &limit_min, &limit);

			chromosome[which_gene] = get_connection_gene(limit_min, limit);
		}
	}
	else /* it is an output gene */
	{
		chromosome[which_gene] = get_output_gene();
	}
}

/* carry out num_mutations mutations on the chromosome */
void mutate_a_chromosome(int *chromosome, int num_mutations)
{
	int i;

	for (i = 0; i < num_mutations; i++)
		mutate_a_gene(chromosome);
}

/* (1+lambda evolutionary strategy where lamda = population size -1 */
void generate_new_pop_es(int** chromosomes,
	int*  best_chromosome)
{
	int j, k;
	int num_mutant;

	num_mutant = get_num_mutant(num_genes, per_cent_mutate);

	/* copy best_chromosome into last member of chromosome array */
	for (j = 0; j<num_genes; j++)
		chromosomes[population_size - 1][j] = best_chromosome[j];

	/* generate new population by mutating all but last */
	for (k = 0; k<population_size - 1; k++)
	{
		for (j = 0; j<num_genes; j++) /* copy best chromosome */
			chromosomes[k][j] = best_chromosome[j];

		/* mutate the chromosome */
		mutate_a_chromosome(chromosomes[k], num_mutant);
	}
}

/* allocate space for 2 dimensional array
e.g. a population of chromosomes */
int** create_2d_space(int num_horizontals, int num_verticals)
{
	int i;
	int **array2d = NULL;

	/* create space for pointers to int pointers */
	array2d = (int**)calloc(num_verticals, sizeof(int*));
	if (array2d == NULL)
	{
		printf("ERROR.Can not allocate space for %d many int pointers\n", num_verticals);
		exit(0);
	}

	/* create array of pointers to ints  */
	for (i = 0; i < num_verticals; i++)
	{
		array2d[i] = create_1d_space(num_horizontals);
		if (array2d[i] == NULL)
		{
			printf("ERROR.Not enough memory for int pointer arrays of length %d\n", num_horizontals);
			exit(0);
		}
	}

	return array2d;
}

/* allocate space for a 1d array of with size items */
int* create_1d_space(int size)
{
	int* array1d = NULL;

	array1d = (int*)calloc(size, sizeof(int));

	if (array1d == NULL)
	{
		printf("ERROR.Not enough memory for a 1d array of length %d\n", size);
		exit(0);
	}
	return array1d;
}


/* release memory */
/* if this is for chromosomes then num_verticals is population_size */
void free_array2d(int num_verticals, int** array2d)
{
	int i;

	/* free 1darray of pointers  */
	for (i = 0; i < num_verticals; i++)
		free(array2d[i]);

	free(array2d);
}

void write_generation_to_screen(int gen_index)
{
	if (gen_index % report_interval == 0)
		printf("\nGENERATION is %d", gen_index);
}

/* writes generation and fitness of best in population */
void write_progress_info_to_screen(int generation, double fit)
{
	printf("\nGENERATION is %d Best fitness is now %8.5lf that is %8.5lf%% percent.", generation, fit, fit / perfect * 100);
}


/* writes out chromosome to file defined by string prog */
void write_progress_info_to_file(char prog[MAX_NUM_LETTERS],
	int gen, double best_fit,
	int* best_chromosome)
{
	FILE* fp;

	fp = fopen(prog, "a");
	fprintf(fp, "\nGENERATION is %u     Best fitness is now %8.5lf", gen, best_fit);
	fprintf(fp, "\nThe chromosome is\n");
	fclose(fp);
	fprint_a_chromosome(best_chromosome, prog, 1);
	fprint_active_genes(best_chromosome, prog);
}

/* checks to see if the best fitness in the population has improved.
writes the generation, the new best fitness and the improved chromosome
to the progress report (if progress report is set to 1)
*/
void check_if_improvement(double best_fit, double* previous_best_fit, int* best_gen, int gen,
	char prog[MAX_NUM_LETTERS],
	int* best_chromosome)
{
	if (best_fit > *previous_best_fit) /* we have an improvement */
	{
		if (progress_report)
		{
			write_progress_info_to_screen(gen, best_fit);
			write_progress_info_to_file(prog, gen, best_fit, best_chromosome);
		}
		*best_gen = gen;				/* update the best generation */
		*previous_best_fit = best_fit;	/* update previous best fitness */
	}
}


/* report on results of evolutionary run in cgp.txt */
void write_result_of_EA_to_file(int run, int bgen, double best_fit,
	int* best_chromosome)
{
	FILE* fp;

	fp = fopen("cgp.txt", "a");
	fprintf(fp, "Run %d and gen %d achieved fitness %6.2lf\n", run, bgen, best_fit);
	fprintf(fp, "Here is the chromosome\n");
	fprintf(fp, "Perfect = %d \n", perfect);
	fclose(fp);

	fprint_a_chromosome(best_chromosome, "cgp.txt", 1);
	fprint_active_genes(best_chromosome, "cgp.txt");
}

/* Do a run of the EA */
double  EA(int *gen_of_best, int* num_nodes_active, int run, double best_of_best_fit,
	char prog[MAX_NUM_LETTERS])
{
	int		gen, best_gen;
	int		nodes_to_process[MAX_NUM_NODES];
	int**	chromosomes;
	int*	best_chromosome;
	double	best_fit, best, previous_best_fit = -1.0;
	int flag = 0;
	//best_fit = 0;
	chromosomes = create_2d_space(num_genes, population_size);
	best_chromosome = create_1d_space(num_genes);
	initialise(chromosomes);
	if (per_cent_mutate == 0.0) flag = 1;

	for (gen = 1; gen <= num_generations; gen++)
	{
		write_generation_to_screen(gen);
		if (flag == 1)
		{
			/*factor = (double) gen/num_generations;
			per_cent_mutate = max_mutate * (1.0 - factor);*/
			//best = best_fit / num_tests;
			if (gen < num_generations/4) per_cent_mutate = 20.0;
			else if (gen < num_generations / 2) per_cent_mutate = 10.0;
			else if (gen < 3 * num_generations / 4) per_cent_mutate = 5.0;
			else per_cent_mutate = 1.0;
		}
		best_fit = get_best_chromosome(chromosomes,
			best_chromosome,
			previous_best_fit, gen);

		check_if_improvement(best_fit, &previous_best_fit, &best_gen, gen, prog,
			best_chromosome);

		/* jump out of run if maximum fitness acheived */
		if ((best_fit == perfect) && (shrink_phenotype == 0))
			break;
		else /* create a new population */
			generate_new_pop_es(chromosomes, best_chromosome);
	}

	*num_nodes_active = get_nodes_to_process(best_chromosome, nodes_to_process);

	write_result_of_EA_to_file(run, best_gen, best_fit, best_chromosome);

	*gen_of_best = best_gen;

	/* write the raw best chromosome to cgp.chr */
	if (best_fit> best_of_best_fit)// Added by Gulrukh. best_of_best_fit added to parameters.
	{
		fprint_a_raw_chromosome(best_chromosome, "cgp.chr", 0);
		//fprint_active_genes2(best_chromosome, "best.chr");
		//get_cgp_outputs(best_chromosome);

	}

//	free_array2d(population_size, chromosomes);
	free(best_chromosome);
	return best_fit;
}

void setup_report_files(int run, char prog[MAX_NUM_LETTERS])
{
	char runstring[MAX_NUM_LETTERS];
	FILE* fp;

	sprintf(runstring, "%d", run); /* store run as characters */
	if (progress_report > 0)
	{
		strcpy(prog, "cgp");
		strcat(prog, runstring);
		strcat(prog, ".prg"); /* create .prg file name */
		fp = fopen(prog, "w");  /* create empty .prg file */
		fclose(fp);
	}
}

void report_final_results(double av_fitness, double st_dev,
	double av_num_nodes_active,
	double best_of_best_fit,
	double worst_of_best_fit,
	double av_best_gen,
	int    num_perfect)
{
	FILE* best;
	best = fopen("cgp.txt", "a");
	fprintf(best, "\naverage fitness  %6.4lf\n", av_fitness);
	fprintf(best, "\nstd dev          %6.4lf\n\n", st_dev);
	fprintf(best, "\naverage number of active nodes is  %6.4lf\n", av_num_nodes_active);
	fprintf(best, "\nThe best fitnes of all runs  is  %6.4lf\n", best_of_best_fit);
	fprintf(best, "\nThe worst fitness of all runs is  %6.4lf\n", worst_of_best_fit);
	if ((num_perfect > 0) && (shrink_phenotype == 0))
	{
		fprintf(best, "\nNumber of perfect solutions is %d\n", num_perfect);
		fprintf(best, "\nOf perfect solutions found, the average number of generations is  %6.4lf\n", av_best_gen);
	}
	fclose(best);
}

double average(int num_items, double items[MAX_NUM_RUNS])
{
	int		i;
	double	av;

	av = 0.0;
	for (i = 0; i < num_items; i++)
		av = av + items[i];

	av = av / num_items;

	return av;
}

double get_standard_deviation(int num_items, double average,
	double items[MAX_NUM_RUNS])
{
	int i;
	double temp, st_dev;

	st_dev = 0.0;
	for (i = 0; i < num_items; i++)
	{
		temp = (items[i] - average);
		temp = temp*temp;
		st_dev = st_dev + temp;
	}

	st_dev = st_dev / (double)num_items;
	st_dev = sqrt(st_dev);

	return st_dev;
}

/* do mutiple runs of EA and write out results */
void run_EA(int num_runs_total)
{
	int		best_gen, run, num_perfect = 0;
	int		num_nodes_active;
	double  active_nodes_in_run[MAX_NUM_RUNS];
	double	generations_for_best[MAX_NUM_RUNS];
	char	prog[MAX_NUM_LETTERS];
	double  worst_of_best_fit, best_of_best_fit;
	double  fitness_final;
	double	st_dev, av_fitness, av_best_gen;
	double  av_num_nodes_active;
	double	fitnesses[MAX_NUM_RUNS];

	best_of_best_fit = 0.0;
	worst_of_best_fit = perfect + num_nodes;
	for (run = 0; run < num_runs_total; run++)
	{

		setup_report_files(run, prog);

		fitness_final = EA(&best_gen, &num_nodes_active, run, best_of_best_fit, prog);

		fitnesses[run] = fitness_final;
		active_nodes_in_run[run] = num_nodes_active;
		generations_for_best[run] = best_gen;

		if (fitness_final < worst_of_best_fit)
			worst_of_best_fit = fitness_final;

		if (fitness_final > best_of_best_fit)
			best_of_best_fit = fitness_final;

		if ((fitness_final == perfect) && (shrink_phenotype == 0))
			num_perfect++;
	}
	av_best_gen = average(num_runs_total, generations_for_best);
	av_num_nodes_active = average(num_runs_total, active_nodes_in_run);
	av_fitness = average(num_runs_total, fitnesses);
	st_dev = get_standard_deviation(num_runs_total, av_fitness, fitnesses);
	report_final_results(av_fitness, st_dev, av_num_nodes_active,
		best_of_best_fit, worst_of_best_fit,
		av_best_gen, num_perfect);
}

void run_test(char testfile[MAX_NUM_LETTERS], char outfile[MAX_NUM_LETTERS])
{
	int* chrom;
	FILE* fp;
	int i;
	read_data(testfile);
	chrom = create_1d_space(num_genes);

	fp = fopen("cgp.chr", "r");
	if (!fp)
	{
		puts("Missing file cgp.chr (contains a chromosome)");
		exit(1);
	}
	else
	{
		i = 0;
		do
		{
			fscanf(fp, "%d", &chrom[i]);
			i++;
		} while (!feof(fp));

		if (i != num_genes)
		{
			puts("ERROR. Number of genes in cgp.chr does not match the expected number");
			printf("\nnum_genes required is %d, num_genes read is %d", num_genes, i);
			puts("Check the number of genes in the .par file");
			exit(0);
		}

	}
	fclose(fp);
	get_cgp_outputs(chrom, outfile);

}

/* Gulrukh
This function gets the output of the best chromosome into the global variable array cgpout2. It then prints the cgp out put as well as corressponding data output.
*/
void get_cgp_outputs(int* best_chromosome, char file[MAX_NUM_LETTERS])
{
	FILE* fo;
	int x, y, i;
	int positive[MAX_NUM_OUTPUTS], negative[MAX_NUM_OUTPUTS], ones[MAX_NUM_OUTPUTS], zeroes[MAX_NUM_OUTPUTS];
	int			fitness_test;
	int			num_nodes_to_process;
	int			nodes_to_process[MAX_NUM_NODES];
	data_type	cgp_outputs[MAX_NUM_OUTPUTS], final_outputs[MAX_NUM_DATA][MAX_NUM_OUTPUTS];
	double fit[MAX_NUM_OUTPUTS], error[MAX_NUM_OUTPUTS], sq_error[MAX_NUM_OUTPUTS], tot_error[MAX_NUM_OUTPUTS];


	/* find out how many nodes there are in the phenotype */
	num_nodes_to_process = get_nodes_to_process(best_chromosome, nodes_to_process);

	fo = fopen(file, "w");

	/* apply all fitness tests */
	for (fitness_test = 0; fitness_test < num_tests; fitness_test++)
	{
		decode_cgp(best_chromosome, data_inputs, cgp_outputs,
			num_nodes_to_process, nodes_to_process, fitness_test);

		for (i = 0; i < num_outputs; i++)
		{
			fit[i] = error[i] = sq_error[i] = tot_error[i] = 0.0;
			positive[i] = negative[i] = ones[i] = zeroes[i] = 0;
			final_outputs[fitness_test][i] = cgp_outputs[i];
		}
	}
	//prints output to file
	/*for (y = 0; y<num_tests; y++)
	{
		for (x = 0; x < num_outputs; x++)
		{
			fprintf(fo, "%lf	 %lf\t", cgpout2[y][x], final_outputs[y][x]);
			if (data_outputs[y][x] != final_outputs[y][x])
			{
				error[x] = error[x] + 1;
				if (data_outputs[y][x] > final_outputs[y][x]) positive[x] = positive[x] + 1; else negative[x] = negative[x] + 1;
			}
			if (data_outputs[y][x] == 1) ones[x] = ones[x] + 1; else zeroes[x] = zeroes[x] + 1;
			fit[x] = fit[x] + correctness_test(data_outputs[y][x], final_outputs[y][x]);
			sq_error[x] = sq_error[x] + ((data_outputs[y][x] - final_outputs[y][x]) * (data_outputs[y][x] - final_outputs[y][x])) / 2;
			tot_error[x] = tot_error[x] + (data_outputs[y][x] - final_outputs[y][x]);
		}
		fprintf(fo, "\t\t\t");
		for (x = 0; x<num_outputs; x++)fprintf(fo, " %lf\t\t", data_outputs[y][x]);
		fprintf(fo, "\n");
	}
	for (x = 0; x<num_outputs; x++)fprintf(fo, "Nodes used = %d	fitness = %lf	correctness= %lf	errors = %lf\n	positive = %d	negative = %d	ones =%d	zeroes = %d\t\t", num_nodes_to_process, fit[x], 100 * (num_tests - error[x])/num_tests, error[x], positive[x], negative[x], ones[x], zeroes[x]);

	fclose(fo);*/
	for (y = 0; y<num_tests; y++)
	{
		for (x = 0; x < num_outputs; x++)
		{
			fprintf(fo, "%lf	 %lf\t", cgpout2[y][x], final_outputs[y][x]);
			if (data_outputs[y][x] != final_outputs[y][x])
			{
				error[x] = error[x] + 1;
				if ((data_outputs[y][x] - final_outputs[y][x])> 0) positive[x] += 1; else negative[x] += 1;
			}
			fit[x] = fit[x] + correctness_test(data_outputs[y][x], final_outputs[y][x]);
			sq_error[x] = sq_error[x] + ((data_outputs[y][x] - final_outputs[y][x]) * (data_outputs[y][x] - final_outputs[y][x])) / 2;
			tot_error[x] = tot_error[x] + abs(data_outputs[y][x] - cgpout2[y][x]);
		}
		fprintf(fo, "\t\t\t");
		for (x = 0; x < num_outputs; x++)
		{
			fprintf(fo, " %lf\t\t", data_outputs[y][x]);
			if (data_outputs[y][x] == 1) ones[x] = ones[x] + 1; else zeroes[x] = zeroes[x] + 1;
		}
		fprintf(fo, "\n");
	}
	for (x = 0; x<num_outputs; x++)fprintf(fo, "Nodes used = %d	fitness = %lf(%lf)	\n errors = %lf	sum of squared error = %lf	total error = %lf\nfitness percentage = %lf\n positive = %d	negative =%d	ones=%d	zeroes =%d", num_nodes_to_process, fit[x], (num_tests - error[x]), error[x], sq_error[x], tot_error[x], 100 * (num_tests - error[x]) / num_tests, positive[x], negative[x], ones[x], zeroes[x]);

	fclose(fo);
}
