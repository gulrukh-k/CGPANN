/* cgp.c
Julian F. Miller (c), 2009
version 1.1 of first public release on 20-July-2009
Dept. of Electronics, University of York, UK
*/

#include <stdio.h>
#include "cgp.h" 


int main(int argc, char* argv[])
{
	char parfile[MAX_NUM_LETTERS], plufile[MAX_NUM_LETTERS];

	char testfile[MAX_NUM_LETTERS] = "iristrain1.dat";

	char trainfile[MAX_NUM_LETTERS] = "iristest1.dat";

	char testout[MAX_NUM_LETTERS] = "iristestout.txt";

	char trainout[MAX_NUM_LETTERS] = "iristrainout.txt";
	
	validate_command_line(argc, argv, parfile, plufile);

	get_parameters(parfile, plufile);

	write_cgp_info(argv[0], plufile);

	run_EA(num_runs_total);

	puts("");
	puts("*********    CGP COMPLETED      *********");

	run_test(trainfile, trainout);

	run_test(testfile, testout);
	return 0;
}

