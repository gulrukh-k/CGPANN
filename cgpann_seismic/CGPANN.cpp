/* cgp.c
Julian F. Miller (c), 2009
version 1.1 of first public release on 20-July-2009
Dept. of Electronics, University of York, UK
*/

#include <stdio.h>
#include "CGPANN.h"


int main(int argc, char* argv[])
{
	char parfile[MAX_NUM_LETTERS], plufile[MAX_NUM_LETTERS];

	//char trainfile[MAX_NUM_LETTERS] = "noisetrain.txt";
	char trainfile[MAX_NUM_LETTERS] = "trainallavg.txt";
	
	//char testfile[MAX_NUM_LETTERS] = "noisetest.txt";                                                                                     
	char testfile[MAX_NUM_LETTERS] = "testallavg.txt";

	char testout[MAX_NUM_LETTERS] = "avgtestout2.txt";

	char trainout[MAX_NUM_LETTERS] = "avgtrainout2.txt";
	char datafile[MAX_NUM_LETTERS] = "noisedata.txt";
	char outfile[MAX_NUM_LETTERS] = "out2.txt";

	validate_command_line(argc, argv, parfile, plufile);

	get_parameters(parfile, plufile);

	write_cgp_info(argv[0], plufile);

	run_EA(num_runs_total);

	puts("");
	puts("*********    CGP COMPLETED      *********");
	run_test(trainfile, trainout);

	run_test(testfile, testout);

	//run_volcanotest(datafile, outfile);

	return 0;
}

