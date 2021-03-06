/* node-function.c
Julian F. Miller (c) 2009
version 1.1 of first public release on 20-July-2009
Dept. of Electronics, University of York, UK

node_type calculates the output of a node given
the data provided in the array in
*/

#include "cgp.h"
#include <math.h>
#include <stdlib.h>



#ifdef DATA_IS_DOUBLE

data_type  node_type(data_type in[MAX_NUM_GENES_PER_NODE],
	int function_gene, int arity)
{

	data_type sum, result, weight;
	int i;
	result = 0; sum = 0;
	for (i = 0; i < arity; i++)
	{
		weight = in[((i * 2) + 1)] / (RAND_MAX/ 2);
		sum = sum + in[i * 2] * weight;
	}

	switch (function_gene)
	{
	case 0:  /*   sigmoid  */
		result = 1 / (1 + exp(-1 * sum));
		break;
	case 1:  /*   tanh  */
		result = tanh(sum);
		break;
	}
	return result;
}

#endif
