/* Code skeleton of a Scilab gateway function.
 */

#include "api_scilab.h"
#include "Scierror.h"
#include "malloc.h"
#include "sciprint.h"
#include <localization.h>

/* Place extra `#include's here.
 */
#include "skeleton.h"

/* Replace `skeleton' with your C-function name.
 */
int sci_skeleton(char *fname, unsigned long fname_len)
{
    SciErr      err;
    int         num_of_in_args;
    int         num_of_out_args;
    int         *pvaraddrs = NULL;
    int         *pvartypes = NULL;
    int         i;
    
    double      *in_var1;
    int         var1_rows, var1_cols;
    double      out_var1;
    
    /* Constants used by step 1. */
    const int   MIN_NUM_OF_IN_ARGS      = 1;
    const int   MAX_NUM_OF_IN_ARGS      = 1;
    const int   MIN_NUM_OF_OUT_ARGS     = 1;
    const int   MAX_NUM_OF_OUT_ARGS     = 1;
    
    /* Step 1.  Check the number of input / output arguments presented in the calling Scilab
     *          function.
     */
    CheckInputArgument(pvApiCtx, MIN_NUM_OF_IN_ARGS, MAX_NUM_OF_IN_ARGS);
    CheckOutputArgument(pvApiCtx, MIN_NUM_OF_OUT_ARGS, MAX_NUM_OF_OUT_ARGS);
    num_of_in_args = nbInputArgument(pvApiCtx);
    num_of_out_args = nbOutputArgument(pvApiCtx);

    /* Step 2.  Get addresses of input arguments and place them into `pvaraddrs[]'.
     */
    pvaraddrs = (int*)malloc(sizeof(int) * num_of_in_args);
    if (pvaraddrs == NULL)
    {
        Scierror(999, "Failed to allocate memory. \n\n");
        return 0;
    }

    for (i = 0; i < num_of_in_args; ++i)
    {
        err = getVarAddressFromPosition(pvApiCtx, i + 1, &pvaraddrs[i]);
        if (err.iErr)
        {
            printError(&err, 0);
            return 0;
        }
    }

    /* Step 3.  Check the data type of each variable in `pvaraddrs[]'.
     */
    pvartypes = (int*)malloc(sizeof(int) * num_of_in_args);
    if (pvartypes == NULL)
    {
        Scierror(999, "Failed to allocate memory. \n\n");
        return 0;
    }

    for (i = 0; i < num_of_in_args; ++i)
    {
        err = getVarType(pvApiCtx, &pvaraddrs[i], &pvartypes[i]);
        if (err.iErr)
        {
            printError(&err, 0);
            return 0;
        }
    }

    /* Step 4.  Check if the relevant input argument is complex or not, if needed. While
     *          dealing with integer, further checks should be done on the precision of
     *          the integer. Check the size of matrices and vectors. Check consistency of
     *          some input variables.
     */
    if (isVarComplex(pvApiCtx, pvaraddrs[0]))
    {
        Scierror(999, "Input argument 1 must be a real number. \n\n");
        return 0;
    }
    
    err = getMatrixOfDouble(pvApiCtx, pvaraddrs[0], &var1_rows, &var1_cols, &in_var1);
    if (err.iErr)
    {
        printError(&err, 0);
        return 0;
    }
    if ((var1_rows != 1) || (var1_cols != 1))
    {
        Scierror(999, "Wrong size of input argument 1. \n\n");
        return 0;
    }

    /* Step 5.  Your application code here.
     */
    out_var1 = skeleton(*in_var1);
    sciprint("out_var1 = %g \n\n", out_var1);

    /* Step 6.  Create the output arguments for Scilab engine and assign them.
     *
     *          NOTE: The position of the n_th output variable is nbInputArgument(pvApiCtx) + n ,
     *          (or equivalently, num_of_in_args + n ,) where n counts up from 1 .
     */
    err = createMatrixOfDouble(pvApiCtx, num_of_in_args + 1, 1, 1, &out_var1);
    AssignOutputVariable(pvApiCtx, 1) = num_of_in_args + 1;
    
    /* Step 7.  Return the output arguments to the Scilab engine.
     */
    ReturnArguments(pvApiCtx);

    return 0;
}
