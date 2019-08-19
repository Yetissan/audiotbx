/* autbx_emd function prototype:
 *
 * function y = autbx_emd(x, order, iters, locality);
 */

#include "api_scilab.h"
#include "Scierror.h"
#include "malloc.h"
#include "sciprint.h"
#include <localization.h>

#include "EmpiricalModeDecomposition.h"

int sci_autbx_emd(char *fname, unsigned long fname_len)
{
    SciErr      err;
    int         num_of_in_args;
    int         num_of_out_args;
    int         **pvaraddrs = NULL;
    int         **pvartypes = NULL;
    int         i;
    
    double      *x;
    int         x_rows, x_cols, x_len;
    double      *order;
    int         order_rows, order_cols;
    double      *iters;
    int         iters_rows, iters_cols;
    double      *locality;
    int         locality_rows, locality_cols;
    
    emdData     data;
    
    double      *y;
        
    /* Constants used in step 1. */
    const int   MIN_NUM_OF_IN_ARGS      = 4;
    const int   MAX_NUM_OF_IN_ARGS      = 4;
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
    pvaraddrs = (int**)malloc(sizeof(int*) * num_of_in_args);
    if (pvaraddrs == NULL)
    {
        Scierror(999, "Failed to allocate memory. \n\n");
        goto cleanup_1;
    }

    for (i = 0; i < num_of_in_args; ++i)
    {
        err = getVarAddressFromPosition(pvApiCtx, i + 1, &pvaraddrs[i]);
        if (err.iErr)
        {
            printError(&err, 0);
            goto cleanup_1;
        }
    }

    /* Step 3.  Check the data type of each variable in `pvaraddrs[]'.
     */
    pvartypes = (int**)malloc(sizeof(int*) * num_of_in_args);
    if (pvartypes == NULL)
    {
        Scierror(999, "Failed to allocate memory. \n\n");
        goto cleanup_2;
    }

    for (i = 0; i < num_of_in_args; ++i)
    {
        err = getVarType(pvApiCtx, pvaraddrs[i], pvartypes[i]);
        if (err.iErr)
        {
            printError(&err, 0);
            goto cleanup_2;
        }
    }
    
    for (i = 0; i < num_of_in_args; ++i)
    {
        if (*pvartypes[i] != sci_matrix)
        {
            Scierror(999, "Input argument #%d has an invalid datatype %d, sci_matrix expected. \n\n", i + 1, pvartypes[i]);
            goto cleanup_2;
        }
    }

    /* Step 4.  Check if the relevant input argument is complex or not, if needed. While
     *          dealing with integer, further checks should be done on the precision of
     *          the integer. Check the size of matrices and vectors. Check consistency of
     *          some input variables.
     */
    for (i = 0; i < num_of_in_args; ++i)
    {
        if (isVarComplex(pvApiCtx, pvaraddrs[i]))
        {
            Scierror(999, "Input argument %d must be a real matrix of doubles. \n\n", i + 1);
            goto cleanup_2;
        }
    }
    
    /* Check `x' */
    err = getMatrixOfDouble(pvApiCtx, pvaraddrs[0], &x_rows, &x_cols, &x);
    if (err.iErr)
    {
        printError(&err, 0);
        goto cleanup_2;
    }
    if ((x_rows == 0) || (x_cols == 0))
    {
        Scierror(999, "Wrong size of input argument 1. \n\n");
        goto cleanup_2;
    }
    x_len = x_rows * x_cols;

    /* Check `order' */
    err = getMatrixOfDouble(pvApiCtx, pvaraddrs[1], &order_rows, &order_cols, &order);
    if (err.iErr)
    {
        printError(&err, 0);
        goto cleanup_2;
    }
    if ((order_rows != 1) || (order_cols != 1))
    {
        Scierror(999, "Wrong size of input argument 2. \n\n");
        goto cleanup_2;
    }
    
    /* Check `iterations' */
    err = getMatrixOfDouble(pvApiCtx, pvaraddrs[2], &iters_rows, &iters_cols, &iters);
    if (err.iErr)
    {
        printError(&err, 0);
        goto cleanup_2;
    }
    if ((iters_rows != 1) || (iters_cols != 1))
    {
        Scierror(999, "Wrong size of input argument 3. \n\n");
        goto cleanup_2;
    }
    
    /* Check `locality' */
    err = getMatrixOfDouble(pvApiCtx, pvaraddrs[3], &locality_rows, &locality_cols, &locality);
    if (err.iErr)
    {
        printError(&err, 0);
        goto cleanup_2;
    }
    if ((locality_rows == 0) || (locality_cols == 0))
    {
        Scierror(999, "Wrong size of input argument 4. \n\n");
        goto cleanup_2;
    }

    /* Step 5.  Your application code here.
     */
    pvartypes = (int**)malloc(sizeof(int*) * num_of_in_args);
    if (pvartypes == NULL)
    {
        Scierror(999, "Failed to allocate memory. \n\n");
        goto cleanup_2;
    }

    emdCreate(&data, x_len, (int)*order, (int)*iters, (int)*locality);
    emdDecompose(&data, (float*)x);


    /* Step 6.  Create the output arguments for Scilab engine and assign them.
     *
     *          NOTE: The position of the n_th output variable is nbInputArgument(pvApiCtx) + n ,
     *          (or equivalently, num_of_in_args + n ,) where n counts upward from 1 .
     */
    err = createMatrixOfDouble(pvApiCtx, num_of_in_args + 1, data.order, data.size, y);
    if (err.iErr)
    {
        printError(&err, 0);
        goto cleanup_2;
    }
    
    for (i = 0; i < data.order; ++i)
    {
        memcpy(y + i * (data.size), (data.imfs) + i * (data.size), data.size);
    }
    
    AssignOutputVariable(pvApiCtx, 1) = num_of_in_args + 1;
    emdClear(&data);
    
    /* Step 7.  Return the output arguments to the Scilab engine.
     */
    ReturnArguments(pvApiCtx);

cleanup_2:
    free(pvartypes);
cleanup_1:
    free(pvaraddrs);
    return 0;
}

