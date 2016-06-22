#include <mex.h>
#include <dlib/optimization.h>

/* This is a MEX-file for MATLAB */

/* ==== License ====
   
   Copyright (c) [2014] [Karlsruhe Institute of Technology
                         Institute of Engineering Mechanics]
   
   Permission is hereby granted, free of charge, to any person obtaining a
   copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to permit
   persons to whom the Software is furnished to do so, subject to the
   following condition:
   
     * The above copyright notice and this permission notice shall be
       included in all copies or substantial portions of the Software.
   
   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
   NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
   USE OR OTHER DEALINGS IN THE SOFTWARE.
   
   
   For detailed documentation regarding execution refer to the MATLAB file "bobyqa.m".
*/

/* 
   In dlib, the general purpose solvers optimize functions that
   take a column vector as input and return a double.
*/
typedef dlib::matrix<double,0,1> column_vector;
 
 /* Objective function */
class cObj
{
public:
    cObj ( char sFunName, int nN )
    {
        *sObjFunName = sFunName;
		nN_opt = nN;
    }

    double operator() ( const column_vector vX_cv ) const
    {
		mxArray *lhs[1], *rhs[2];
		
		lhs[0] = mxCreateDoubleScalar(0.0);
		rhs[0] = mxCreateString(sObjFunName);
		rhs[1] = mxCreateDoubleMatrix(nN_opt, 1, mxREAL);
		
		double *vX_dbl = mxGetPr(rhs[1]);
		for (int i=0; i < nN_opt; i++)
		{
			vX_dbl[i] = vX_cv(i);
		}
		
		mexCallMATLAB(1, lhs, 2, rhs, "bobyqa");
		
		return *mxGetPr(lhs[0]);
    }

private:
    char *sObjFunName;
	int  nN_opt;
}; 

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

    /* Input variables */
    char   *sObjFunName = mxArrayToString(prhs[0]);
	double *vX_in       = mxGetPr(prhs[1]);
	double *nN_dbl      = mxGetPr(prhs[2]);
	double *nN_pt       = mxGetPr(prhs[3]);
	double *vX_in_l     = mxGetPr(prhs[4]);
	double *vX_in_u     = mxGetPr(prhs[5]);
    double *nRho_beg    = mxGetPr(prhs[6]);
    double *nRho_end    = mxGetPr(prhs[7]);
    double *nMaxFunEval = mxGetPr(prhs[8]);
	
	/* Number of elements in optimization vector */
	*nN_dbl = *nN_dbl + 0.5;
	const int nN = (int)*nN_dbl;
	
	/* Output variables */
    plhs[0] = mxCreateDoubleScalar(0.0);
    plhs[1] = mxCreateDoubleMatrix(nN, 1, mxREAL);
	
    double *nF_opt = mxGetPr(plhs[0]);
	double *vX_out = mxGetPr(plhs[1]);
	
	/* Local variables */
	column_vector vX_opt(nN);
	column_vector vX_l(nN);
	column_vector vX_u(nN);
	
	/* Input */
    for (int i=0; i < nN; i++)
	{
		vX_opt(i) = vX_in[i];
		vX_l(i)   = vX_in_l[i];
		vX_u(i)   = vX_in_u[i];
	}
	
    /* Call BOBYQA algorithm */
	cObj fObj = cObj(*sObjFunName, nN);
	
	try
	{
		*nF_opt = dlib::find_min_bobyqa(fObj,          // objective function
										vX_opt,        // optimization vector
										*nN_pt,        // number of interpolation points
										vX_l,          // lower bound constraint
										vX_u,          // upper bound constraint
										*nRho_beg,     // initial trust region radius
										*nRho_end,     // stopping trust region radius
										*nMaxFunEval   // max number of objective function evaluations
										);
	}
	
	catch(dlib::bobyqa_failure& eBobyqa_failure)
    {
		const char *sBobyqa_failure = eBobyqa_failure.what();
		
		*nF_opt = fObj(vX_opt);
		
		mexPrintf("\n");
		mexPrintf(sBobyqa_failure);
		mexPrintf("\n\n");
	}
	
	/* Output */
    for (int i=0; i < nN; i++)
	{
		vX_out[i] = vX_opt(i);
	}
}