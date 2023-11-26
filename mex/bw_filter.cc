// *************************************
// mex routine to read single trial dataset and return the average waveforms for all sensor channels
//
// calling syntax is:
// [fdata] = bw_filter( data, hipass, lowpass, sampleRate, order, bidirectional );
//
// returns
//      fdata = [1 x nsamples] vector of filtered data
//
//		(c) Douglas O. Cheyne, 2004-2010  All rights reserved.
//
//		1.0  - first version
//		1.2	 - modified to be consistent with ctf_BWFilter.cc - fixed order and adds bandreject option
//
//      version 3.3 Dec 2016 - modified to adjust order for bidirectional - more consistent with CTF DataEditor filter
//
//      version 4 - Nov 2023 - modified to filter multiple channels at once
// ************************************

#include "mex.h"
#include "string.h"
#include "../../bw_libraries/ctflib/headers/datasetUtils.h"
#include "../../bw_libraries/ctflib/headers/BWFilter.h"
#include "bw_version.h"

double	*buffer;
double  *bufferOut;

extern "C" 
{
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] )
{ 
    
	double			*data;
	double			*fdata;
	
	double			*dataPtr;
	double          *val;
	
	double			sampleRate;
	double          highPass;
	double          lowPass;
	int				numSamples;
    int             numChannels;
    
	int				filterOrder = 4;
    int             maxOrder = 8;           // will limit coeffs to 4th order
	bool			bidirectional = true;
	bool			bandreject = false;
    int             idx;
    int             idx2;
    
	filter_params 	fparams;

	/* Check for proper number of arguments */
	int n_inputs = 3;
	int n_outputs = 1;
	if ( nlhs != n_outputs | nrhs < n_inputs)
	{
		mexPrintf("bw_filter ver. %.1f (%s) (c) Douglas Cheyne, PhD. 2010. All rights reserved.\n", BW_VERSION, BW_BUILD_DATE); 
		mexPrintf("Incorrect number of input or output arguments\n");
		mexPrintf("Usage:\n"); 
		mexPrintf("   [fdata] = bw_filter( data, sampleRate, [hipass lowpass], {options}  )\n");
		mexPrintf("   [data] must be nsamples x nchannels array\n");
		mexPrintf("   [sampleRate] sample rate of data.\n");
		mexPrintf("   [highPass lowPass] high  and low pass cutoff frequency in Hz for bandpass. Enter 0 for highPass for lowPass only.\n");
		mexPrintf("Options:\n");
		mexPrintf("   [order]           - specify filter order. (4th order recommended)\n");
		mexPrintf("   [bidirectional]   - if true filter is bidirectional (two-pass non-phase shifting). Default = true\n");
		mexPrintf("   [bandreject]      - if true filter is band-reject. Default = band-pass\n");
        mexPrintf("   [fdata]           - output is nsamples x nchannels array\n");
		mexPrintf(" \n");
		return;
	}
       
    // input array must be nsamples x nchannels
    // in order to loop across channels for filtering
    
    numChannels = (int)mxGetN(prhs[0]);
    numSamples = (int)mxGetM(prhs[0]);
    
    data = mxGetPr(prhs[0]);
    
	val = mxGetPr(prhs[1]);
	sampleRate = *val;
	
	if (mxGetM(prhs[2]) != 1 || mxGetN(prhs[2]) != 2)
		mexErrMsgTxt("Input [2] must be a row vector [hipass lowpass].");
	dataPtr = mxGetPr(prhs[2]);
	highPass = dataPtr[0];
	lowPass = dataPtr[1];

	if (nrhs > 3)
	{
		val = mxGetPr(prhs[3]);
		filterOrder = (int)*val;
	}
	
	if (nrhs > 4)
	{
		val = mxGetPr(prhs[4]);
		bidirectional = (int)*val;
	}

	if (nrhs > 5)
	{
		val = mxGetPr(prhs[5]);
		bandreject = (int)*val;
	}

    // create output array - has to be filled as if flat array
	plhs[0] = mxCreateDoubleMatrix(numSamples, numChannels, mxREAL);
	fdata = mxGetPr(plhs[0]);
	
//	mexPrintf("Read %d samples, BW = %g %g Hz, sampleRate = %g, order = %d, bidirectional = %d\n", 
//			  numSamples, highPass, lowPass, sampleRate, filterOrder, bidirectional);
	
	fparams.enable = true;
	
	if ( highPass == 0.0 && lowPass == 0.0)
	{
		mexPrintf("invalid filter settings\n");
		return;
	}
	
	if ( highPass == 0.0 )
		fparams.type = BW_LOWPASS;
	else
		fparams.type = BW_BANDPASS;
	
	if (bandreject)
	{
		if ( highPass == 0.0 )
		{
			mexPrintf("high-pass frequency must be specified for band-reject filter.\n");
			return;
		}
		else
			fparams.type = BW_BANDREJECT;
	}
	else
	{
		if ( highPass == 0.0 )
			fparams.type = BW_LOWPASS;
		else if ( lowPass == 0.0 )
			fparams.type = BW_HIGHPASS;
		else
			fparams.type = BW_BANDPASS;
	}
    
    if ( filterOrder > maxOrder)
	{
		mexPrintf("filter order too high...\n");
		return;
	}
	
    if (bidirectional )
    {
        double t = filterOrder / 2.0;
        filterOrder = round(t);   // in case non multiple of 2
        if (filterOrder < 1)
            filterOrder = 1;
    }
    
    
	fparams.bidirectional = bidirectional;
	fparams.hc = lowPass;
	fparams.lc = highPass;
	fparams.fs = sampleRate;
	fparams.order = filterOrder;
	fparams.ncoeff = 0;
	
	
	if (build_filter (&fparams) == -1)
	{
		mexPrintf("memory allocation failed for trial array\n");
		return;
	}
	
	
	buffer = (double *)malloc( sizeof(double) * numSamples );
	if (buffer == NULL)
	{
		mexPrintf("memory allocation failed for buffer array");
		return;
	}
    
    bufferOut = (double *)malloc( sizeof(double) * numSamples );
    if (bufferOut == NULL)
    {
        mexPrintf("memory allocation failed for bufferOut array");
        return;
    }
    
    idx = 0;
    idx2 = 0;
    for (int j=0; j<numChannels; j++)
    {
        for (int k=0; k< numSamples; k++)
            buffer[k] = data[idx++] ;  //  multidimensional arrays are passed to mex function as one dimensinal arrays
        
        
        applyFilter( buffer, bufferOut, numSamples, &fparams);
        for (int k=0; k< numSamples; k++)
            fdata[ idx2++ ] = bufferOut[k];
    }
    
    
	free(buffer);
    free(bufferOut);
    
	return;
         
}
    
}


