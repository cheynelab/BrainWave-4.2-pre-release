/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//		bw_fitDipole
//
//		Mex function to do a single dipole fit based on dipoleFit.cc, added options for filtering and channel selection
//
//		(c) Douglas O. Cheyne, 2022 All rights reserved.
//
//		revisions:
//      Version 2.0 - update to be able to do multiple dipole fits. Currently limited to max 10 dipoles (= 60 free parameters )
//                    practically more than 3 or 4 will work if start parameters are optimal.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "mex.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#if defined _WIN64 || defined _WIN32
	#include <pthread.h>//pthread library for windows, added by zhengkai
#endif

#include "../../bw_libraries/ctflib/headers/datasetUtils.h"
#include "../../bw_libraries/ctflib/headers/sourceUtils.h"
#include "../../bw_libraries/ctflib/headers/BWFilter.h"
#include "../../bw_libraries/ctflib/headers/path.h"
#include "bw_version.h"

// TOLERANCE values for simplex termination
const int       DEFAULT_MAX_ITER = 300;
const int    	DEFAULT_NUM_PASSES = 2;
const double    DEFAULT_TOLERANCE = 0.001;     // tolerance after NUM_TOL improvements for terminating fit

const double    TOL_SINGULAR = 1.0e-4;          // tolerance for LDL decomp in multiple linear fit
const int       MAX_DIPS = 10;                  // max. number of dipoles

double			**aveTrialData;
ds_params		dsParams;
ds_params		t_dsParams;

// global variables that have to be accessed by the errorFunction
double          *g_measured_field;
double      	*g_forward;
double          **g_field_patterns;

double      	g_totalSumOfSquares;
bool			g_magnetic = false;
bool			g_dewarCoords = false;
bool            g_includeReferenceChannels;

int             g_numDips;
int             g_numMEGChannels;
int             g_selectedGradient;
int             g_constraintFlags[MAX_DIPS];
vectorCart     	g_sphereOrigin;

// arrays to hold dipole positions / orientations
double          g_input_pos[3 * MAX_DIPS];
double          g_input_ori[3 * MAX_DIPS];

extern "C" 
{
	
double 			getDipoleFit( double *fitted_params );
dip_params 		makeDipoleTangential(dip_params dipole, vectorCart sphereOrigin);
double 			fitLinearMoment( double *forward, double *measured, int numChannels );
bool            fitLinearMomentMultiple(int ndips, double *moments, double **field_patterns, double *measured, int numChannels );

	
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] )
{ 
    
	double			*data;
	double			*dataPtr;
	char			*dsName;
	
	int				buflen;
	int				status;
  	unsigned int	m;
	unsigned int	n; 
	char			msg[256];
	double			*val;
    int             ival;
	
	
	double          highPass;
	double          lowPass;
	bool			bidirectional = true;
	double			latency;
	double 			baselineWindowStart = 0.0;
	double 			baselineWindowEnd = 0.0;
	int				useBaseline = 0;
	
	double			moment;
	
	int				numBadChannels = 0;
	char			**badChannelNames;
    bool			badChannelIndex[MAX_CHANNELS];
	double			*buffer;
	
    // array to hold ndips * 6 parameters for Simplex
	double      	param_array[6 * MAX_DIPS];
	double      	delta_array[6 * MAX_DIPS];

    
	double			percentError;
	dip_params  	dipole;
	dip_params		startDip;
    dip_params      previous_dip;
    
	vectorCart  	pvec;
	vectorCart  	orient;
    
	int				sample = 0;
    int             numDips = 1;
    int             idx;
    
	filter_params 	fparams;
	
	int 			max_iterations = DEFAULT_MAX_ITER;
	double 			tolerance = DEFAULT_TOLERANCE;
	int 			num_passes = DEFAULT_NUM_PASSES;
	
	// return values
	double 			*forward;
	double			*position;
	double			*orientation;
	double			*fitMoment;
	double			*fitError;
    double          *moments;
	
	int n_inputs = 9;
	int n_outputs = 5;
	mexPrintf("bw_fitDipole ver. %.1f (%s) (c) Douglas Cheyne, PhD. 2022.\n", BW_VERSION, BW_BUILD_DATE);
	
	if ( nlhs != n_outputs | nrhs < n_inputs)
	{
		mexPrintf("\nincorrect number of input (%d) or output (%d) arguments for bw_makeVS  ...\n", nrhs, nlhs);
		mexPrintf("\nCalling syntax:\n");
		mexPrintf("[forward, positions, orientations, moment, error] = bw_fitDipole(dsName, filter, latency, numDips,\n");
		mexPrintf("       start_positions, start_orientations, sphere, constraintFlags, [numPasses], [tolFactor], [baseline], [badChannelList] )\n");
		mexPrintf("dsName:            - name of raw data (single trial) CTF dataset (will average across trials)\n");
		mexPrintf("filter:            - [highpass lowpass] bandpass in Hz to filter data\n");
        mexPrintf("latency:           - latency to apply fit to (in seconds)\n");
        mexPrintf("numDips:           - number of dipoles to fit\n");
		mexPrintf("start_positions:   - [x, y, z, x2, y2, z2, ... xN, yN, zN] start positions (in cm) for fit (N = numdips)\n");
		mexPrintf("start_orientations:- [xo, yo, zo, xo2, yo2, zo2, ... xoN, yoN, zoN] start orientations for fit (N = numdips)\n");
		mexPrintf("sphere:            - [x, y, z] sphere origin for spherical model (Sarvas)\n");
        mexPrintf("constraintFlags:   - [1 x numDips] integer array where (0=unconstrained, 1=fixed position, 2=mirror symmetric to dipole n-1).\n");
        mexPrintf("dipoleType:        - 1 = current Dipole, 2 = magnetic Dipole\n");
		mexPrintf("\nOptions: (pass empty argument for default)\n");
		mexPrintf("number of passes:  - [numPasses] - number of times to restart the simplex fit (default: 2)\n");
		mexPrintf("tolerance:         - [tolFactor] - tolerance for terminating Simplex fit (default: 0.001)  \n");
		mexPrintf("baselineWindow:    - [bstart bend] - window boundaries for pre-stimulus baseline subtract \n");
		mexPrintf("badChanneList:     - [nchannels x 5 chars] character array of MEG channel names to be excluded.\n");
		mexPrintf("\n returns: [forwardSolution, fittedPositions, fittedOrientations, fittedMoments, percentError] \n");
		return;
	}
    
    int M;
    int N;
    
	///////////////////////////////////
	// get datasest name 
  	if (mxIsChar(prhs[0]) != 1)
		mexErrMsgTxt("Input [0] must be a string.");
 	if (mxGetM(prhs[0]) != 1)
		mexErrMsgTxt("Input [0] must be a row vector.");
   	// Get the length of the input string.
  	buflen = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
  	dsName = (char *)mxCalloc(buflen, sizeof(char));
  	status = mxGetString(prhs[0], dsName, buflen);  	// Copy the string into a C string
 	if (status != 0) 
		mexErrMsgTxt("Not enough space for dsName. String is truncated.");
	
	/////////////////
	// get input parameters

	if (mxGetM(prhs[1]) != 1 || mxGetN(prhs[1]) != 2)
		mexErrMsgTxt("Input [1] must be a row vector [hipass lowpass].");
	dataPtr = mxGetPr(prhs[1]);
	highPass = dataPtr[0];
	lowPass = dataPtr[1];

	val = mxGetPr(prhs[2]);
	latency = *val;
    
    val = mxGetPr(prhs[3]);
    numDips = (int)*val;
    
    if (numDips > MAX_DIPS)
        mexErrMsgTxt("Maximum number of dipoles exceeded (max = 10).");
    
    if (mxGetM(prhs[4]) != 3 || mxGetN(prhs[4]) != numDips)
		mexErrMsgTxt("Input [4] must be a 3 x numdips array.");
	dataPtr = mxGetPr(prhs[4]);
    for (int k=0; k<numDips * 3; k++)
        g_input_pos[k] = dataPtr[k];

	if (mxGetM(prhs[5]) != 3 || mxGetN(prhs[5]) != numDips)
        mexErrMsgTxt("Input [5] must be 3 x numdips array.");
    dataPtr = mxGetPr(prhs[5]);
    for (int k=0; k<numDips * 3; k++)
       g_input_ori[k] = dataPtr[k];
	
	if (mxGetM(prhs[6]) != 1 || mxGetN(prhs[6]) != 3)
		mexErrMsgTxt("Input [6] must be a row vector [sphereX sphereY sphereZ].");
	dataPtr = mxGetPr(prhs[6]);
	g_sphereOrigin.x = dataPtr[0];
	g_sphereOrigin.y = dataPtr[1];
	g_sphereOrigin.z = dataPtr[2];
    
    M = mxGetM(prhs[7]);
    if (M > 0)
    {
        dataPtr = mxGetPr(prhs[7]);
        for (int k=0; k<M; k++)
        {
            g_constraintFlags[k] = (int)dataPtr[k];
            if (g_constraintFlags[k] == 1)
                mexPrintf("* constraining dipole %d to fixed position *\n", k+1);
            if (g_constraintFlags[k] == 2)
                mexPrintf(" * constraining dipole %d to be mirror symmetric to dipole %d *\n", k+1, k);
        }
    }
    
    val = mxGetPr(prhs[8]);
    ival = (int)*val;
    if (ival == 2)
    {
        g_magnetic = true;
        mexPrintf("Fitting magnetic dipole...\n");
    }
    else
    {
        g_magnetic = false;
        mexPrintf("Fitting current dipole...\n");
    }
	// *** optional inputs ***

	if (nrhs > 9)
	{
        M = mxGetM(prhs[9]);
        if (M > 0)
        {
            val = mxGetPr(prhs[9]);
            num_passes = (int)*val;
            mexPrintf("Num passes = %d\n",num_passes);
        }
	}
	
	if (nrhs > 10)
	{
        M = mxGetM(prhs[10]);
        if (M > 0)
        {
            val = mxGetPr(prhs[10]);
            tolerance = *val;
            mexPrintf("Setting Tolerance to: %g\n",tolerance);
        }
	}

    if (nrhs > 11)
    {
        M = mxGetM(prhs[11]);
        if (M > 0)
        {
            dataPtr = mxGetPr(prhs[11]);
            baselineWindowStart = dataPtr[0];
            baselineWindowEnd = dataPtr[1];
            useBaseline = 1;
            mexPrintf("Baselining data from %g to: %g s\n",baselineWindowStart, baselineWindowEnd);
        }
    }
                
   if (nrhs > 12)
   {
        M = mxGetM(prhs[12]);
        if (M > 0)
        {
            if (mxIsCell(prhs[12]) != 1)
                mexErrMsgTxt("channel list input [12] must be cell string array");
            
            mexPrintf("Using custom channel set...\n");
            numBadChannels =mxGetNumberOfElements(prhs[12]);
            badChannelNames = (char **)mxCalloc(numBadChannels,sizeof(char*));
            for (int i=0; i<numBadChannels;i++)
            {
                if (!mxIsChar( mxGetCell(prhs[12],i)))
                    mexErrMsgTxt("channel list must be cell string array");
                else
                {
                    badChannelNames[i] = mxArrayToString(  mxGetCell(prhs[12],i));
    				mexPrintf("excluding bad channel %s \n", badChannelNames[i]);
                }
            }
        }
	}

    
    // add flag to enable constraints (e.g., symmetric dipoles)
    
    // *** assign to numDips to global variable for error function
    g_numDips = numDips;
	
	// Step 1.  Read all the MEG average data. For now saved gradient only.
	
    if ( !readMEGResFile( dsName, dsParams) )
	{
		mexPrintf("Error reading res4 file for %s/n", dsName);
		return;
	}

	mexPrintf("dataset:  %s, (%d trials, %d samples, %d sensors, epoch time = %g to %g s)\n", 
			  dsName, dsParams.numTrials, dsParams.numSamples, dsParams.numSensors, dsParams.epochMinTime, dsParams.epochMaxTime);
	
	// sanity checks here...
	// Get sensor data at latency t
	sample = round( dsParams.sampleRate * (latency + fabs(dsParams.epochMinTime)) );
	
	if (sample < 0 || sample > dsParams.numSamples)
	{
		mexPrintf("Sample %d exceeds data boundaries...\n", sample);
		return;
	}
	
	// setup filter for data sample rate....
	fparams.enable = true;
	if ( highPass == 0.0 )
		fparams.type = BW_LOWPASS;
	else
		fparams.type = BW_BANDPASS;
	fparams.bidirectional = bidirectional;
	fparams.hc = lowPass;
	fparams.lc = highPass;
	fparams.fs = dsParams.sampleRate;
	fparams.order = 4;	//
	fparams.ncoeff = 0;				// init filter
	
	if (build_filter (&fparams) == -1)
	{
		mexPrintf("Could not build filter.  Exiting\n");
		return;
	}

	////////////////////////////////////////////
	// Allocate memory and read average
	////////////////////////////////////////////
	
	
	// disable bad channels - set isSensor to false so it is excluded from fit and update the sensor count
	// have to do this before allocating memory for data arrays
	if (numBadChannels > 0)
	{
		for (int k=0; k < dsParams.numChannels; k++)
		{
			for (int j=0; j<numBadChannels; j++)
			{
				if ( !strncmp( dsParams.channel[k].name, badChannelNames[j], strlen(badChannelNames[j])) )
				{
					dsParams.channel[k].isSensor = false;
					break;
				}
			}
		}
		dsParams.numSensors -= numBadChannels;
		mexPrintf("excluding %d channels from fit (number of remaining sensors = %d)...\n", numBadChannels, dsParams.numSensors);
	}

	
	// allocate memory to read average for all channels - need both sensors and reference channels for filtering.
	aveTrialData = (double **)malloc( sizeof(double *) * dsParams.numChannels );
	if (aveTrialData == NULL)
	{
		mexPrintf("memory allocation failed for trial array");
		return;
	}
	for (int i = 0; i < dsParams.numChannels; i++)
	{
		aveTrialData[i] = (double *)malloc( sizeof(double) * dsParams.numSamples);
		if ( aveTrialData[i] == NULL)
		{
			mexPrintf( "memory allocation failed for trial array" );
			return;
		}
	}
	
	buffer = (double *)malloc( sizeof(double) * dsParams.numSamples );
	if (buffer == NULL)
	{
		mexPrintf("memory allocation failed for buffer array");
		return;
	}
	
    moments = (double *)malloc( sizeof(double) * numDips );
    if (moments == NULL)
    {
        mexPrintf("memory allocation failed for moments array");
        return;
    }
    
    g_field_patterns = (double **)malloc( sizeof(double *) * numDips );
    if (g_field_patterns == NULL)
    {
        mexPrintf("memory allocation failed for g_field_patterns array");
        return;
    }
    for (int i = 0; i < numDips; i++)
    {
        g_field_patterns[i] = (double *)malloc( sizeof(double) * dsParams.numChannels);
        if ( g_field_patterns[i] == NULL)
        {
            mexPrintf( "memory allocation failed for g_field_patterns array" );
            return;
        }
    }
    
	// get average need to read the references as well since
	// forward calculations will include reference gradient corrections
	
	// read average - use saved gradient
	// note readMEGDataAverage will override the local (modified) dsParams struct if passed...
	if ( !readMEGDataAverage( dsName, t_dsParams, aveTrialData, -1, 0) )
	{
		mexPrintf("Error returned from getSensorDataAverage...\n");
		return;
	}
	
	g_selectedGradient = -1;
	
	// Step 2.  Here have to loop through data and filter average
    //			and eliminate channels in bad channels list

	// need to initialize dsParams channel structure with correct head model (sphere origins) since they are used by computeForwardSolution()
	// ** note: need to set sphere origin(s) also for the balancing reference channels even though we don't use them in the fit,
	// since they are used to correct the forward solution for the selected gradient
	
	mexPrintf("Applying filter from %g to %g Hz (bidirectional)\n", highPass, lowPass);

	int bstart = 0;
	int bend = 0;
	int bpts = 0;
	if (useBaseline)
	{
		bstart = dsParams.numPreTrig + (baselineWindowStart * dsParams.sampleRate);
		bend = dsParams.numPreTrig + (baselineWindowEnd * dsParams.sampleRate);
		bpts = bend - bstart;
		mexPrintf("Removing baseline from %g to %g s (sample %d to sample %d) \n", baselineWindowStart, baselineWindowEnd, bstart, bend);
		if (bstart < 0 || bend > dsParams.numSamples || bpts < 1)
		{
			mexPrintf("Error *** invalid baseline parameters ***\n");
			return;
		}
	}
	
	for (int i=0; i < dsParams.numChannels; i++)
	{
		if ( dsParams.channel[i].isSensor || dsParams.channel[i].isBalancingRef )
		{
			// filter data and set sphere for all sensors and reference channels !
			applyFilter( aveTrialData[i], buffer, dsParams.numSamples, &fparams);
			for (int k=0; k< dsParams.numSamples; k++)
				aveTrialData[i][k] = buffer[k];
			
			// remove baseline
			if (useBaseline)
			{
				double mean = 0.0;
				for (int k=0; k< dsParams.numSamples; k++)
					mean += aveTrialData[i][k];
				mean /= (double)bpts;
				for (int k=0; k< dsParams.numSamples; k++)
					aveTrialData[i][k] -= mean;
			}
			
			dsParams.channel[i].sphereX = g_sphereOrigin.x;
			dsParams.channel[i].sphereY = g_sphereOrigin.y;
			dsParams.channel[i].sphereZ = g_sphereOrigin.z;

		}
	}
	
	g_numMEGChannels = dsParams.numSensors;
	g_selectedGradient = dsParams.gradientOrder;

	// set flag to not include references channels in fit (still applies data gradient correction in the forward solution ! )
	g_includeReferenceChannels = false;
	
	g_measured_field = (double *)malloc( sizeof(double) * g_numMEGChannels );
	if ( g_measured_field == NULL)
	{
		mexPrintf( "memory allocation failed for measured_field \n" );
		return;
	}
	
	g_forward = (double *)malloc( sizeof(double) * g_numMEGChannels );
	if ( g_forward == NULL)
	{
		mexPrintf( "memory allocation failed for g_forward \n" );
		return;
	}
	
	mexPrintf("Fitting %d dipole(s) at latency = %lf s (sample %d)\n", numDips, latency, sample);
	
	int channelCount = 0;
	for (int i=0; i < dsParams.numChannels; i++)
	{
		if ( dsParams.channel[i].isSensor )
		{
			g_measured_field[channelCount++] = aveTrialData[i][sample];
		}
	}
	
    // ************************************************************************
    // initialize vertex (param_array) in way that error function can separate into dipoles...
    // param_array[x1,y1,z1, xo1, yo1, zo1,  x2, y2, z2, xo2, yo2, zo2 ...]
    
    int nparams = numDips * 6;
    idx = 0;
    for (int k=0; k<numDips; k++)
    {
        startDip.xpos = g_input_pos[k*3];
        startDip.ypos = g_input_pos[k*3+1];
        startDip.zpos = g_input_pos[k*3+2];
        startDip.xori = g_input_ori[k*3];
        startDip.yori = g_input_ori[k*3+1];
        startDip.zori = g_input_ori[k*3+2];
        
        // also make sure dipoles are tangential for sphere model
        if (g_magnetic)
            dipole = startDip;
        else
        {
            mexPrintf("Adjusting orientation for spherical model...\n");
            dipole = makeDipoleTangential( startDip, g_sphereOrigin);
        }
        
        if (g_constraintFlags[k] == 0)   // add x, y, z dipole position for dipole k to param_array
        {
            param_array[idx++] = dipole.xpos;
            param_array[idx++] = dipole.ypos;
            param_array[idx++] = dipole.zpos;
        }
        else if (g_constraintFlags[k] == 1)   // else adjust num parameters to fit
        {
            nparams -= 3;
        }
        else if (k > 0 && g_constraintFlags[k] == 2)   // else remove y param only from fit
        {
            dipole.ypos = (2.0 * g_sphereOrigin.y) - startDip.ypos;   // make sure starting y is mirrored to previous in list
            // only add x and z to simplex fit
            param_array[idx++] = dipole.xpos;
            param_array[idx++] = dipole.zpos;
            nparams -= 1;
        }
        
        param_array[idx++] = dipole.xori;
        param_array[idx++] = dipole.yori;
        param_array[idx++] = dipole.zori;
        mexPrintf("Starting parameters for Dipole %d: position: %g %g %g, orientation: %g %g %g \n",
                  k+1, dipole.xpos, dipole.ypos, dipole.zpos, dipole.xori, dipole.yori, dipole.zori);
    }
    
    // initialize the vertex delta values, small variations seem to work best if initial guess is close
    idx = 0;
    double posDelta = 2.0;
    double oriDelta = 0.2;
 
    for (int k=0; k<numDips; k++)
    {
        if (g_constraintFlags[k] == 0)   // add x, y, z position for dipole k to fit params
        {
            delta_array[idx++] = posDelta;
            delta_array[idx++] = posDelta;
            delta_array[idx++] = posDelta;
        }
        else if (k > 0 && g_constraintFlags[k] == posDelta)   // else remove y param from fit
        {
            delta_array[idx++] = posDelta;
            delta_array[idx++] = posDelta;
        }

        delta_array[idx++] = oriDelta;
        delta_array[idx++] = oriDelta;
        delta_array[idx++] = oriDelta;
    }
    
	// need global denominator for error term
	g_totalSumOfSquares= 0.0;
	for (int k=0; k<g_numMEGChannels; k++)
        g_totalSumOfSquares += (g_measured_field[k] * g_measured_field[k]);
	
    percentError = getDipoleFit(param_array);
    
    mexPrintf("Fitting data (initial SS = %g, initial error = %g %%) ...\n", g_totalSumOfSquares, percentError);
    
    // ************************************************************************
    // do Simplex fit of all non-linear parameters

	for (int pass=0; pass<num_passes; pass++)
	{
		int iterCount = runSimplexFit(nparams, param_array, delta_array, getDipoleFit, max_iterations, tolerance);
		percentError = getDipoleFit(param_array);
        mexPrintf("   ...pass %d, (%d iterations, final error %g %%)\n", pass+1, iterCount, percentError);
	}
	
    // ************************************************************************
    // ** after fitting non-linear paramters
    // have to recompute the moments since this isn't returned in param_array
    idx = 0;
    for (int k=0; k<g_numDips; k++)
    {
        if (g_constraintFlags[k] == 0)   // add x, y, z position for dipole k to fit params
        {
            dipole.xpos = param_array[idx++];
            dipole.ypos = param_array[idx++];
            dipole.zpos = param_array[idx++];
        }
        else if (g_constraintFlags[k] == 1)   // else use starting x, y, z position
        {
            dipole.xpos = g_input_pos[k*3];
            dipole.ypos = g_input_pos[k*3+1];
            dipole.zpos = g_input_pos[k*3+2];
        }
        else if (k > 0 && g_constraintFlags[k] == 2)   // else use mirrored y pos
        {
            dipole.xpos = param_array[idx++];
            dipole.ypos = (2.0 * g_sphereOrigin.y) - previous_dip.ypos;
            dipole.zpos = param_array[idx++];
        }
        orient.x = param_array[idx++];
        orient.y = param_array[idx++];
        orient.z = param_array[idx++];
        pvec = unitVector(orient);  // make sure orientations are unit vectors
        dipole.xori = pvec.x;
        dipole.yori = pvec.y;
        dipole.zori = pvec.z;
        
        // fitted params may not reflect adjusted orientations - not returned in getDipoleFit ????
        if (!g_magnetic)
        {
            startDip = dipole;
            mexPrintf("Adjusting orientation for spherical model...\n");
            dipole = makeDipoleTangential( startDip, g_sphereOrigin);
        }
        
        previous_dip = dipole;      // remember for constraints
        
        dipole.moment = 1;
        computeForwardSolution(dsParams, dipole, g_field_patterns[k], g_includeReferenceChannels, g_selectedGradient, g_magnetic, g_dewarCoords);
    }

    // refit moments
    fitLinearMomentMultiple(g_numDips, moments, g_field_patterns, g_measured_field, g_numMEGChannels );

    // create output arrays for dipoles
    plhs[1] = mxCreateDoubleMatrix(3*numDips, 1, mxREAL);
    position = mxGetPr(plhs[1]);
    
    plhs[2] = mxCreateDoubleMatrix(3*numDips, 1, mxREAL);
    orientation = mxGetPr(plhs[2]);
    
    plhs[3] = mxCreateDoubleMatrix(numDips, 1, mxREAL);
    fitMoment = mxGetPr(plhs[3]);
    
    // ************************************************************************
    // recompute all forward solutions with fitted moments, correct for polarity
    // and put results in output array
    
    mexPrintf("Fitted dipole parameters: \n");
    idx = 0;
    int posIdx = 0;
    int oriIdx = 0;
    for (int k=0; k<g_numDips; k++)
    {
        if (g_constraintFlags[k] == 0)   // add x, y, z position for dipole k to fit params
        {
            dipole.xpos = param_array[idx++];
            dipole.ypos = param_array[idx++];
            dipole.zpos = param_array[idx++];
        }
        else if (g_constraintFlags[k] == 1)   // else use starting x, y, z position
        {
            dipole.xpos = g_input_pos[k*3];
            dipole.ypos = g_input_pos[k*3+1];
            dipole.zpos = g_input_pos[k*3+2];
        }
        else if (k > 0 && g_constraintFlags[k] == 2)   // else use mirrored y pos
        {
            dipole.xpos = param_array[idx++];
            dipole.ypos = (2.0 * g_sphereOrigin.y) - previous_dip.ypos;
            dipole.zpos = param_array[idx++];
        }
        orient.x = param_array[idx++];
        orient.y = param_array[idx++];
        orient.z = param_array[idx++];
        
        pvec = unitVector(orient);
        dipole.xori = pvec.x;
        dipole.yori = pvec.y;
        dipole.zori = pvec.z;
        
        if (!g_magnetic)
        {
            startDip = dipole;

            dipole = makeDipoleTangential( startDip, g_sphereOrigin);
        }
        
        dipole.moment = moments[k];
 
        // if fitted moment is negative flip dipoles to have positive moments
        if (dipole.moment < 0.0)
        {
            dipole.moment *= -1.0;
            dipole.xori *= -1.0;
            dipole.yori *= -1.0;
            dipole.zori *= -1.0;
        }
        previous_dip = dipole;          // save k-1 dipole parameters for mirror constraint...
   
        // need forward solution for each dipole with correct moment
        computeForwardSolution(dsParams, dipole, g_field_patterns[k], g_includeReferenceChannels, g_selectedGradient, g_magnetic, g_dewarCoords);
        
        mexPrintf("Dipole %d: position: %lf %lf %lf cm, orientation: %lf %lf %lf,  moment: %lf nAm\n",
                    k+1, dipole.xpos, dipole.ypos, dipole.zpos, dipole.xori, dipole.yori, dipole.zori, dipole.moment);
        
        // have to return dipole parameters in one-dimensional output arrays
        position[posIdx++] = dipole.xpos;
        position[posIdx++] = dipole.ypos;
        position[posIdx++] = dipole.zpos;
        orientation[oriIdx++] = dipole.xori;
        orientation[oriIdx++] = dipole.yori;
        orientation[oriIdx++] = dipole.zori;
        fitMoment[k] = dipole.moment;
    }
    
    // sum forward solutions over all dipoles and save in g_forward
    for (int chan=0; chan<g_numMEGChannels; chan++)
        g_forward[chan] = 0.0;
    for (int k=0; k<g_numDips; k++)
    {
        for (int chan=0; chan<g_numMEGChannels; chan++)
            g_forward[chan] += g_field_patterns[k][chan];
    }

	// returns fitted field for valid channels
	plhs[0] = mxCreateDoubleMatrix(g_numMEGChannels, 1, mxREAL);
	forward = mxGetPr(plhs[0]);
    // return fitted field and optimized parameters
    for (int k=0; k< g_numMEGChannels; k++)
        forward[k] = g_forward[k];
    
    // return error term
	plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL);
	fitError = mxGetPr(plhs[4]);
    fitError[0] = percentError;
	
	///////////////////////////////////
	// free temporary arrays for this routine
	mxFree(dsName);
	
	for (int k=0; k< dsParams.numChannels; k++)
		free(aveTrialData[k]);
	free(aveTrialData);
    
    for (int k=0; k< numDips; k++)
        free(g_field_patterns[k]);
    free(g_field_patterns);
    
	free(g_measured_field);
	free(g_forward);
    free(buffer);
    free(moments);
	
	return;
         
}

// error function to pass to Simplex
double getDipoleFit( double *fitted_params )
{
	dip_params      dip;
	dip_params		corrected_dip[MAX_DIPS];
    dip_params      previous_dip;
    
	vectorCart      pvec;
	vectorCart		pos;
	vectorCart      orient;
	double 			percentError;
    double          dipMoments[MAX_DIPS];
    int             idx;

    // compute forward solution for each dipole
    idx = 0;
    for (int k=0; k<g_numDips; k++)
    {
        if (g_constraintFlags[k] == 0)   // if fitting position get x, y, z position from fitted params
        {
            dip.xpos = fitted_params[idx++];
            dip.ypos = fitted_params[idx++];
            dip.zpos = fitted_params[idx++];
        }
        else if (g_constraintFlags[k] == 1)   // else use starting x, y, z position
        {
            dip.xpos = g_input_pos[k*3];
            dip.ypos = g_input_pos[k*3+1];
            dip.zpos = g_input_pos[k*3+2];
        }
        else if (k > 0 && g_constraintFlags[k] == 2)   // else use mirrored y pos
        {
            dip.xpos = fitted_params[idx++];
            dip.ypos = (2.0 * g_sphereOrigin.y) - previous_dip.ypos;
            dip.zpos = fitted_params[idx++];
        }

        orient.x = fitted_params[idx++];
        orient.y = fitted_params[idx++];
        orient.z = fitted_params[idx++];

        pvec = unitVector(orient);
        dip.xori = pvec.x;
        dip.yori = pvec.y;
        dip.zori = pvec.z;
        
        previous_dip = dip; // save k-1 dipole parameters for mirror constraint...

        if (g_magnetic)
            corrected_dip[k] = dip;
        else
            corrected_dip[k] = makeDipoleTangential( dip, g_sphereOrigin);
        
        // get normalized forward solution (in nanoAmp-meters) using current dip parameters
        // patterns returned in array g_field_patterns
        corrected_dip[k].moment = 1;
        computeForwardSolution(dsParams, corrected_dip[k], g_field_patterns[k], g_includeReferenceChannels, g_selectedGradient, g_magnetic, g_dewarCoords);
    }

    // fit all moments as linear weights
    fitLinearMomentMultiple(g_numDips, dipMoments, g_field_patterns, g_measured_field, g_numMEGChannels );

    // recompute all forward solutions with fitted moments
    for (int k=0; k<g_numDips; k++)
    {
        corrected_dip[k].moment = dipMoments[k];
        computeForwardSolution(dsParams, corrected_dip[k], g_field_patterns[k], g_includeReferenceChannels, g_selectedGradient, g_magnetic, g_dewarCoords);
    }

    // sum forward solutions over all dipoles and save in g_forward
    for (int chan=0; chan<g_numMEGChannels; chan++)
        g_forward[chan] = 0.0;
    for (int k=0; k<g_numDips; k++)
    {
        for (int chan=0; chan<g_numMEGChannels; chan++)
            g_forward[chan] += g_field_patterns[k][chan];
    }
    
	// compute error term normalized to percent of total sum of squares
	double SSError = 0.0;
	for (int chan=0; chan<g_numMEGChannels; chan++)
	{
		double err = g_forward[chan] - g_measured_field[chan];
		SSError += err * err;
	}
    
	percentError = (SSError / g_totalSumOfSquares ) * 100.0;

	return( percentError );
}

// fit moment as linear parameter
// for one dipole this simplifies to simple linear regression where
// moment = forward.measured / forward.forward
double fitLinearMoment( double *forward, double *measured, int numChannels )
{
	double sumForward = 0.0;
	for (int chan=0; chan<numChannels; chan++)
		sumForward += forward[chan] * forward[chan];
	
	double sumForwardMeasured = 0.0;
	for (int chan=0; chan<numChannels; chan++)
		sumForwardMeasured += forward[chan] * measured[chan];
	
	return( sumForwardMeasured/sumForward );
}

// fit moment as linear parameter for multiple dipoles using LDL' decomposition for multiple linear equations.
// Based on original SDIP code .. Boolean TFit_Object::DoLinearFitOfMoments(void))
// * note for spatiotemporal fits the original code just looped over all calculations for t samples
// (i.e., for spatiotemporal fit one could just call this routine in a loop over t samples to get moment over time).

bool fitLinearMomentMultiple(int ndips, double *moments, double **field_patterns, double *measured, int numChannels )
{
    double  c[MAX_DIPS+1][MAX_DIPS+1];
    double  sum_field;
    double  temp;
    double  sum;
    double  b[MAX_DIPS];
    
    c[ndips][ndips] = g_totalSumOfSquares;
    for (int i=0; i< ndips; i++)
    {
        sum_field = 0.0;
        for (int chan=0; chan<numChannels; chan++)
            sum_field += field_patterns[i][chan]  * measured[chan];
        
        c[ndips][i] = sum_field;
        
        for (int k=0; k<=i; k++)
        {
            sum_field = 0.0;
            for (int chan=0; chan<numChannels; chan++)
                sum_field += field_patterns[i][chan]  * field_patterns[k][chan];
            c[i][k] = sum_field;
            
        }
        b[i] = sum_field * TOL_SINGULAR;
    }
    
    // do in-place LDL' decomposition of matrix c
    for (int i=0; i< ndips; i++)
    {
        temp = c[i][i];
        // check if c is singular
        if (temp <= b[i])
        {
            for (int i=0; i<ndips; i++)
                moments[i] = 0.0;
            return(false);
        }
        for (int j=i+1; j<ndips+1; j++)
        {
            sum = c[j][i];
            c[j][i] = sum / temp;
            for (int k=i+1; k<=j; k++)
                c[j][k] -= sum * c[k][i];
        }
    }
    
    b[ndips-1] = c[ndips][ndips-1];     // weight for first dipole
    
    // if more than one dipole
    if (ndips > 1)
    {
        for (int k=ndips-2; k>=0; k--)
        {
            sum = 0.0;
            for (int j=k+1; j<ndips; j++)
                sum += b[j] * c[j][k];
            b[k] = c[ndips][k] - sum;
        }
    }
    
    for (int i=0; i<ndips; i++)
        moments[i] = b[i];
    
    return (true);
}

dip_params makeDipoleTangential(dip_params dipole, vectorCart sphereOrigin)
{
	vectorCart	a;
	vectorCart	pos;
	vectorCart	pvec;
	vectorCart	orient;
	
	dip_params	result;
	
	pos.x = dipole.xpos;
	pos.y = dipole.ypos;
	pos.z = dipole.zpos;
	orient.x = dipole.xori;
	orient.y = dipole.yori;
	orient.z = dipole.zori;
	
	// make dipole tangential in sphere (orthogonal to position vector)
	a = subtractVectors( pos, sphereOrigin);   	// dipole position vector relative to sphere origin
	pvec = makeOrthogonalTo( orient, a );		// make orient orthognal to a
	pvec = unitVector( pvec );
	
	result.xpos = pos.x;
	result.ypos = pos.y;
	result.zpos = pos.z;
	result.xori = pvec.x;
	result.yori = pvec.y;
	result.zori = pvec.z;
	
	return(result);
}

	
	
}


