// *****************************************************************************************************
// mex routine to read a segment of data for all MEG channels from a single trial dataset
//
// calling syntax is:
//  
//		data = bw_getCTFData(datasetName, startSample, numSamples);
//
//		datasetName:	name of CTF dataset
//		startSample:	offset from beginning of trial (1st sample = zero!).
//		numSamples:		number of samples to return;
//
// returns
//      data = [numSamples x numSensors] matrix of data in Tesla with gradient of saved data...
//		** this returns primary sensor data only ***
//             
//		(c) Douglas O. Cheyne, 2010-2012  All rights reserved.
//
//		revisions:
//
// ****************************************************************************************************

#include "mex.h"
#include "string.h"
#include "../../bw_libraries/ctflib/headers/datasetUtils.h"
#include "../../bw_libraries/ctflib/headers/BWFilter.h"
#include "../../bw_libraries/ctflib/headers/path.h"
#include "bw_version.h"

double	*chanBuffer;
int		*sampleBuffer;

ds_params		CTF_Data_dsParams;

extern "C" 
{
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] )
{ 
    
	double			*data;
	double			*dPtr;
	char			*dsName;
	char            channelName[256];
	char			megName[4096];
	char			baseName[256];
	char			s[256];
		
	int				buflen;
	int				status;
  	unsigned int 	m;
	unsigned int	n; 
  	char			msg[256];
	
	int				idx;
	int				startSample = 0;
	int				numSamples = 0;
    int             trialNo = 0;
    int             nchans = 1;
    int             channelIndices[MAX_CHANNELS];
    
	double          *val;
	
	int             nbytes;
    int             bytesToSkip;
	double			*dataPtr;	
	

	/* Check for proper number of arguments */
	int n_inputs = 5;
	int n_outputs = 1;
	if ( nlhs != n_outputs | nrhs < n_inputs)
	{
		mexPrintf("bw_getCTFData ver. %.1f (%s) (c) Douglas Cheyne, PhD. 2012. All rights reserved.\n", BW_VERSION, BW_BUILD_DATE); 
		mexPrintf("Incorrect number of input or output arguments\n");
		mexPrintf("Usage:\n"); 
		mexPrintf("   data = bw_getCTFData(datasetName, startSample, numSamples, trialNo, channelIndices) \n");
		mexPrintf("   [datasetName]        - name of dataset\n");
		mexPrintf("   [startSample]        - sample from beginning of trial (1st sample = zero!)\n");
		mexPrintf("   [numSamples]         - sample length to read \n");
        mexPrintf("   [trialNo]            - trial number to read (1st trial = zero!)\n");
		mexPrintf("   [channelIndices]     - [1 x N] integer array of channel indices to read\n");
		
		mexPrintf(" \n");
		return;
	}

	/* get file name */

  	/* Input must be a string. */
  	if (mxIsChar(prhs[0]) != 1)
    		mexErrMsgTxt("Input must be a string.");

  	/* Input must be a row vector. */
  	if (mxGetM(prhs[0]) != 1)
    		mexErrMsgTxt("Input must be a row vector.");

  	/* Get the length of the input string. */
  	buflen = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
  
	/* Allocate memory for input and output strings. */
  	dsName = (char *)mxCalloc(buflen, sizeof(char));

  	/* Copy the string data from prhs[0] into a C string input_buf. */
  	status = mxGetString(prhs[0], dsName, buflen);
  	if (status != 0)
        mexWarnMsgTxt("Not enough space. String is truncated.");

	val = mxGetPr(prhs[1]);
	startSample = (int)*val;	

	val = mxGetPr(prhs[2]);
	numSamples = (int)*val;
    
    val = mxGetPr(prhs[3]);
    trialNo = (int)*val;
    
    if (mxGetM(prhs[4]) != 1 || mxGetN(prhs[4]) < 1)
        mexErrMsgTxt("Input [4] must be a 1 x N row vector of channel numbers.");

    nchans = mxGetN(prhs[4]);
    if (nchans > MAX_CHANNELS)
    {
        mexPrintf("Number of channels %d exceeds maximum number %d ...\n", nchans, MAX_CHANNELS);
        return;
    }

    val = mxGetPr(prhs[4]);
    for (int k=0; k<nchans; k++)
        channelIndices[k] = (int)val[k] - 1;
	
//    mexPrintf("reading %d samples,  %d channels, %d trials  ...\n", numSamples, trialNo, nchans);
    
	// get dataset info
    if ( !readMEGResFile( dsName, CTF_Data_dsParams ) )
    {
		mexPrintf("Error reading res4 file ...\n");
		return;
    }
                  
    if ( nchans >  CTF_Data_dsParams.numChannels)
    {
      mexPrintf("Number of channels %d exceed number in dataset %d ...\n", nchans, CTF_Data_dsParams.numSamples);
      return;
    }
                  
	if ( startSample < 0 ||  startSample + numSamples > CTF_Data_dsParams.numSamples )
	{
		mexPrintf("valid sample range is 0 to %d ...\n", CTF_Data_dsParams.numSamples );
		return;
	}
    
    if ( trialNo+1 > CTF_Data_dsParams.numTrials )
    {
        mexPrintf("valid trial range is 0 to %d ...\n", CTF_Data_dsParams.numTrials );
        return;
    }
    
    if ( trialNo+1 > CTF_Data_dsParams.numTrials )
    {
        mexPrintf("valid trial range is 0 to %d ...\n", CTF_Data_dsParams.numTrials );
        return;
    }

			
	plhs[0] = mxCreateDoubleMatrix(numSamples, nchans, mxREAL);
	data = mxGetPr(plhs[0]);
	
		// mexPrintf("getting data from %s (sample %d to %d)\n", dsName, startSample, startSample+numSamples-1);

	chanBuffer = (double *)malloc( sizeof(double) * CTF_Data_dsParams.numSamples );
	if (chanBuffer == NULL)
	{
		mexPrintf("memory allocation failed for chanBuffer array");
		return;
	}

	sampleBuffer = (int *)malloc( sizeof(int) * CTF_Data_dsParams.numSamples );
	if (sampleBuffer == NULL)
	{
		mexPrintf("memory allocation failed for sampleBuffer array");
		return;
	}
	
	removeFilePath( dsName, baseName);
	baseName[strlen(baseName)-3] = '\0';
	
	sprintf(megName, "%s%s%s.meg4", dsName, FILE_SEPARATOR, baseName );
		
	FILE *fp;
	// open data file and start reading...
	//
	if ( ( fp = fopen( megName, "rb") ) == NULL )
	{
		return;
	}
	
	nbytes = fread( s, sizeof( char ), 8, fp );
	if ( strncmp( s, "MEG4CPT", 7 ) && strncmp( s, "MEG41CP", 7 )  )
	{
		mexPrintf("%s does not appear to be a valid CTF meg4 file\n", megName);
		return;
	}
	
	// num trial bytes per channel 
	int numBytesPerChannel = CTF_Data_dsParams.numSamples * sizeof(int);
    int numBytesPerTrial = CTF_Data_dsParams.numChannels * numBytesPerChannel;
        
    // go to trial offset
    bytesToSkip = trialNo * numBytesPerTrial;
    fseek(fp, bytesToSkip, SEEK_CUR);

    idx=0;    // output array index
	for (int k=0; k<CTF_Data_dsParams.numChannels; k++)
	{
        bool includeChannel = false;
        for (int j=0; j<nchans; j++)
        {
            if (channelIndices[j] == k)
            {
                includeChannel = true;
//                mexPrintf("Reading channel %d ...\n", k );
                break;
            }
        }
        if (includeChannel)
		{
			double thisGain =  CTF_Data_dsParams.channel[k].gain;
			
			// go to sample offset
			int bytesToStart = startSample * sizeof(int);
			fseek(fp, bytesToStart, SEEK_CUR);
			
			nbytes = fread( sampleBuffer, sizeof(int), numSamples, fp);
			
			for (int j=0; j<numSamples; j++)
			{
				 double d = ToHost( (int)sampleBuffer[j] );
				 data[idx++] = d / thisGain;
			}
			bytesToSkip = (CTF_Data_dsParams.numSamples - numSamples - startSample) * sizeof(int);
			fseek(fp, bytesToSkip, SEEK_CUR);			
		}
		else
			fseek(fp, numBytesPerChannel, SEEK_CUR);
		
	}
	


	fclose(fp);
	
	free(chanBuffer);
	free(sampleBuffer);
	
	mxFree(dsName);
	 
	return;
         
}
    
}


