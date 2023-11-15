//////////////////////////////////////////////////////////////////////////////////////////////////////
//		cleanDs
//		program to fix problem with datasets from new Acq with M and E test channel names
//		(currently assumes these are in a contiguous block after MEG channels)
//
// 		D. Cheyne  Dec, 2020
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include "string.h"
#include "../../bw_libraries/ctflib/headers/datasetUtils.h"
#include "../../bw_libraries/ctflib/headers/BWFilter.h"
#include "../../bw_libraries/ctflib/headers/path.h"
#include "bw_version.h"

const double	VERSION_NO = 4.2;					// Version
char const      *RELEASE_DATE = "Oct-14-2023";


ds_params         dsParams;
ds_params         newParams;


double            *channelData;
double            *channelData2;
double            **DE_trialArray;
double            **DE_trialArray2;
extern "C"
{
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] )
{
    
    char            *dsName;
    char            *newDsName;
    
    char            dsBaseName[4096];
    char            newDsBaseName[4096];
    
    double          *badChannels;
    int             badChannelIndices[MAX_CHANNELS];
    int             numBadChannels;
    
    double          *badTrials;
    int             badTrialIndices[MAX_CHANNELS];
    int             numBadTrials;
    
    int             buflen;
    int             status;
    unsigned int    m;
    unsigned int    n;
    char            msg[256];
    char            tStr[256];
    char            cmd[256];

    
    int             vectorLen;
    double          *val;
    double          *dataPtr;
    
    double          *err;
    FILE            *fp2;

    double          highPass;
    double          lowPass;
    int             downSample = 1;
    bool            filterData = false;
    bool            useAllData = true;
    int             startSample;
    int             endSample;
    

    filter_params   fparams;
    
    if ( nlhs != 1 | nrhs < 6)
    {
        mexPrintf("bw_CTFNewDs ver. %.1f (%s) (c) Douglas Cheyne, PhD. 2010. All rights reserved.\n", BW_VERSION, BW_BUILD_DATE);
        mexPrintf("Incorrect number of input or output arguments\n");
        mexPrintf("Usage:\n");
        mexPrintf("   err = bw_CTFNewDs(dsName, newDsName, filterData, [highPass lowPass], downSample, {badChannelList} )\n\n");
        mexPrintf("   [dsName]                  - name of raw data (single trial) CTF dataset to epoch\n");
        mexPrintf("   [newDsName]               - output name for epoched data (use *.ds extension!) \n");
        mexPrintf("   [filterData]              - flag indicating whether to filter data prior to saving. if set to 0 next arguments are ignored\n");
        mexPrintf("   [highPass lowPass]        - row vector specifying prefilter data with high pass and low pass filter in Hz. (enter dummy values if not filtering)\n");
        mexPrintf("   [badChannels]             - row vector of channel indices (first channel = 0) specifying channels to exclude. Enter [] if no bad channels to exclude\n");
        mexPrintf("   [badTrials]               - row vector of trial indices (first trial = 0) to exclude. Enter [] if no trials to exclude\n");
        mexPrintf("   [startSample endSample]   - row vector to specify start and end sample for each trial. Enter [] to save all samples.\n");
        mexPrintf(" \n");
        return;
    }

	bool gotData = false;
	
    /* ================== */
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
    
      if (mxIsChar(prhs[1]) != 1)
        mexErrMsgTxt("Input must be a string.");
    
      /* Input must be a row vector. */
      if (mxGetM(prhs[1]) != 1)
        mexErrMsgTxt("Input must be a row vector.");
    
      /* Get the length of the input string. */
      buflen = (mxGetM(prhs[1]) * mxGetN(prhs[1])) + 1;
    
    /* Allocate memory for input and output strings. */
      newDsName = (char *)mxCalloc(buflen, sizeof(char));
    
    /* Copy the string data from prhs[0] into a C string input_buf. */
    status = mxGetString(prhs[1], newDsName, buflen);
    if (status != 0)
        mexWarnMsgTxt("Not enough space. String is truncated.");

    dataPtr = mxGetPr(prhs[2]);
    filterData = (int)dataPtr[0];

    if (mxGetM(prhs[3]) != 1 || mxGetN(prhs[3]) != 2)
        mexErrMsgTxt("Input [3] must be a row vector [hipass lowpass].");
    dataPtr = mxGetPr(prhs[3]);
    highPass = dataPtr[0];
    lowPass = dataPtr[1];

    // optional arguments
    
    numBadChannels = 0;
    if (mxGetM(prhs[4]) == 1 )
    {
        badChannels = mxGetPr(prhs[4]);
        numBadChannels = mxGetN(prhs[4]);
        for (int k=0; k<numBadChannels; k++)
        {
            double dval = badChannels[k];
            badChannelIndices[k]= (int)dval;
        }
    }
                  
    numBadTrials = 0;
    if (mxGetM(prhs[5]) == 1 )
    {
      badTrials = mxGetPr(prhs[5]);
      numBadTrials = mxGetN(prhs[5]);
      for (int k=0; k<numBadTrials; k++)
      {
          double dval = badTrials[k];
          badTrialIndices[k]= (int)dval;
      }
    }
    
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    err = mxGetPr(plhs[0]);
    err[0] = 0;
	
    // first make sure we can create the new dataset folder
    
    sprintf(msg,"creating new dataset %s...\n", newDsName );
    mexPrintf(msg);
        
    int result;
#if _WIN32||WIN64
    result = mkdir(newDsName);
#else
    result = mkdir(newDsName, S_IRUSR | S_IWUSR | S_IXUSR );
#endif

    int errCode;

    if ( result != 0 )
    {
        mexPrintf("** overwriting existing directory %s ...\n", newDsName);
        sprintf(cmd,"rm -r %s",newDsName);
        errCode = system(cmd);
#if _WIN32||WIN64
        mkdir(newDsName);
#else
        mkdir(newDsName, S_IRUSR | S_IWUSR | S_IXUSR );
#endif
    }
    
    // make sure directory is readable
    sprintf(cmd,"chmod a+rX %s",newDsName);
    errCode = system(cmd);
    
        
    // get dataset info
    if ( !readMEGResFile( dsName, dsParams ) )
    {
        mexPrintf("Error reading res4 file ...\n");
        err[0] = -1;
        mxFree(dsName);
        mxFree(newDsName);
        return;
    }
    
    
    // make copy of original res4 - this might be modified
    if ( !readMEGResFile( dsName, newParams ) )
    {
        mexPrintf("Error reading res4 file ...\n");
        err[0] = -1;
        mxFree(dsName);
        mxFree(newDsName);
        return;
    }
    
    startSample = 0;
    endSample = dsParams.numSamples;
    
    // optional specify subset of sample range
    if (mxGetM(prhs[6]) == 1 && mxGetN(prhs[6]) == 2)
    {
        dataPtr = mxGetPr(prhs[6]);
        startSample = (int)dataPtr[0];
        endSample = (int)dataPtr[1];
    }
        
    // check valid range
    if (startSample < 0 || endSample > dsParams.numSamples || startSample > endSample)
    {
        mexPrintf("Invalid values for start and end samples. Exiting\n");
        err[0] = -1;
        mxFree(dsName);
        mxFree(newDsName);
        return;
    }
        
    // otherwise input and output wil be same size...
    if (startSample > 1 || endSample < dsParams.numSamples)
    {
        double dwel = (1.0 / newParams.sampleRate);
        newParams.numSamples = endSample - startSample + 1;
        newParams.trialDuration = newParams.numSamples * dwel;
         
        // modify epoch min max time based on truncation - however,
        // epochMinTime and epochMaxTime don't get saved in res4 - can only change # of preTrigPts
        // also CTF datasets can't have negative preTrigPts - i.e., epochMinTime always starts at zero or is negative.
        // so in principle startTime can't exceed original preTrigPts (t = 0.0)
        // instead of flagging error here should be prevented by calling routine or markerTimes will be incorrect.
        newParams.numPreTrig = dsParams.numPreTrig - startSample;
        if (newParams.numPreTrig < 0)   // not allowed in CTF .ds
            newParams.numPreTrig = 0;
        
        // only relevant for printout below - these values don't get saved in .res4...
        newParams.epochMinTime = -newParams.numPreTrig * dwel;
        newParams.epochMaxTime = (newParams.numSamples - newParams.numPreTrig - 1) * dwel;

    }
    
    sprintf(msg, "Reading dataset: %s\n", dsName);
    mexPrintf(msg);
    sprintf(msg,"Number of Channels %d\n", dsParams.numChannels);
    mexPrintf(msg);
    sprintf(msg,"Number of Primary MEG Channels %d\n", dsParams.numSensors);
    mexPrintf(msg);
    sprintf(msg,"Number of samples: %d\n", dsParams.numSamples);
    mexPrintf(msg);
    sprintf(msg,"Sample Rate: %.1f\n", dsParams.sampleRate);
    mexPrintf(msg);
    sprintf(msg,"Number of trials: %d\n", dsParams.numTrials);
    mexPrintf(msg);
    sprintf(msg,"Trial Duration (s): %.4f\n", dsParams.trialDuration);
    mexPrintf(msg);
    sprintf(msg,"Epoch Min Time (s): %.4f\n", dsParams.epochMinTime);
    mexPrintf(msg);
    sprintf(msg,"Epoch Max Time (s): %.4f\n", dsParams.epochMaxTime);
    mexPrintf(msg);
    fparams.enable = false;
    if (filterData)
    {
        sprintf(msg,"Filtering data from %.1f to %.1f Hz (bidirectional, order = 8)\n", highPass, lowPass);
        mexPrintf(msg);
        if ( highPass == 0.0 )
            fparams.type = BW_LOWPASS;
        else
            fparams.type = BW_BANDPASS;
    
        fparams.bidirectional = true;
        fparams.hc = lowPass;
        fparams.lc = highPass;
        fparams.fs = dsParams.sampleRate;
        fparams.order = 4;    //
        fparams.ncoeff = 0;                // init filter
        
        fparams.enable = true;
        if (build_filter (&fparams) == -1)
        {
            mexPrintf("Could not build filter.  Exiting\n");
            err[0] = -1;
            mxFree(dsName);
            mxFree(newDsName);
            return;
        }
    }
 
    // memory allocation
    
    // make buffer to read one trial of data at a time
    
    DE_trialArray = (double **)malloc( sizeof(double *) * dsParams.numChannels );
    if (DE_trialArray == NULL)
    {
        printf("memory allocation failed for DE_trialArray ");
        return;
    }
    for (int i = 0; i < dsParams.numChannels; i++)
    {
        DE_trialArray[i] = (double *)malloc( sizeof(double) * dsParams.numSamples);
        if ( DE_trialArray[i] == NULL)
        {
            printf( "memory allocation failed for DE_trialArray" );
            return;
        }
    }
        
    // make output array same number of channels (rows) as original
    // but set to number of samples in saved data to potentially save memory
        
    DE_trialArray2 = (double **)malloc( sizeof(double *) * dsParams.numChannels );
    if (DE_trialArray2 == NULL)
    {
        printf("memory allocation failed for DE_trialArray2 ");
        return;
    }
    for (int i = 0; i < dsParams.numChannels; i++)
    {
        DE_trialArray2[i] = (double *)malloc( sizeof(double) * newParams.numSamples);
        if ( DE_trialArray2[i] == NULL)
        {
            printf( "memory allocation failed for DE_trialArray2" );
            return;
        }
    }
        
    
    // allocate data arrays for one channel of data for filtering etc.
    
    channelData = (double *)malloc( sizeof(double) * dsParams.numSamples);
    if ( channelData == NULL)
    {
        mexPrintf("memory allocation failed for channel data buffer\n");
        err[0] = -1;
        mxFree(dsName);
        mxFree(newDsName);
        return;
    }
    
    channelData2 = (double *)malloc( sizeof(double) * dsParams.numSamples);
    if ( channelData2 == NULL)
    {
        mexPrintf("memory allocation failed for channel data 2 buffer\n");
        err[0] = -1;
        mxFree(dsName);
        mxFree(newDsName);
        return;
    }
    
    
    // create the ds directory and meg4 file with header
    if ( !createMEG4File( newDsName ) )
    {
        mexPrintf("Error creating new dataset \n");
        err[0] = -1;
        mxFree(dsName);
        mxFree(newDsName);
        return;
    }
    
    if (numBadTrials > 0)
        newParams.numTrials = dsParams.numTrials - numBadTrials;
        
    // adjust numChannels and counters once here rather than every trial loop
    for (int k=0; k<dsParams.numChannels; k++)
    {
        for (int j=0; j<numBadChannels; j++)
        {
            if (badChannelIndices[j] == k)
            {
                newParams.numChannels--;
                // adjust other counters
                if (dsParams.channel[k].isSensor)
                    newParams.numSensors--;
                if (dsParams.channel[k].isReference)
                    newParams.numReferences--;
                if (dsParams.channel[k].isBalancingRef)
                    newParams.numBalancingRefs--;
                sprintf(msg, "Excluding bad channel %d, (%s)\n", k, dsParams.channel[k].name);
                mexPrintf(msg);
                mexEvalString("drawnow");
            }
        }
    }
    
    sprintf(msg, "Writing new dataset: %s\n", newDsName);
    mexPrintf(msg);
    sprintf(msg,"Number of Channels %d\n", newParams.numChannels);
    mexPrintf(msg);
    sprintf(msg,"Number of Primary MEG Channels %d\n", newParams.numSensors);
    mexPrintf(msg);
    sprintf(msg,"Number of samples: %d\n", newParams.numSamples);
    mexPrintf(msg);
    sprintf(msg,"Sample Rate: %.1f\n", newParams.sampleRate);
    mexPrintf(msg);
    sprintf(msg,"Number of trials: %d\n", newParams.numTrials);
    mexPrintf(msg);
    sprintf(msg,"Trial duration: %.4f\n", newParams.trialDuration);
    mexPrintf(msg);
    sprintf(msg,"Epoch Min Time (s): %.4f\n", newParams.epochMinTime);
    mexPrintf(msg);
    sprintf(msg,"Epoch Max Time (s): %.4f\n", newParams.epochMaxTime);
    mexPrintf(msg);

    // write modified dataset
    
    
    for (int trial=0; trial <dsParams.numTrials; trial++)
    {

        bool includeTrial = true;
        if (numBadTrials > 0)
        {
           for (int k=0; k<numBadTrials; k++)
           {
               if (badTrialIndices[k] == trial)
               {
                   mexPrintf("*** Excluding trial %d ***\n", trial);
                   includeTrial = false;
               }
           }

        }
        else
            mexPrintf("Reading trial %d\n", trial);
            
        mexEvalString("drawnow");
        
        if (includeTrial)
        {
            // get one trial of data n channels x m samples - returns in 2D double array with gains applied
            if ( !readMEGTrialData( dsName, dsParams, DE_trialArray, trial, -1, 0) )
            {
                printf("Error reading .meg4 file\n");
                err[0] = -1;
                mxFree(dsName);
                mxFree(newDsName);
                return;
            }

            // now have all data in converted floating point values. Do post-processing
            // filter data (even for excluded channels...)
            if (filterData)
            {
                for (int k=0; k<dsParams.numChannels; k++)
                {
                    // if not digital replace data with filtered data
                    if (dsParams.channel[k].isSensor || dsParams.channel[k].isReference || dsParams.channel[k].isEEG )
                    {
                        for (int j=0; j<dsParams.numSamples; j++)
                            channelData[j] = DE_trialArray[k][j];
                        applyFilter( channelData, channelData2, dsParams.numSamples, &fparams);
                        for (int j=0; j< dsParams.numSamples; j++)
                            DE_trialArray[k][j] = channelData2[j];
                    }
                 }
            }
            
            // copy channel data and compress the number of rows in data array in place if dropping channels
            // may have less channels than original but extra rows will be ignored by writeMEGTrialData
            int newChannelCount = 0;
            for (int chan=0; chan<dsParams.numChannels; chan++)
            {
                bool includeChannel = true;
                if (numBadChannels > 0)
                {
                    for (int k=0; k<numBadChannels; k++)
                        if (badChannelIndices[k] == chan)
                            includeChannel = false;
                }
            
                if (includeChannel)
                {
                    // overwrite channel record
                    newParams.channel[newChannelCount] = dsParams.channel[chan];
                    
                    // copy modified data to second array - can be a subset of samples from original data array
                    int sampleCount = 0;
                    for (int k=startSample; k<endSample; k++)
                       DE_trialArray2[newChannelCount][sampleCount++] = DE_trialArray[chan][k];
                    newChannelCount++;  // increment counter or else this channel will get overwritten with next
                }
            }
            
            mexPrintf("Writing trial %d\n", trial);
            mexEvalString("drawnow");
            
            // append this trial to the meg4 file
            if ( !writeMEGTrialData( newDsName, newParams, DE_trialArray2) )
            {
                printf("Error occurred writing new .meg4 file\n");
                exit(0);
            }
        }

    } // next trial
    
    
    if ( !writeMEGResFile(newDsName, newParams) )
    {
        mexPrintf("WARNING: error occurred writing new res4 file... dataset may be invalid \n");
        err[0] = -1;
    }
    
            
    // copy over other necessary files
            
    //get filenames without path or ext..
    removeFilePath( dsName, dsBaseName);
    dsBaseName[strlen(dsBaseName)-3] = '\0';
    removeFilePath( newDsName, newDsBaseName);
    newDsBaseName[strlen(newDsBaseName)-3] = '\0';

    // copy head coil file
    #if _WIN32||WIN64
        sprintf(cmd, "copy %s%s%s.hc %s%s%s.hc", dsName, FILE_SEPARATOR, dsBaseName, newDsName, FILE_SEPARATOR, newDsBaseName );
    #else
        sprintf(cmd, "cp %s%s%s.hc %s%s%s.hc", dsName, FILE_SEPARATOR, dsBaseName, newDsName, FILE_SEPARATOR, newDsBaseName );
    #endif
    //printf("executing %s\n", cmd);
    errCode = system(cmd);

    // copy marker file
    #if _WIN32||WIN64
        sprintf(cmd, "copy %s%sMarkerFile.mrk %s%sMarkerFile.mrk", dsName, FILE_SEPARATOR, newDsName, FILE_SEPARATOR);
    #else
        sprintf(cmd, "cp %s%sMarkerFile.mrk %s%sMarkerFile.mrk", dsName, FILE_SEPARATOR, newDsName, FILE_SEPARATOR);
    #endif
    //printf("executing %s\n", cmd);
    errCode = system(cmd);
    
    mexPrintf("\n");
    mexEvalString("drawnow");
    mexPrintf("cleaning up...\n");
    
    // avoid memory leak
    for (int j = 0; j < dsParams.numChannels; j++)
        free(DE_trialArray[j]);
    free(DE_trialArray);
        
    mxFree(dsName);
    mxFree(newDsName);

    mexPrintf("done...\n");
    
    return;
}
}
