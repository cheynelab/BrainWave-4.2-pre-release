%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   getCTFHeadPosition
%
%   function [na le re] = getCTFHeadPosition(dsName, startSample, numSamples);
%
%   DESCRIPTION: Read CHL channels of a CTF dataset and return mean 
%                head position (fiducials in dewar coordinates) over the
%                specifed sample range
%
%   Version 4.0 updated Dec, 2023 with new mex function, no longer needs
%   mex function bw_CTFGetChannelLabels
%
% (c) D. Cheyne, 2014. All rights reserved.
% This software is for RESEARCH USE ONLY. Not approved for clinical use.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [na, le, re] = bw_getCTFHeadPosition(dsName, startSample, numSamples)

    CHL_channels = {'HLC0011'; 'HLC0012'; 'HLC0013'; 'HLC0021'; 'HLC0022'; 'HLC0023'; 'HLC0031'; 'HLC0032'; 'HLC0033'};
    
    na = [7 7 -23];
    le = [-7 7 -23];
    re = [7 -7 -23];
                        
    % Version 4 - simpler / faster version using new version of bw_getCTFData()
    header = bw_CTFGetHeader(dsName);
    longnames = {header.channel.name};
    channelNames = cleanChannelNames(longnames);
    idx = ismember(channelNames, CHL_channels);
    CHL_idx = find(idx == 1);

    data = bw_getCTFData(dsName, startSample, numSamples, 0, CHL_idx)';   

    % return mean x, y z position data for each fiducial in cm
    pos = mean(data,2) * 100;

    na = pos(1:3)';
    le = pos(4:6)';
    re = pos(7:9)';
end

% strip sensor version number from channel names for CTF
% make standalone function?
function channelNames = cleanChannelNames(names) 
    channelNames = [];
    if iscellstr(names)
        names = char(names);
    end
    for k=1:length(names) 
        s = names(k,:);
        ss = deblank(s);        % remove trailing whitespaces
        channelNames{k} = strtok(ss,'-');
    end
end


