%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function bw_dataEditor(dsName)
% GUI to view raw data - based on bw_eventMarker.m
%
%
% Input variables:
% dsName :         dataset path
% 
% (c) D. Cheyne. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bw_dataEditor(dsName)

versionNo = 4.2;
warning('off', 'MATLAB:linkaxes:RequireDataAxes');
    
tStr = sprintf('Data Editor(ver %.1f)', versionNo);

fprintf(tStr);
fprintf('\n(c) D. Cheyne (2023) Hospital for Sick Children. Version %.2f\n', versionNo);

timeVec = [];
dataarray = [];
markerFileName = [];

numMarkers = 0;
currentMarkerIndex = 1;
markerNames = {};
markerLatencies = {};
markerTrials = {};
currentMarker = 1;
enableMarking = 0;
amplitudeLabels = [];
amplitudeUnits = {};

channelName = [];
channelMenuItems = [];

eventList = [];  
currentEvent = 0;
numEvents = 0;

threshold = 0;
maxAmplitude = 0.0;
minAmplitude = 0.0;
minSeparation = 0.0;
showAmplitudes = 0;

currentScaleMenuIndex = 1;
% set all data ranges to autoscale first time
maxRange = [NaN NaN NaN NaN NaN NaN];
minRange = [NaN NaN NaN NaN NaN NaN];

% channel types according to CTF sensorType
MEG_CHANNELS = [0 1 2 3 4 5];
ADC_CHANNELS = 18;
TRIGGER_CHANNELS = 19;
DIGITAL_CHANNELS = 20;
EEG_CHANNELS = [9 21];  % 9 = unipolar, 21 = bipolar
OTHER_CHANNELS = [6 7 8 10 11 12 13 14 15 16 17 22 23 24 25 26 27 28 29];    

header = [];
channelNames = [];
channelTypes = [];
selectedChannelList = [];
channelMenuIndex = 1;
customChannelList = 1;
selectedMask = 1;
badChannelMask = 0;
badTrialMask = 0;
% set defaults
rectify = false;
envelope = false;
notchFilter = false;
notchFilter2 = false;
invertData = false;
differentiate = false;
removeOffset = false;
epochStart = 0.0;
epochTime = 0.0;
trialNo = 1;


bandPass = [1 50];
minDuration = 0.00;
filterOff = true;
minScale = NaN;
maxScale = NaN;
epochSamples = 0;
enableMarking = 0;
reverseScan = false;

cursorHandles = [];
plotAxes = [];

cursorLatency = 0.0;
overlayPlots = 0;
numColumns = 1;


showMap = 0;

% default channel sets.
channelSets = {'<HTML><FONT COLOR="black">Custom</HTML>'; ...
        '<HTML><FONT COLOR="blue">MEG Sensors</HTML>'; ...
        '<HTML><FONT COLOR="rgb(26,179,26)">ADC Channels</HTML>';...
        '<HTML><FONT COLOR="rgb(204,128,26)">Trigger Channel</HTML>';...
        '<HTML><FONT COLOR="black">Digital Channels</HTML>';...
        '<HTML><FONT COLOR="cyan">EEG/EMG Channels</HTML>'};
darkGreen = [0.1 0.7 0.1];
orange = [0.8 0.5 0.1];
gray = [0.5 0.5 0.5];

% set defaults
% Draw arrows - calls:  uparrow.m and downarrow.m - %%%% ADDED BY CECILIA %%%%
uparrow_im=draw_uparrow;
downarrow_im=draw_downarrow;
rightarrow_im=draw_rightarrow;
leftarrow_im=draw_leftarrow;

defaultsFile = [];
dsPath = [];

% in case launched from command line
tpath = which('bw_dataEditor');
pathparts = strsplit(tpath,filesep);
s = pathparts(1:end-1);
DE_PATH = strjoin(s,filesep);
externalPath = strcat(DE_PATH,filesep,'external');
dirpath=strcat(externalPath,filesep,'topoplot');
if exist(dirpath,'dir') ~= 7   % should not happen as folder is part of BW
    fprintf('error: topoplot folder is missing...\n');
else
    addpath(dirpath);
end

% attempt to load all data - needed for faster plotting and topoplot
all_data = [];
meg_idx = [];

sx = 200;
sy = 800;
swidth = 1600;
sheight = 1200;

left2 = sx + swidth;
top = sy + sheight;

fh = figure('numbertitle','off','position',[sx, sy, swidth, sheight],...
    'Name','Data Editor', 'Color','white','menubar','none','WindowButtonUpFcn',@stopdrag,'WindowButtonDownFcn',@buttondown, 'keypressfcn',@capturekeystroke);

if ispc
    movegui(fh, 'center');
end

filemenu=uimenu('label','File');
uimenu(filemenu,'label','Open CTF Dataset','accelerator','O','callback',@openFile_callback)
importMenu = uimenu(filemenu,'label','Import MEG Data');
uimenu(importMenu,'label','Neuromag/MEGIN data (*.fif)...','callback',@import_fif_data_callback)
uimenu(importMenu,'label','KIT data (*.con)...','callback',@import_kit_data_callback)
uimenu(filemenu,'label','Save Dataset As...','separator','on','callback',@saveFile_callback)

uimenu(filemenu,'label','Open layout','accelerator','L','separator','on','callback',@open_layout_callback)
uimenu(filemenu,'label','Save layout','accelerator','S','callback',@save_layout_callback)

uimenu(filemenu,'label','Close','accelerator','W','separator','on','callback',@quit_filemenu_callback)

channelMenu=uimenu('label','ChannelSets');

markerMenu=uimenu('label','Edit Events');

loadMenu = uimenu(markerMenu,'label','Load Events');
uimenu(loadMenu,'label','from MarkerFile','callback',@load_marker_events_callback)
uimenu(loadMenu,'label','Import from Text File','callback',@load_events_callback)
uimenu(loadMenu,'label','Import from KIT Event File','callback',@load_KIT_events_callback)


saveMenu = uimenu(markerMenu,'label','Save Events');
uimenu(saveMenu,'label','to MarkerFile','callback',@save_marker_callback)
uimenu(saveMenu,'label','Export to Text File','callback',@save_events_callback)
uimenu(markerMenu,'label','Create Conditional Event...','callback',@create_event_callback)
uimenu(markerMenu,'label','Edit Markers...','separator','on','callback',@edit_markers_callback)
uimenu(markerMenu,'label','Export Markers to Excel...','callback',@save_marker_as_excel_callback)

set(markerMenu,'enable','off');

% build channel menu once
for kk=1:numel(channelSets)
    s = char(channelSets(kk));
    channelMenuItems(kk) = uimenu(channelMenu,'label',s,'callback',@channel_menu_callback); 
end
channelMenuItems(end+1) = uimenu(channelMenu,'label','Edit Custom',...
        'separator','on','callback',@editChannelSet_callback); 


% +++++++++++++ set plot windows +++++++++++++


mapbox = [0.46 0.005 0.2 0.2];
mapLocs = [];

plotbox = [0.05 0.31 0.9 0.65];
subplot('position',plotbox);
    

function openFile_callback(~,~)
    dsName = uigetdir('.ds', 'Select CTF dataset ...');
    if dsName == 0
        return;
    end

    initData;
    drawTrial;

end

function saveFile_callback(~,~)
    
    applyFilter = 0;
    excludeChans = 0;

    if ~filterOff
        r = questdlg('Apply current filter setting to saved data?,','Data Editor','Yes','No','Cancel','No');
        if strcmp(r,'Cancel')
            return;
        end 
        if strcmp(r,'Yes')
            applyFilter = 1;
        end
    end
    
    % check bad channels
    badChanIdx = find(badChannelMask == 1);
    badTrialIdx = find(badTrialMask == 1);
    
    if ~isempty(badChanIdx) || ~isempty(badTrialIdx)
        r = questdlg('Exclude bad channels and bad trials?','Data Editor','Yes','No','Cancel','No');
        if strcmp(r,'Cancel')
            return;
        end 
        if strcmp(r,'Yes')
            excludeChans = 1;
        end
    end
    
            
    s1 = num2str(header.epochMinTime);
    s2 = num2str(header.epochMaxTime);
    s3 = '1';
    s4 = num2str(header.gradientOrder);
    defaultanswer={s1,s2,s3,s4};
    answer=inputdlg({'Start time (s):','End Time (s)','Downsample Factor (1 = no downsampling):',...
        'Gradient Order: (0=raw, 1=1st, 2=2nd, 3=3rd, 4=3rd+adaptive)'},'Data Parameters',[1 100; 1 100; 1 100; 1 100],defaultanswer);

    if isempty(answer)
        return;
    end
 
    t1 = str2num(answer{1});
    t2 = str2num(answer{2});
    tt = str2num(answer{3});
    ds = round(tt);
    tt = str2num(answer{4});
    gradient = round(tt);
    
    % get sample range - base 0 for mex function
    startSample = round( (t1 - header.epochMinTime) * header.sampleRate );
    endSample  = round( (t2 - header.epochMinTime) * header.sampleRate );  
    
    if startSample < 0 || endSample > header.numSamples || startSample > endSample
        errordlg('Invalid time range entered');
        return;
    end  
    
    newFs = header.sampleRate;
    if ds > 1
        newFs = header.sampleRate / ds;
        if newFs < bandPass(2) * 2.0
            s = sprintf('New sample rate (%.1f Samples/s) is too low for current lowpass filter (%.1f Hz)', newFs, bandPass(2));
            errordlg(s);
            return;         
        end
    end
    
    if gradient ~= header.gradientOrder
        if ~header.hasBalancingCoefs
            errordlg('No gradient reference channels for this dataset');
            return;
        end
        if gradient < 0 || gradient > 4
            s = sprintf('Invalid Gradient (%d) (0=raw, 1=1st, 2=2nd, 3=3rd, 4=3rd+adaptive)', gradient); 
            errordlg(s);
            return;
        end
    end
    
    sampleRange = [startSample endSample];
        
    % auto modify save name to avoid duplications...
    appendStr = '_copy';
    if filterOff == 0
        appendStr = sprintf('_f%d_%dHz', round(bandPass));
    end
    
    if ~isempty(badChanIdx) || ~isempty(badTrialIdx)
        appendStr = strcat(appendStr, '_e');
    end
    
    if gradient ~= header.gradientOrder
        s = sprintf('_g%d',gradient);
        appendStr = strcat(appendStr, s);     
    end

    [~,n,~] = fileparts(dsName);
  
    name = sprintf('%s%s.ds',n,appendStr);
    if exist(name,'dir')
        name = sprintf('%s%s%s.ds',n,appendStr,appendStr);
    end
    [n,p] = uiputfile('*.*','Save as:', name);
    if isequal(n,0)
        return;
    end         
 
    newDsName = fullfile(p,n);
    fprintf('Saving data in: %s \n', newDsName);

    wbh = waitbar(0,'Saving data...');

    % sanity check for mex function    
    if filterOff || applyFilter == 0
        filterFlag = 0;
    else
        filterFlag = 1;
    end       
    
    badChans = [];
    badTrials = [];
    if excludeChans
        if ~isempty(badChanIdx) 
            % for mex function channel numbers start at zero.
            badChans = badChanIdx-1;
        end

        if ~isempty(badTrialIdx)
            % for mex function channel numbers start at zero.
            badTrials = badTrialIdx-1;
        end
    end   
       
    waitbar(0.5,wbh,'Writing data...');

    err = bw_CTFNewDs(dsName, newDsName, filterFlag, bandPass, badChans, badTrials, sampleRange, ds, gradient);  
    if err ~= 0
        errordlg('bw_CTFNewDs returned error');
        return;
    end

    waitbar(1.0,wbh,'Saving markers...');
   
    % if time zero has been shifted > 0.0 seconds have to correct marker
    % latencies. If dropped trials have to correct trial numbers !

    markerFileName = sprintf('%s%smarkerFile.mrk',dsName,filesep);
    hasMarkers = exist(markerFileName,'file');
    
    % if has markerFile and needs correcting....
    if hasMarkers && ( t1 > 0.0 || ~isempty(badTrialIdx) )
        
        % note bw_readCTFMarkerFile adds 1 to the trial numbers ...
        [markerNames, markerData] = bw_readCTFMarkerFile(markerFileName);     
        
        % for each marker correct trial number and/or latency
        for k=1:numel(markerData)            
            newMarkerData(k).ch_name = char(markerNames{k}); 
           
            markerTimes = markerData{k};
            trials = markerTimes(:,1);      % numbering starts at 1
            latencies = markerTimes(:,2);
                       
            % remove deleted trials and renumber
            if ~isempty(badTrialIdx)     
                fprintf('correcting marker trial numbers...');
                trials(badTrialIdx) = [];       % compress lists to valid trials only
                latencies(badTrialIdx) = [];
                
                % renumber the trials - tricky!               
                for j=1:length(trials)
                    n = trials(j);
                    idx = find(badTrialIdx < n);   % which trials preceding this one deleted?
                    shift = length(idx);           % how many? shift = 0 if idx = []
                    trials(j) = n - shift;                
                end
            end
            
            % correct latencies if time zero was shifted forward
            if t1 > 0.0 
                fprintf('correcting marker latencies (subtracting %.4f seconds)\n', t1);
                latencies = latencies - t1;
            end
            
            newMarkerData(k).trials = trials - 1;      % correct back to base zero before writing
            newMarkerData(k).latencies = latencies;          
        end
        
        % write corrected markerFile to the new dataset (overwrite
        bw_writeNewMarkerFile(dsName, newMarkerData);     
        
     end

    delete(wbh);

    dsName = newDsName;
    initData;
    drawTrial;
     
end

function quit_filemenu_callback(~,~)
      
    response = questdlg('Quit Data Editor?','BrainWave','Yes','No','No');
    if strcmp(response,'No')    
        return;
    end       

    response = questdlg('Save current layout?','BrainWave','Yes','No','No');
    if strcmp(response,'Yes')    
        saveDefaults(defaultsFile);
    end       

    close(fh);
end

function save_layout_callback(~,~)
   
    s = fullfile(dsName,'myLayout');
    [name,path,~] = uiputfile('*.mat','Save current layout in File:',s);
    if isequal(name,0)
        return;
    end         
    
    saveFile = fullfile(path,name);
    fprintf('Saving layout in: %s \n', saveFile);
    saveDefaults(saveFile);

end

    
function open_layout_callback(~,~)
   
    [name, path, ~] = uigetfile( ...
        {'*.mat','Layout File (*.mat)'}, ...
           'Select Layout File', defaultsFile);
    if isequal(name,0) 
       return;
    end    
    layoutFile =fullfile(path,name);

    readDefaults(layoutFile);
    drawTrial;

end

function readDefaults(defaultsFile)

    if ~exist(defaultsFile,'file')
        return;
    end

    fprintf('Loading previous layout from %s ...\n', defaultsFile)
    params = load(defaultsFile);
    channelMenuIndex = params.channelMenuIndex;
    selectedChannelList = params.selectedChannelList;
    selectedMask = zeros(1,numel(selectedChannelList));
    badChannelMask = zeros(1,header.numChannels);
    badTrialMask = zeros(1,header.numTrials);

    overlayPlots = params.overlayPlots;
    if isfield(params,'epochTime')
        epochTime = params.epochTime; 
        s = sprintf('%.4f', epochTime);
        set(epochDurationEdit,'string', s);
    end
    
    if isfield(params,'numColumns')
        numColumns = params.numColumns;
        set(numColumnsMenu,'value',numColumns); 
    end
    
    set(overlayPlotsCheck,'value',overlayPlots);

    updateChannelMenu;   
    
end

function saveDefaults(defaultsFile)
    params.channelMenuIndex = channelMenuIndex;
    params.selectedChannelList = selectedChannelList;
    params.overlayPlots = overlayPlots;
    params.epochTime = epochTime;
    params.numColumns = numColumns;
    
    save(defaultsFile, '-struct', 'params');
end

function initData
        
    if ~exist(dsName,'file')
        fprintf('Could not find file %s...\n', dsName);
        return;      
    end

    [dsPath,~,~] = fileparts(dsName);
    
    if ~isempty(dsPath)
        cd(dsPath);
    end
    
    defaultsFile = sprintf('%s%sdataEditor.mat', dsName, filesep);

    eventList = [];  
    currentEvent = 0;
    numEvents = 0;
    set(eventCtl(:),'enable','off');
    s = sprintf('Events: %d of %d', currentEvent, numEvents);
    set(numEventsTxt,'String',s);
    
    enableMarking = 0;

    threshold = 0;
    maxAmplitude = 0.0;
    minAmplitude = 0.0;
    minSeparation = 0.0;
    
    header = bw_CTFGetHeader(dsName);

    meg_idx = [];  
    all_data = [];  % force re-read of raw data (all channels);

    % clear map plot
    showMap = 0;
    mapLocs = [];
    set(showMapCheck,'value',0);
   
    currentAx = gca;
    subplot('Position',mapbox);
    cla;
    axis off
    axes(currentAx);
%     subplot('Position',plotbox);    

    [~,n,e] = fileparts(dsName);
    tStr = sprintf('Data Editor: %s', [n e]);
    set(fh,'Name', tStr);

    timeVec = [];
    dataarray = [];

    numMarkers = 0;
    currentMarkerIndex = 1;
    currentMarker = 1;
    set(marker_Popup,'value',1)
    markerLatencies = {};
    markerTrials = {};
    set(markerIncButton,'enable','off')
    set(markerDecButton,'enable','off')
    set(markerStartButton,'enable','off')
    set(markerEndButton,'enable','off')

    % forces autoscaling 
    maxRange = [NaN NaN NaN NaN NaN NaN];
    minRange = [NaN NaN NaN NaN NaN NaN];
    
    epochStart = header.epochMinTime;
    epochTime = 10.0;
    if epochTime > header.epochMaxTime
        epochTime = header.epochMaxTime;
    end
    s = sprintf('%.4f', epochTime);
    set(epochDurationEdit,'string', s);

    trialNo = 1;
    if header.numTrials > 1
        set(trialIncButton,'enable','on');
        set(trialDecButton,'enable','on');                    
        set(trialStartButton,'enable','on');
        set(trialEndButton,'enable','on');
        
    else
        set(trialIncButton,'enable','off');
        set(trialDecButton,'enable','off');
        set(trialStartButton,'enable','off');
        set(trialEndButton,'enable','off');
    end
    s = sprintf('Trial: %d of %d', trialNo, header.numTrials);
    set(trialNumTxt,'string',s);

    cursorLatency = header.epochMinTime + (epochTime/2);
    epochSamples = round(epochTime * header.sampleRate);

    % when loading new dataset reset bandpass and turn filter off

    set(filt_hi_pass,'string',bandPass(1),'enable','off')  
    set(filt_low_pass,'string',bandPass(2),'enable','off')
    filterOff = 1;
    set(filterCheckBox,'value',0);
    
    fprintf('Loading data...\n\n');
    
    wbh = waitbar(0,'Loading data...');
            
    
    markerFileName = strcat(dsName,filesep,'MarkerFile.mrk');
    
    markerNames = {'none'};
    if exist(markerFileName,'file')
        loadMarkerFile(markerFileName);
    else

       fprintf('no marker file found...\n'); 
       set(marker_Popup,'string','None')
       numMarkers = 0;
    end
   
    % ** need to rebuild channel set menu with valid channel types *** 

    s = sprintf('Sample Rate: %4.1f S/s',header.sampleRate);
    set(sampleRateTxt,'string',s);
 
    s = sprintf('Total Samples: %d',header.numSamples);
    set(totalSamplesTxt,'string',s);
   
    s = sprintf('# of Trials: %d',header.numTrials);
    set(numTrialsTxt,'string',s);    
    
    s = sprintf('Trial Duration: %.4f s',header.trialDuration);
    set(totalDurationTxt,'string',s);    
    
    
    s = sprintf('Total Channels: %d',header.numChannels);
    set(numChannelsTxt,'string',s);
    
    s = sprintf('MEG Sensors: %d',header.numSensors);
    set(numSensorsTxt,'string',s);
   
    s = sprintf('MEG References: %d',header.numReferences);
    set(numReferencesTxt,'string',s);
    
    longnames = {header.channel.name};   
    channelNames = cleanChannelNames(longnames);
    channelTypes = [header.channel.sensorType];
    idx = ismember(channelTypes,ADC_CHANNELS);
    numAnalog = length(find(idx == 1));
    s = sprintf('ADC Channels: %d',numAnalog);
    set(numAnalogTxt,'string',s);
    
    % set default display on opening to first MEG sensor...
    meg_idx = find(channelTypes == 5);
    if ~isempty(meg_idx)
        selectedChannelList = meg_idx(1);
    else
        selectedChannelList = 1;
    end
    
    customChannelList = selectedChannelList;
    channelMenuIndex = 1;

    selectedMask = zeros(1,numel(selectedChannelList));
    badChannelMask = zeros(1,header.numChannels);
    
    badTrialMask = zeros(1,header.numTrials);
    if header.numTrials < 2
        set(setBadTrialButton, 'enable','off');
    else
        set(setBadTrialButton, 'enable','on');
    end        
    
    % override default settings
    readDefaults(defaultsFile);  

    waitbar(0.5,wbh,'Loading data...');

    
    loadData;   
    
    waitbar(1.0,wbh,'done...');
  
    delete(wbh);
    
    drawTrial;
    updateMap;

    updateSlider;
    updateCursors;

    updateChannelMenu;

    set(markerMenu,'enable','on');
    updateMarkerControls;

    maxScale = maxRange(currentScaleMenuIndex);
    minScale = minRange(currentScaleMenuIndex);
    s = sprintf('%.3g', maxScale);
    set(max_scale,'string',s);
    s = sprintf('%.3g', minScale);
    set(min_scale,'string',s);


end

function loadMarkerFile(markerFileName)
        fprintf('reading marker file %s\n', markerFileName); 
        [names, markerData] = bw_readCTFMarkerFile(markerFileName);       
        % drop trial numbers for now
        numMarkers = size(names,1);
        fprintf('dataset has %d markers ...\n', numMarkers); 
        markerNames = {'none'};
        if (numMarkers > 0)
            for j = 1:numMarkers
                x = markerData{j}; 
                markerLatencies{j} = x(:,2);
                markerTrials{j} = x(:,1);
                markerNames{j+1} = names{j};
            end
            markerNames{numMarkers+2} = 'All Markers';
        end
        set(marker_Popup,'string',markerNames)
end

function updateChannelMenu
        
    menuItems = get(channelMenu,'Children');
    set(menuItems(:),'Checked','off');

    % Note: menuItems indices are always in reverse order to their order in the menu !
    menuIdx = numel(menuItems):-1:1; 
    idx = menuIdx(channelMenuIndex);
    set(menuItems(idx),'Checked','on');

    for k=2:numel(channelSets)
        % turn off default channel types that don't exist
        idx = menuIdx(k);  % flipped index order for menu
        switch k 
            case 2
                if ~ismember(channelTypes,MEG_CHANNELS)
                    set(menuItems(idx),'enable','off');
                end
            case 3
                if ~ismember(channelTypes,ADC_CHANNELS)
                    set(menuItems(idx),'enable','off');
                end    
            case 4
                if ~ismember(channelTypes,TRIGGER_CHANNELS)
                    set(menuItems(idx),'enable','off');
                end                     
            case 5
                if ~ismember(channelTypes,DIGITAL_CHANNELS)
                    set(menuItems(idx),'enable','off');
                end         
            case 6
                if ~ismember(channelTypes,EEG_CHANNELS)
                    set(menuItems(idx),'enable','off');
                end
        end 
    end

end

function channel_menu_callback(src,~)

    if isempty(channelNames)
        return;
    end
    
    selection = get(src,'position');
    nchans = numel(channelNames);
    channelExcludeFlags = ones(1,nchans);             

    switch selection 
        case 1
            channelExcludeFlags(customChannelList) = 0;
        case 2
            idx = find(ismember(channelTypes,5));
            channelExcludeFlags(idx) = 0;
        case 3
            idx = find(ismember(channelTypes,ADC_CHANNELS));
            channelExcludeFlags(idx) = 0;
        case 4
            idx = find(ismember(channelTypes,TRIGGER_CHANNELS));
            channelExcludeFlags(idx) = 0;
        case 5
            idx = find(ismember(channelTypes,DIGITAL_CHANNELS));
            channelExcludeFlags(idx) = 0;              
        case 6
            idx = find(ismember(channelTypes,EEG_CHANNELS));
            channelExcludeFlags(idx) = 0;
    end        

    selectedChannelList = find(channelExcludeFlags == 0);
    selectedMask = zeros(1,numel(selectedChannelList));
        
    menuItems = get(channelMenu,'Children');
    set(menuItems(:),'Checked','off');
    if selection < numel(menuItems)      % don't check last menu item
        channelMenuIndex = selection;
        set(src,'Checked','on')
    end

    processData;
    
    updateMarkerControls; 

    autoScale_callback;  % calls drawTrial;

      
end

function editChannelSet_callback(~,~)  
       
    % get new custom set
    [selected] = bw_channelSelector(header, selectedChannelList, badChannelMask);
    if isempty(selected)
        return;
    end
    
    customChannelList = selected;
    selectedChannelList = customChannelList;
    selectedMask = zeros(1,numel(selectedChannelList));

    channelMenuIndex = 1;
    % get handles to menu items - indices always returned in reverse order
    % to which they appear in the menu ...
    channelMenuItems = get(channelMenu,'Children');
    set(channelMenuItems(:),'Checked', 'off')
    set(channelMenuItems(end), 'Checked','on')  

    processData;
    
    drawTrial;
    autoScale_callback;
    updateMarkerControls;

end

%+++++++++ event detection controls +++++++++  

annotation('rectangle',[0.75 0.02 0.2 0.18],'EdgeColor','blue');

% uicontrol('style','text','units','normalized','position',[0.81 0.185 0.12 0.025],...
%     'string','Mark Events','backgroundcolor','white','foregroundcolor','blue','fontweight','bold',...
%     'FontSize',11);
   
enable_marking_check = uicontrol('style','checkbox','units','normalized','fontsize',11,'position',[0.79 0.188 0.12 0.025],...
    'enable','off','string','Enable Threshold Marking','Foregroundcolor','blue','backgroundcolor','white','callback',@enable_marking_callback);

function updateMarkerControls

    if numel(selectedChannelList) == 1 && header.numTrials == 1
        set(enable_marking_check,'enable','on');
    else
        set(enable_marking_check,'enable','off');
        set(enable_marking_check,'value',0);
        enableMarking = 0; 
        setMarkingCtls('off');
    end
end

function setMarkingCtls(state)
    set(threshold_text,'enable',state)
    set(threshold_edit,'enable',state)
    set(min_duration_text,'enable',state)
    set(min_duration_edit,'enable',state)
    set(min_amplitude_edit,'enable',state)
    set(max_amplitude_edit,'enable',state)
    set(min_amp_text,'enable',state)
    set(min_amp_text2,'enable',state)
    set(min_sep_text,'enable',state)
    set(min_sep_edit,'enable',state)
    set(min_sep_text,'enable',state)
    set(forward_scan_radio,'enable',state)
    set(reverse_scan_radio,'enable',state)
    set(find_events_button,'enable',state)    
end

find_events_button = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.9 0.04 0.04 0.04],...
    'enable','off','string','Scan','Foregroundcolor','blue','backgroundcolor','white','callback',@find_events_callback);

forward_scan_radio = uicontrol('style','radiobutton','units','normalized','position',[0.77 0.165 0.06 0.02],...
    'enable','off','string','rising edge','backgroundcolor','white','value',~reverseScan,'FontSize',11,'callback',@forward_check_callback);
reverse_scan_radio = uicontrol('style','radiobutton','units','normalized','position',[0.84 0.165 0.06 0.02],...
    'enable','off','string','falling edge','backgroundcolor','white','value',reverseScan,'FontSize',11,'callback',@reverse_check_callback);

threshold_text = uicontrol('style','text','units','normalized','position',[0.77 0.135 0.05 0.02],...
    'enable','off','string','Threshold:','fontsize',11,'backgroundcolor','white','horizontalalignment','left');
threshold_edit=uicontrol('style','edit','units','normalized','position',[0.85 0.14 0.035 0.02],...
    'enable','off','FontSize', 11, 'BackGroundColor','white','string',threshold,'callback',@threshold_callback);
    function threshold_callback(src,~)
        threshold =str2double(get(src,'string'));
        drawTrial;
    end

min_amp_text = uicontrol('style','text','units','normalized','position',[0.77 0.1 0.05 0.03],...
    'enable','off','string','Amp. Range:','fontsize',11,'backgroundcolor','white','horizontalalignment','left');
min_amplitude_edit=uicontrol('style','edit','units','normalized','position',[0.85 0.105 0.035 0.02],...
    'enable','off','FontSize', 11, 'BackGroundColor','white','string',minAmplitude,...
    'callback',@min_amplitude_callback);
    function min_amplitude_callback(src,~)
        minAmplitude=str2double(get(src,'string'));
        drawTrial;
    end
max_amplitude_edit=uicontrol('style','edit','units','normalized','position',[0.9 0.105 0.035 0.02],...
    'enable','off','FontSize', 11, 'BackGroundColor','white','string',maxAmplitude,...
    'callback',@max_amplitude_callback);
min_amp_text2 = uicontrol('style','text','units','normalized','position',[0.89 0.1 0.01 0.02],...
    'enable','off','string','to','fontsize',11,'backgroundcolor','white','horizontalalignment','left');
    function max_amplitude_callback(src,~)
        maxAmplitude=str2double(get(src,'string'));
        drawTrial;
    end

min_duration_text = uicontrol('style','text','units','normalized','position',[0.77 0.065 0.08 0.02],...
    'enable','off','string','Min. duration (s):','fontsize',11,'backgroundcolor','white','horizontalalignment','left');
min_duration_edit=uicontrol('style','edit','units','normalized','position',[0.85 0.07 0.035 0.02],...
    'enable','off','FontSize', 11, 'BackGroundColor','white','string',minDuration,...
    'callback',@min_duration_callback);

    function min_duration_callback(src,~)
        minDuration=str2double(get(src,'string'));
    end

min_sep_text = uicontrol('style','text','units','normalized','position',[0.77 0.03 0.08 0.02],...
    'enable','off','string','Min. separation (s):','fontsize',11,'backgroundcolor','white','horizontalalignment','left');
min_sep_edit = uicontrol('style','edit','units','normalized','position',[0.85 0.035 0.035 0.02],...
    'enable','off','FontSize', 11, 'BackGroundColor','white','string',minSeparation,...
    'callback',@min_separation_callback);

    function min_separation_callback(src,~)
        minSeparation=str2double(get(src,'string'));
        drawTrial;
    end


s = sprintf('Events: %d of %d', currentEvent, numEvents);
numEventsTxt = uicontrol('style','text','units','normalized','position',[0.63 0.963 0.1 0.025],...
    'string',s,'fontsize',12,'backgroundcolor','white', 'foregroundcolor', 'black','horizontalalignment','left');
eventCtl(1) = uicontrol('style','pushbutton','units','normalized','fontsize',14,'fontweight','bold','position',[0.84 0.97 0.02 0.025],...
    'enable','off','string','<<','Foregroundcolor','black','backgroundcolor','white','callback',@event_first_callback);
eventCtl(2) = uicontrol('style','pushbutton','units','normalized','fontsize',14,'fontweight','bold','position',[0.865 0.97 0.02 0.025],...
    'enable','off','String','<','Foregroundcolor','black','backgroundcolor','white','callback',@event_dec_callback);
eventCtl(3) = uicontrol('style','pushbutton','units','normalized','fontsize',14,'fontweight','bold','position',[0.89 0.97 0.02 0.025],...
    'enable','off', 'String', '>','Foregroundcolor','black','backgroundcolor','white','callback',@event_inc_callback);
eventCtl(4) = uicontrol('style','pushbutton','units','normalized','fontsize',14,'fontweight','bold','position',[0.915 0.97 0.02 0.025],...
    'enable','off', 'string','>>','Foregroundcolor','black','backgroundcolor','white','callback',@event_last_callback);
eventCtl(5) = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.69 0.97 0.04 0.025],...
    'enable','off','string','Insert','Foregroundcolor','black','backgroundcolor','white','callback',@add_event_callback);
eventCtl(6) = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.735 0.97 0.04 0.025],...
    'enable','off','string','Delete','Foregroundcolor','black','backgroundcolor','white','callback',@delete_event_callback);
eventCtl(7) = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.78 0.97 0.05 0.025],...
    'enable','off','string','Delete All','Foregroundcolor','black','backgroundcolor','white','callback',@delete_all_callback);


    function event_inc_callback(~,~)  
        if numEvents < 1
            return;
        end

        if currentEvent < numEvents
            currentEvent = currentEvent + 1;
            epochStart = eventList(currentEvent) - (epochTime / 2);
            cursorLatency = eventList(currentEvent);
            
            drawTrial;     
            updateCursors;
            % adjust slider position
            val = (epochStart + header.epochMinTime) / header.trialDuration;
            if val < 0, val = 0.0; end
            if val > 1.0, val = 1.0; end
            set(latency_slider, 'value', val);
            s = sprintf('Events: %d of %d', currentEvent, numEvents);
            set(numEventsTxt,'string',s);

        end
    end
    function event_dec_callback(~,~)  
        if numEvents < 1
            return;
        end

        if currentEvent > 1
            currentEvent = currentEvent - 1;
            epochStart = eventList(currentEvent) - (epochTime / 2);
            cursorLatency = eventList(currentEvent);

            drawTrial;
            updateCursors;
            % adjust slider position
            val = (epochStart + header.epochMinTime) / header.trialDuration;
            if val < 0, val = 0.0; end
            if val > 1.0, val = 1.0; end
            set(latency_slider, 'value', val);
            s = sprintf('Events: %d of %d', currentEvent, numEvents);
            set(numEventsTxt,'string',s);
        end
    end

    function event_first_callback(~,~)  
        if numEvents < 1
            return;
        end

        currentEvent = 1;
        epochStart = eventList(currentEvent) - (epochTime / 2);
        cursorLatency = eventList(currentEvent);
        
        drawTrial;     
        updateCursors;
        % adjust slider position
        val = (epochStart + header.epochMinTime) / header.trialDuration;
        if val < 0, val = 0.0; end
        if val > 1.0, val = 1.0; end
        set(latency_slider, 'value', val);
        s = sprintf('Events: %d of %d', currentEvent, numEvents);
        set(numEventsTxt,'string',s);
    end

    function event_last_callback(~,~)  
        if numEvents < 1
            return;
        end

        currentEvent = numEvents;
        epochStart = eventList(currentEvent) - (epochTime / 2);
        cursorLatency = eventList(currentEvent);
        
        drawTrial;     
        updateCursors;
        % adjust slider position
        val = (epochStart + header.epochMinTime) / header.trialDuration;
        if val < 0, val = 0.0; end
        if val > 1.0, val = 1.0; end
        set(latency_slider, 'value', val);
        s = sprintf('Events: %d of %d', currentEvent, numEvents);
        set(numEventsTxt,'string',s);
    end



    function add_event_callback(~,~)
        % manually add event

        response = questdlg('Add new event at cursor latency?','BrainWave','Yes','No','Yes');
        if strcmp(response,'No')
            return;
        end

        latency = cursorLatency;
        if isempty(eventList)
            eventList = latency;
            currentEvent = 1;
        else
            % add to list and resort - idx tells me where the last value
            % moved to in the sorted list
            eventList = [eventList latency];
            [eventList, idx] = sort(eventList);
            currentEvent = find(idx == length(eventList));
        end
        numEvents = length(eventList);
        s = sprintf('Events: %d of %d', currentEvent, numEvents);
        set(numEventsTxt,'String',s);

        drawTrial;
    end

    function delete_event_callback(~,~)

       if numEvents < 1
           errordlg('No events to delete ...');
           return;
       end
       % make sure we have currentEvent in window
       if ( eventList(currentEvent) >= epochStart && eventList(currentEvent) < (epochStart+epochTime) )        
           s = sprintf('Delete event #%d? (Cannot be undone)', currentEvent);
           response = questdlg(s,'Mark Events','Yes','No','No');
           if strcmp(response,'No')
               return;
           end
           eventList(currentEvent) = [];
           if currentEvent > 1
                currentEvent = currentEvent - 1;   
           end
           numEvents = length(eventList);  
           drawTrial;
           s = sprintf('Events: %d of %d', currentEvent, numEvents);
           set(numEventsTxt,'String',s);
       end
    end


    function delete_all_callback(~,~)
        if numEvents < 1
           fprintf('No events to delete ...\n');
           return;
        end

        response = questdlg('Clear all events?','Event Marker','Yes','No','No');
        if strcmp(response','No')
           return;
        end

        eventList = [];
        numEvents = 0;
        currentEvent = 0;             

        drawTrial;
        s = sprintf('Events: %d of %d', currentEvent, numEvents);
        set(numEventsTxt,'string',s)
        
    end

    % now works on normalized scale
    function enable_marking_callback(src,~)
        enableMarking = get(src,'value');
      
        if enableMarking
            % reset params for marking        
            setMarkingCtls('on');

            maxAmplitude = 0.8;
            minAmplitude = 0.0;
            threshold = 0.1;

            s = sprintf('%.2g', maxAmplitude);
            set(max_amplitude_edit,'string',s);
            s = sprintf('%.2g', minAmplitude);
            set(min_amplitude_edit,'string',s);
            s = sprintf('%.2g', threshold);
            set(threshold_edit,'string',s);
                       
            % make sure data is plotted showing normalized range and
            % disable scaling
            set(scaleUpArrow,'enable','off');
            set(scaleDownArrow,'enable','off');
            set(autoScaleButton,'enable','off');
            set(min_scale,'enable','off');
            set(max_scale,'enable','off');
            
        else
            setMarkingCtls('off');
            % enable scaling
            set(scaleUpArrow,'enable','on');
            set(scaleDownArrow,'enable','on');
            set(autoScaleButton,'enable','on');
            set(min_scale,'enable','on');
            set(max_scale,'enable','on');
        end
        
        drawTrial;
        autoScale_callback;
        updateCursors;
             
    end



% +++++++++ amplitude and time controls +++++++++

cursor_text = uicontrol('style','text','fontsize',12,'units','normalized','position',[0.82 0.31 0.1 0.02],...
     'string','Cursor =','BackgroundColor','white','foregroundColor',[0.7,0.41,0.1]);

latency_slider = uicontrol('style','slider','units', 'normalized',...
    'position',[0.05 0.25 0.9 0.02],'min',0,'max',1,'Value',0,...
    'sliderStep', [0 1],'BackGroundColor',[0.8 0.8 0.8],'ForeGroundColor',...
    'white'); 

% this callback is called everytime slider value is changed. 
% replaces slider callback function

addlistener(latency_slider,'Value','PostSet',@slider_moved_callback);

    function slider_moved_callback(~,~)   
        % slider from 0 to 1 but data may be negative
       val = get(latency_slider,'Value');
       epochStart = header.epochMinTime + (val * header.trialDuration);
       
       if (epochStart + epochTime > header.epochMaxTime)
            epochStart = header.epochMaxTime - epochTime;
       end
       
       if epochStart < header.epochMinTime
           epochStart = header.epochMinTime;
       end
       drawTrial;
    end


scaleMenuItems = {'<HTML><FONT COLOR="blue">MEG</HTML>'; ...
        '<HTML><FONT COLOR="rgb(26,179,26)">ADC</HTML>';...
        '<HTML><FONT COLOR="rgb(204,128,26)">Trigger</HTML>';...
        '<HTML><FONT COLOR="black">Digital</HTML>';...
         '<HTML><FONT COLOR="cyan">EEG/EMG</HTML>';...       
        '<HTML><FONT COLOR="magenta">Other</HTML>'};

% ++++++++++++ scale menu and controls ...


setGoodButton = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.05 0.97 0.05 0.02],...
    'Foregroundcolor','black','string','Set Good','backgroundcolor','white','callback',@setGood);

    function setGood(~,~)
        idx = find(selectedMask == 1);
        
        if idx < 1
            return;
        end
        chanIdx = selectedChannelList(idx);
        badChannelMask(chanIdx) = 0;
        % deselect so we can see channel status change.
        selectedMask(:) = 0;
        drawTrial;
    end

setBadButton = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.11 0.97 0.05 0.02],...
    'Foregroundcolor','black','string','Set Bad','backgroundcolor','white','callback',@setBad);
    function setBad(~,~)
        idx = find(selectedMask == 1);
        if idx < 1
            return;
        end
        chanIdx = selectedChannelList(idx);
        badChannelMask(chanIdx) = 1;
        % deselect so we can see channel status change.
        selectedMask(:) = 0;
        
        drawTrial;
    end


setBadTrialButton = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.18 0.97 0.1 0.02],...
    'Foregroundcolor','black','string','Set Trial Good/Bad','backgroundcolor','white','callback',@setTrialBad);
    function setTrialBad(~,~)

        if header.numTrials < 2
            return;
        end        
        badTrialMask(trialNo) = ~badTrialMask(trialNo);
        
        drawTrial;
    end


overlayPlotsCheck = uicontrol('style','checkbox','units','normalized','position',[0.31 0.97 0.1 0.02],...
    'string','Overlay Plots','backgroundcolor','white','value',overlayPlots,'FontSize',11,'callback',@overlay_plots_callback);

    function overlay_plots_callback(src,~)
        overlayPlots=get(src,'value');
        drawTrial;
    end

uicontrol('style','checkbox','units','normalized','position',[0.39 0.97 0.1 0.02],...
    'string','Show Amplitude','backgroundcolor','white','value',showAmplitudes,'FontSize',11,'callback',@show_amplitudes_callback);

    function show_amplitudes_callback(src,~)
        showAmplitudes=get(src,'value');
        drawTrial;
    end


numColumnsMenu = uicontrol('style','popupmenu','units','normalized','fontsize',11,'position',[0.47 0.95 0.08 0.04],...
  'Foregroundcolor','black','string',{'1 Column';'2 Columns';'3 Columns'; '4 Columns'},'value',...
            numColumns,'backgroundcolor','white','callback',@column_number_callback);
             
    function column_number_callback(src,~)
        numColumns = get(src,'value');
        drawTrial;      
    end

uicontrol('style','popupmenu','units','normalized','fontsize',11,'position',[0.04 0.21 0.06 0.03],...
  'Foregroundcolor','black','string',scaleMenuItems,'value',...
            currentScaleMenuIndex,'backgroundcolor','white','callback',@scaleMenu_callback);
             
    function scaleMenu_callback(src,~)
        
        % update fields
        currentScaleMenuIndex = get(src,'value');

        maxScale = maxRange(currentScaleMenuIndex);
        minScale = minRange(currentScaleMenuIndex);
        s = sprintf('%.3g', maxScale);
        set(max_scale,'string',s);
        s = sprintf('%.3g', minScale);
        set(min_scale,'string',s);
        
    end

uicontrol('style','text','units','normalized','position',[0.1 0.22 0.1 0.02],...
    'string','Max:','fontsize',11,'backgroundcolor','white','horizontalalignment','left');
max_scale=uicontrol('style','edit','units','normalized','position',[0.12 0.22 0.04 0.025],...
    'FontSize', 11, 'BackGroundColor','white','string',maxScale,...
    'callback',@max_scale_callback);

    function max_scale_callback(src,~)
        val =str2double(get(src,'string'));

        maxScale = maxRange(currentScaleMenuIndex);
        minScale = minRange(currentScaleMenuIndex);
        if val < minScale
            s = sprintf('%.3g', maxScale);
            set(max_scale,'string',s);
            return;
        end

        maxRange(currentScaleMenuIndex) = val;
        drawTrial;
    end

uicontrol('style','text','units','normalized','position',[0.17 0.22 0.08 0.02],...
    'string','Min:','fontsize',11,'backgroundcolor','white','horizontalalignment','left');
min_scale=uicontrol('style','edit','units','normalized','position',[0.19 0.22 0.04 0.025],...
    'FontSize', 11, 'BackGroundColor','white','string',minScale,...
    'callback',@min_scale_callback);

    function min_scale_callback(src,~)
        val = str2double(get(src,'string'));
        maxScale = maxRange(currentScaleMenuIndex);
        minScale = minRange(currentScaleMenuIndex);
        if val > maxScale
            s = sprintf('%.3g', minScale);
            set(min_scale,'string',s);
            return;
        end
        minRange(currentScaleMenuIndex) = val;
        drawTrial;
    end
     
scaleUpArrow = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.24 0.22 0.025 0.025],...
    'CData',uparrow_im,'Foregroundcolor','black','backgroundcolor','white','callback',@scaleUp_callback);

scaleDownArrow = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.27 0.22 0.025 0.025],...
    'CData',downarrow_im,'Foregroundcolor','black','backgroundcolor','white','callback',@scaleDown_callback);

autoScaleButton = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.31 0.22 0.06 0.025],...
    'Foregroundcolor','black','string','Autoscale','backgroundcolor','white','callback',@autoScale_callback);
    
    function scaleUp_callback(~,~)
        if isempty(header)
            return;
        end
        maxScale = maxRange(currentScaleMenuIndex);
        minScale = minRange(currentScaleMenuIndex);

        inc = maxScale * 0.1;
        maxScale = maxScale  - inc;
        minScale = minScale + inc;
        s = sprintf('%.3g', maxScale);
        set(max_scale,'String', s);
        s = sprintf('%.3g', minScale);
        set(min_scale,'String', s);
        
        maxRange(currentScaleMenuIndex) = maxScale;
        minRange(currentScaleMenuIndex) = minScale;
        
        drawTrial;
    end

    function scaleDown_callback(~,~)
        if isempty(header)
            return;
        end
        maxScale = maxRange(currentScaleMenuIndex);
        minScale = minRange(currentScaleMenuIndex);

        inc = maxScale * 0.1;
        maxScale = maxScale  + inc;
        minScale = minScale - inc;
        s = sprintf('%.3g', maxScale);
        set(max_scale,'String', s);
        s = sprintf('%.3g', minScale);
        set(min_scale,'String', s);
               
        maxRange(currentScaleMenuIndex) = maxScale;
        minRange(currentScaleMenuIndex) = minScale;

        drawTrial;
    end


    function autoScale_callback(~,~)
        if isempty(header)
            return;
        end
        maxRange(currentScaleMenuIndex) = NaN;
        minRange(currentScaleMenuIndex) = NaN;
        drawTrial;
    end

trialNumTxt = uicontrol('style','text','fontsize',12,'units','normalized','horizontalalignment','left','position',...
     [0.38 0.21 0.08 0.03],'string','Trial: 1 of 1','BackgroundColor','white','foregroundcolor','black');

trialStartButton = uicontrol('style','pushbutton','units','normalized','fontsize',14,'fontweight','bold','position',[0.43 0.22 0.02 0.025],...
    'enable','off','String','<<','Foregroundcolor','black','backgroundcolor','white','callback',@trial_start_callback);
trialDecButton = uicontrol('style','pushbutton','units','normalized','fontsize',14,'fontweight','bold','position',[0.455 0.22 0.02 0.025],...
    'enable','off','String','<','Foregroundcolor','black','backgroundcolor','white','callback',@trial_dec_callback);
trialIncButton = uicontrol('style','pushbutton','units','normalized','fontsize',14,'fontweight','bold','position',[0.48 0.22 0.02 0.025],...
    'enable','off', 'String', '>','Foregroundcolor','black','backgroundcolor','white','callback',@trial_inc_callback);

trialEndButton = uicontrol('style','pushbutton','units','normalized','fontsize',14,'fontweight','bold','position',[0.505 0.22 0.02 0.025],...
    'enable','off', 'String', '>>','Foregroundcolor','black','backgroundcolor','white','callback',@trial_end_callback);


    function trial_inc_callback(~,~)  
        if isempty(header)
            return;
        end
        trialNo = trialNo + 1;
        if trialNo > header.numTrials
            trialNo = header.numTrials;
        end
        if trialNo < 1 
            trialNo = 1;
        end

        loadData;
        drawTrial;

    end

    function trial_dec_callback(~,~)  
        if isempty(header)
            return;
        end
        
        trialNo = trialNo - 1;
        if trialNo > header.numTrials
            trialNo = header.numTrials;
        end
        if trialNo < 1 
            trialNo = 1;
        end
        loadData;
        drawTrial;
    end

    function trial_start_callback(~,~)  
        if isempty(header)
            return;
        end
        trialNo = 1;
        loadData;
        drawTrial;
    end

    function trial_end_callback(~,~)  
        if isempty(header)
            return;
        end
        trialNo = header.numTrials;
        loadData;
        drawTrial;
    end


function edit_markers_callback(~,~)    
    bw_editCTFMarkers(dsName);
    loadMarkerFile(markerFileName);
end

markerNames = {'none'};
uicontrol('style','text','fontsize',12,'units','normalized','horizontalalignment','left','position',...
     [0.535 0.21 0.08 0.03],'string','Markers:','BackgroundColor','white','foregroundcolor','black');
marker_Popup =uicontrol('style','popup','units','normalized','fontsize',12,'position',[0.57 0.19 0.08 0.05],...
    'string',markerNames,'value',currentMarkerIndex,'Foregroundcolor','black','backgroundcolor','white','callback',@marker_popup_callback);
% 
% markerDecButton = uicontrol('style','pushbutton','units','normalized','fontsize',14,'fontweight','bold','position',[0.42 0.22 0.02 0.025],...
%     'enable','off','String','<<','Foregroundcolor','black','backgroundcolor','white','callback',@trial_start_callback);

markerStartButton = uicontrol('style','pushbutton','units','normalized','fontsize',14,'fontweight','bold','position',[0.65 0.22 0.02 0.025],...
    'enable','off','String','<<','Foregroundcolor','black','backgroundcolor','white','callback',@marker_start_callback);
markerDecButton = uicontrol('style','pushbutton','units','normalized','fontsize',14,'fontweight','bold','position',[0.675 0.22 0.02 0.025],...
    'enable','off','String','<','Foregroundcolor','black','backgroundcolor','white','callback',@marker_dec_callback);
markerIncButton = uicontrol('style','pushbutton','units','normalized','fontsize',14,'fontweight','bold','position',[0.7 0.22 0.02 0.025],...
    'enable','off', 'String', '>','Foregroundcolor','black','backgroundcolor','white','callback',@marker_inc_callback);
markerEndButton = uicontrol('style','pushbutton','units','normalized','fontsize',14,'fontweight','bold','position',[0.725 0.22 0.02 0.025],...
    'enable','off', 'String', '>>','Foregroundcolor','black','backgroundcolor','white','callback',@marker_end_callback);

% trialEndButton = uicontrol('style','pushbutton','units','normalized','fontsize',14,'fontweight','bold','position',[0.495 0.22 0.02 0.025],...
%     'enable','off', 'String', '>>','Foregroundcolor','black','backgroundcolor','white','callback',@trial_end_callback);
% 
    function marker_popup_callback(src,~)    
        currentMarkerIndex = get(src,'value');
        if currentMarkerIndex == 1
            set(markerIncButton,'enable','off')
            set(markerDecButton,'enable','off')
            set(markerStartButton,'enable','off')
            set(markerEndButton,'enable','off')
        else
            set(markerIncButton,'enable','on')
            set(markerDecButton,'enable','on')
            set(markerStartButton,'enable','on')
            set(markerEndButton,'enable','on')
        end
        drawTrial;
    end

    function marker_inc_callback(~,~)  
        if currentMarkerIndex == 1
            return;
        end
        markerTimes = markerLatencies{currentMarkerIndex-1};
        if currentMarker < numel(markerTimes)
            currentMarker = currentMarker + 1;
            epochStart = markerTimes(currentMarker) - (epochTime / 2);
            cursorLatency = markerTimes(currentMarker);

            drawTrial;       
            updateCursors;
            % adjust slider position
            val = (epochStart + header.epochMinTime) / header.trialDuration;
            if val < 0, val = 0.0; end
            if val > 1.0, val = 1.0; end
            set(latency_slider, 'value', val);
        end
    end

    function marker_dec_callback(~,~)  
        if currentMarkerIndex == 1
            return;
        end
        markerTimes = markerLatencies{currentMarkerIndex-1};
        if currentMarker > 1
            currentMarker = currentMarker - 1;
            epochStart = markerTimes(currentMarker) - (epochTime / 2);
            cursorLatency = markerTimes(currentMarker);
            drawTrial;    
            updateCursors;
            % adjust slider position
            val = (epochStart + header.epochMinTime) / header.trialDuration;
            if val < 0, val = 0.0; end
            if val > 1.0, val = 1.0; end
            set(latency_slider, 'value', val);
        end
    end

    function marker_start_callback(~,~)  
        if currentMarkerIndex == 1
            return;
        end
        markerTimes = markerLatencies{currentMarkerIndex-1};
        currentMarker = 1;
        epochStart = markerTimes(currentMarker) - (epochTime / 2);
        cursorLatency = markerTimes(currentMarker);

        drawTrial;       
        updateCursors;
        % adjust slider position
        val = (epochStart + header.epochMinTime) / header.trialDuration;
        if val < 0, val = 0.0; end
        if val > 1.0, val = 1.0; end
        set(latency_slider, 'value', val);
    end

    function marker_end_callback(~,~)  
        if currentMarkerIndex == 1
            return;
        end
        markerTimes = markerLatencies{currentMarkerIndex-1};
        currentMarker = numel(markerTimes);
        epochStart = markerTimes(currentMarker) - (epochTime / 2);
        cursorLatency = markerTimes(currentMarker);

        drawTrial;       
        updateCursors;
        % adjust slider position
        val = (epochStart + header.epochMinTime) / header.trialDuration;
        if val < 0, val = 0.0; end
        if val > 1.0, val = 1.0; end
        set(latency_slider, 'value', val);
    end

% window duration
uicontrol('style','text','units','normalized','position',[0.77 0.22 0.1 0.02],...
    'string','Window Duration (s):','fontsize',11,'backgroundcolor','white','horizontalalignment','left');
epochDurationEdit = uicontrol('style','edit','units','normalized','position',[0.85 0.225 0.05 0.02],...
    'FontSize', 11, 'BackGroundColor','white','string',epochTime,...
    'callback',@epoch_duration_callback);

    function epoch_duration_callback(src,~)
        t =str2double(get(src,'string'));

        if isnan(t) || t >= header.trialDuration
            t = header.trialDuration;
            set(src,'string',t);
            epochStart = header.epochMinTime;
        end    
        
        % minimum = two samples
        if t < (1.0 / header.sampleRate) * 2.0
            t = (1.0 / header.sampleRate) * 2.0;
            set(src,'string',t);
        end

        epochTime = t;
        drawTrial;
        updateSlider;

    end

uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.91 0.225 0.05 0.02],...
    'String','Whole trial','Foregroundcolor','black','backgroundcolor','white','callback',@whole_trial_callback);
    function whole_trial_callback(~,~)
        if isempty(header)
            return;
        end
        epochStart = header.epochMinTime;
        epochTime = header.trialDuration;
        drawTrial;    
        updateSlider;
        set(epochDurationEdit,'string',header.trialDuration);
    end

    function updateSlider
        epochSamples = round(epochTime * header.sampleRate);
        dataRange = epochTime / header.trialDuration;        
        sliderScale =[dataRange/10 dataRange];
        set(latency_slider,'sliderStep',sliderScale);
        drawTrial;
    end

% +++++++++++++++++++++++++

% ++++++++++ plot settings 

annotation('rectangle',[0.05 0.02 0.39 0.18],'EdgeColor','blue');
uicontrol('style','text','fontsize',11,'units','normalized','position',...
     [0.08 0.185 0.1 0.025],'string','Plot Settings','BackgroundColor','white','foregroundcolor','blue','fontweight','b');

sampleRateTxt = uicontrol('style','text','units','normalized','position',[0.06 0.16 0.08 0.02],...
    'string','Sample Rate:','backgroundcolor','white','FontSize',11, 'HorizontalAlignment','Left');

totalSamplesTxt = uicontrol('style','text','units','normalized','position',[0.16 0.16 0.08 0.02],...
    'string','Total Samples:','backgroundcolor','white','FontSize',11, 'HorizontalAlignment','Left');

totalDurationTxt = uicontrol('style','text','units','normalized','position',[0.25 0.16 0.09 0.02],...
    'string','Trial Duration:','backgroundcolor','white','FontSize',11, 'HorizontalAlignment','Left');

numTrialsTxt = uicontrol('style','text','units','normalized','position',[0.35 0.16 0.06 0.02],...
    'string','# of Trials:','backgroundcolor','white','FontSize',11, 'HorizontalAlignment','Left');


numChannelsTxt = uicontrol('style','text','units','normalized','position',[0.06 0.14 0.08 0.02],...
    'string','Total Channels:','backgroundcolor','white','FontSize',11, 'HorizontalAlignment','Left');

numSensorsTxt = uicontrol('style','text','units','normalized','position',[0.16 0.14 0.08 0.02],...
    'string','MEG Sensors:','backgroundcolor','white','FontSize',11, 'HorizontalAlignment','Left');

numReferencesTxt = uicontrol('style','text','units','normalized','position',[0.25 0.14 0.08 0.02],...
    'string','MEG References:','backgroundcolor','white','FontSize',11, 'HorizontalAlignment','Left');

numAnalogTxt = uicontrol('style','text','units','normalized','position',[0.35 0.14 0.06 0.02],...
    'string','ADC Channels:','backgroundcolor','white','FontSize',11, 'HorizontalAlignment','Left');



filterCheckBox = uicontrol('style','checkbox','units','normalized','position',[0.08 0.1 0.04 0.02],...
    'string','Filter','backgroundcolor','white','value',~filterOff,'FontSize',11,'callback',@filter_check_callback);

% replace with 50 / 60 Hz checks - check code 
uicontrol('style','checkbox','units','normalized','position',[0.08 0.07 0.08 0.02],...
    'string','60 Hz Notch','backgroundcolor','white','value',notchFilter,'FontSize',11,'callback',@notch_check_callback);
uicontrol('style','checkbox','units','normalized','position',[0.17 0.07 0.08 0.02],...
    'string','50 Hz Notch','backgroundcolor','white','value',notchFilter2,'FontSize',11,'callback',@notch2_check_callback);

uicontrol('style','checkbox','units','normalized','position',[0.26 0.07 0.1 0.02],...
    'string','Remove Offset','backgroundcolor','white','value',removeOffset,'FontSize',11,'callback',@remove_offset_callback);

uicontrol('style','checkbox','units','normalized','position',[0.08 0.04 0.08 0.02],...
    'string','Invert','backgroundcolor','white','value',invertData,'FontSize',11,'callback',@invert_check_callback);

uicontrol('style','checkbox','units','normalized','position',[0.17 0.04 0.08 0.02],...
    'string','Rectify','backgroundcolor','white','value',rectify,'FontSize',11,'callback',@rectify_check_callback);

uicontrol('style','checkbox','units','normalized','position',[0.26 0.04 0.1 0.02],...
    'string','Differentiate','backgroundcolor','white','value',differentiate,'FontSize',11,'callback',@firstDiff_check_callback);

% uicontrol('style','checkbox','units','normalized','position',[0.34 0.04 0.05 0.02],...
%     'string','Envelope','backgroundcolor','white','value',envelope,'FontSize',11,'callback',@envelope_check_callback);

uicontrol('style','text','units','normalized','position',[0.15 0.1 0.08 0.02],...
    'string','High Pass (Hz):','fontsize',11,'backgroundcolor','white','horizontalalignment','left');
filt_hi_pass=uicontrol('style','edit','units','normalized','position',[0.22 0.1 0.04 0.02],...
    'FontSize', 11, 'BackGroundColor','white','string',bandPass(1),...
    'callback',@filter_hipass_callback);

    function filter_hipass_callback(src,~)
        bandPass(1)=str2double(get(src,'string'));
        processData;
        drawTrial;
    end

uicontrol('style','text','units','normalized','position',[0.28 0.1 0.08 0.02],...
    'string','Low Pass (Hz):','fontsize',11,'backgroundcolor','white','horizontalalignment','left');
filt_low_pass=uicontrol('style','edit','units','normalized','position',[0.35 0.1 0.04 0.02],...
    'FontSize', 11, 'BackGroundColor','white','string',bandPass(2),...
    'callback',@filter_lowpass_callback);

    function filter_lowpass_callback(src,~)
        bandPass(2)=str2double(get(src,'string'));
        processData;
        drawTrial;
    end

if filterOff
    set(filt_hi_pass, 'enable','off');
    set(filt_low_pass, 'enable','off');
else
    set(filt_hi_pass, 'enable','on');
    set(filt_low_pass, 'enable','on');
end


showMapCheck = uicontrol('style','checkbox','units','normalized','position',[0.46 0.185 0.08 0.02],...
    'string','show Topoplot','backgroundcolor','white','value',showMap,'FontSize',11,'callback',@show_map_callback);

    function show_map_callback(src,~)
        showMap = get(src,'value');
        if ~showMap
            currentAx = gca;
            subplot('Position',mapbox);
            cla;
            axes(currentAx);
%             subplot('Position',plotbox);            
        else
            updateMap;
        end
    end


function rectify_check_callback(src,~)
    rectify=get(src,'value');
    processData;
    drawTrial;
    updateMap;
end

function envelope_check_callback(src,~)
    envelope=get(src,'value');
    processData;
    drawTrial;
    updateMap;
end

function invert_check_callback(src,~)
    invertData=get(src,'value');
    processData;
    drawTrial;
    updateMap;

end

function firstDiff_check_callback(src,~)
    differentiate=get(src,'value');
    processData;
    drawTrial;
    updateMap;
end

function notch_check_callback(src,~)
    notchFilter=get(src,'value');
    processData;
    drawTrial;
    updateMap;
end

function notch2_check_callback(src,~)
    notchFilter2=get(src,'value');
    processData;
    drawTrial;
    updateMap;
end

function remove_offset_callback(src,~)
    removeOffset=get(src,'value');
    processData;
    drawTrial;
    updateMap;
end


function filter_check_callback(src,~)
    filterOff=~get(src,'value');
    if filterOff
        set(filt_hi_pass, 'enable','off');
        set(filt_low_pass, 'enable','off');
    else
        bandPass(1)=str2double(get(filt_hi_pass,'string'));
        bandPass(2)=str2double(get(filt_low_pass,'string'));
        set(filt_hi_pass, 'enable','on');
        set(filt_low_pass, 'enable','on');
    end
    
    processData;
    drawTrial;
    updateMap;
end

function find_events_callback(~,~)
    markData;
    drawTrial;
end

function reverse_check_callback(~,~)
    reverseScan = 1;
    set(forward_scan_radio,'value',0)

end

function forward_check_callback(~,~)
    reverseScan = 0;
    set(reverse_scan_radio,'value',0)
end

function loadData
    % load (all) data with current filter settings etc and adjust scale

    readAllChannels = 1;
            
    s = sprintf('Trial: %d of %d', trialNo, header.numTrials);
    set(trialNumTxt,'string',s);

    % read each trial of raw data once...      
    % ** for mex function sample and trial indices begin at zero ! 
    all_data = bw_getCTFData(dsName, 0, header.numSamples, trialNo-1, readAllChannels)';  
    
    timeVec = header.epochMinTime: 1/ header.sampleRate: header.epochMaxTime;
    timeVec = timeVec(1:header.numSamples);  % should be same ...

    dataarray = all_data; 

end

function processData
    
    % reload all data
    if isempty(dataarray)
        return;
    end
   
    % process all displayed channels;
    channelsToProcess = selectedChannelList;
    if size(channelsToProcess,1) > 1
        channelsToProcess = channelsToProcess';
    end
    if size(meg_idx,1) > 1
        meg_idx = meg_idx';
    end
    % if showing map need to process all MEG sensors
    if showMap
        idx = [channelsToProcess meg_idx];
        channelsToProcess = unique(idx);
    end
    
    % only filter analog channels
    aTypes = [MEG_CHANNELS ADC_CHANNELS EEG_CHANNELS];
    chanTypes = channelTypes(channelsToProcess);    % returns list by channel type
    includeMask = ismember(chanTypes,aTypes);       % binary mask of analog channels
    analogChannelsToProcess = channelsToProcess(includeMask);
    
    % start with raw data and apply processing in this order once
    % - will also undo processsing
    dataarray = all_data;
        
    % filter data
    if ~filterOff && ~isempty(analogChannelsToProcess)       
        
        data = dataarray(analogChannelsToProcess,1:header.numSamples);
        nchan = size(data,1);
        
        % input and output of bw_filter has to be nsamples x nchannels
        % need to transpose input and output
        fdata = bw_filter(data', header.sampleRate, bandPass)';
        
        dataarray(analogChannelsToProcess,1:header.numSamples) = fdata(1:nchan,1:header.numSamples);    
            
    end
    
    if notchFilter || notchFilter2 && ~isempty(analogChannelsToProcess)       
            
        lineFilterWidth = 3;        
        if notchFilter 
            lineFilterFreq = 60.0;            
        else
             lineFilterFreq = 50.0;         
        end
        
        data = dataarray(analogChannelsToProcess,1:header.numSamples);
        nchan = size(data,1);

        f1 = -lineFilterWidth;
        f2 = lineFilterWidth;
        fdata = [];
        for j=1:4   % remove fundamental and 1st 3 harmonics
            f1 = f1 + lineFilterFreq;
            f2 = f2 + lineFilterFreq;
            if f2 < header.sampleRate / 2.0               
                fdata = bw_filter(data', header.sampleRate, [f1 f2],4,1,1)';    
                data = fdata;
            end
        end
        
        dataarray(analogChannelsToProcess,1:header.numSamples) = fdata(1:nchan,1:header.numSamples);                       

    end
    
    if differentiate
        data = diff(dataarray(channelsToProcess,1:header.numSamples),1,2);
        dataarray(channelsToProcess,1:header.numSamples) = [data zeros(size(data,1),1)];  % diff drops one sample
    end
    
    if removeOffset
        offset = mean(dataarray(channelsToProcess,1:header.numSamples),2);
        dataarray(channelsToProcess,1:header.numSamples) = dataarray(channelsToProcess,1:header.numSamples) - offset;
    end

    if rectify
        dataarray(channelsToProcess,1:header.numSamples) = abs(dataarray(channelsToProcess,1:header.numSamples));
    end

    if invertData
        dataarray(channelsToProcess,1:header.numSamples) = dataarray(channelsToProcess,1:header.numSamples) * -1.0;
    end
            
    if envelope && ~isempty(analogChannelsToProcess)
        dataarray(analogChannelsToProcess,1:header.numSamples) = abs(hilbert(dataarray(channelsToProcess,1:header.numSamples)));
    end
    

end

    function drawTrial
       
        if isempty(header)
            return;
        end
        %  ** here have to loop over channels and set colours and scales separately ***       
        numChannelsToDisplay = numel(selectedChannelList);
         
        
        % initialize scales for different channels
        for k=1:numChannelsToDisplay
            chanIdx = selectedChannelList(k);
            chanType = channelTypes(chanIdx);
                     
            % autoscale this channel type?
                           
            % get current max for this channel type

            for j=1:numel(maxRange)
                includeChannel = 0;
                switch j
                    case 1
                        includeChannel = ismember(chanType,MEG_CHANNELS);
                    case 2
                        includeChannel = ismember(chanType,ADC_CHANNELS);
                    case 3 
                        includeChannel =  ismember(chanType,TRIGGER_CHANNELS);                     
                    case 4
                        includeChannel =  ismember(chanType,DIGITAL_CHANNELS);                                    
                    case 5
                        includeChannel =  ismember(chanType,EEG_CHANNELS);
                    case 6
                        includeChannel =  ismember(chanType,OTHER_CHANNELS);
                end
                if includeChannel                       
                     mx = maxRange(j);
                     if isnan(mx) || mx == 0.0
                        
                        [timebase, fd] = getTrial(chanIdx, epochStart);                      
                        maxRange(j) = max(abs(fd));
                        minRange(j) = -maxRange(j);        
                     end
                        
                end
            end
            
            
        end

        % normalize plot data for plotting channels together.
        % normalize to 0.5 of full plot range for readability 

        if numChannelsToDisplay < numColumns
            numColumns = 1;
            set(numColumnsMenu,'value',1);
        end
        
        numChannelsPerColumn = floor(numChannelsToDisplay / numColumns);
        if mod(numChannelsToDisplay,numColumns)         % increment by one if not even number
            numChannelsPerColumn = numChannelsPerColumn + 1;
        end
           
        plotMax = 2.0;              
        plotMin = -2.0;

        if ~overlayPlots
            plotMax = plotMax * numChannelsPerColumn;
            plotMin = plotMin * numChannelsPerColumn;
        end
        % add an offset to plotted data to stack plots
        singleRange = ((plotMax-plotMin) / numChannelsPerColumn);
        
        % plot channels
        amplitudeLabels = [];               
        
        switch numColumns
            case 1
                tbox = plotbox;
                cinc = 0;
            case 2
                tbox = [0.05 0.31 0.41 0.65];
                cinc = tbox(3) + 0.07;
            case 3
                tbox = [0.05 0.31 0.25 0.65];
                cinc = tbox(3) + 0.07;            
            case 4
                tbox = [0.05 0.31 0.17 0.65];
                cinc = tbox(3) + 0.07;
        end
       
        plotCount = 0;
        cursorHandles = [];
        
        
        % +++ plot data ++++
        for col=1:numColumns
                           
            start = plotMin + (singleRange / 2.0); 
            if col > 1 
                tbox(1) = tbox(1) + cinc;
            end
            subplot('Position',tbox);                       
            hold off;
                            
            for k=1:numChannelsPerColumn
                              
                plotCount = plotCount + 1;
                
                % if uneven number 
                if plotCount > numChannelsToDisplay 
                    break;                   
                end
                
                chanIdx = selectedChannelList(plotCount);
                chanType = channelTypes(chanIdx);

                [timebase, fd] = getTrial(chanIdx, epochStart); 

                switch chanType
                    case num2cell(MEG_CHANNELS)
                        plotColour = 'blue';
                        amplitudeUnits(plotCount)  = {'Tesla'};
                        maxAmp = maxRange(1);
                    case num2cell(ADC_CHANNELS)
                        plotColour = darkGreen;
                        amplitudeUnits(plotCount) = {'Volts'};
                        maxAmp = maxRange(2);
                    case num2cell(TRIGGER_CHANNELS) 
                        plotColour = orange;
                        amplitudeUnits(plotCount) = {' '};
                        maxAmp = maxRange(3);
                    case num2cell(DIGITAL_CHANNELS)
                        plotColour = 'black';
                        amplitudeUnits(plotCount) = {'Bits'};
                        maxAmp = maxRange(4);
                    case num2cell(EEG_CHANNELS)
                        plotColour = 'cyan';
                        amplitudeUnits(plotCount) = {'Volts'};
                        maxAmp = maxRange(5);
                    otherwise
                        plotColour = 'magenta';
                        amplitudeUnits(plotCount) = {' '};
                        maxAmp = maxRange(6);
                end

                if badChannelMask(chanIdx) == 1 || badTrialMask(trialNo) == 1
                    plotColour = [0.8 0.8 0.8];
                end

                if selectedMask(plotCount) == 1
                    plotColour = 'red';
                end

                pd = fd ./ maxAmp;
                channelName = char( channelNames(chanIdx) );

                % add offset...       
                if ~overlayPlots
                    offset =  start + ((numChannelsPerColumn-k) * singleRange);
                else
                    offset = 0.0;
                end

                pd = pd + offset;
                plot(timebase,pd,'Color',plotColour);

                % draw baselines if stacked or once for overlays

                if ~overlayPlots || k == 1
                    v = [offset offset];
                    h = xlim;            
                    line(h,v, 'color', 'black');
                end
                
                % plot channel Name and amplitude
                if  ~overlayPlots
                    if numColumns == 1
                        coffset = epochTime * 0.03;
                    elseif numColumns == 2
                        coffset = epochTime * 0.07;
                    elseif numColumns == 3
                        coffset = epochTime * 0.11;
                    else
                        coffset = epochTime * 0.13;                        
                    end
                    x = epochStart - coffset;
                    y = offset;
                    if ~enableMarking
                        s = sprintf('%s', channelName);
                        text(x,y,s,'color',plotColour,'interpreter','none','fontsize',8,'UserData',plotCount,'ButtonDownFcn',@clickedOnLabel);
                    end

                    if showAmplitudes
                        sample = round( (cursorLatency - header.epochMinTime - epochStart) * header.sampleRate) + 1;
                        if sample > length(fd), sample = length(fd); end
                        if sample < 1, sample = 1; end
                        amplitude = fd(sample);
                        s = sprintf('%.2g %s', amplitude, char(amplitudeUnits(k)));
                        x = epochStart + epochTime + (epochTime * 0.007);
                        amplitudeLabels(plotCount) = text(x,y,s,'color',plotColour,'interpreter','none','fontsize',8);
                    end

                end
                hold on;                                   
            end            

           % annotate plots
           
           ylim([plotMin plotMax])
           xlim([timebase(1) timebase(end)])
           
           idx = find(selectedMask == 1);
           if isempty(idx) 
               set(setGoodButton,'enable','off');
               set(setBadButton,'enable','off');
           else
               set(setGoodButton,'enable','on');
               set(setBadButton,'enable','on');
           end

            ylim([plotMin plotMax]);

            if enableMarking
                set(gca,'ytick',[-1.0 -0.5 0.5 1.0])
                ylabel('Threshold (normalized)')
            else
                set(gca,'ytick',[])
            end

            % plot xscale and other markers
            xlim([timebase(1) timebase(end)])
            xlabel('Time (sec)', 'fontsize', 12);

            % check if events exist in this window and draw...   
                       
            % draw events
            
            nsamples = length(timebase);
            if enableMarking

                th = ones(nsamples,1) * threshold';

                plot(timebase, th', 'r:', 'lineWidth',1.5);

                th = ones(nsamples,1) * minAmplitude';
                plot(timebase, th', 'g:', 'lineWidth',1.5);

                th = ones(nsamples,1) * maxAmplitude';
                plot(timebase, th', 'c:', 'lineWidth',1.5);
            end
            
            if ~isempty(eventList)
                events = find(eventList > epochStart & eventList < (epochStart+epochTime));

                if ~isempty(events)
                    for i=1:length(events)
                        thisEvent = events(i);
                        t = eventList(thisEvent);
                        h = [t,t];
                        v = ylim;
                        line(h,v, 'color', orange,'linewidth',1);
                        if thisEvent == currentEvent                         
                            latency = eventList(thisEvent);
                            x = t + epochTime * 0.005;
                            y =  v(1) + (v(2) - v(1))*0.05;
                            s = sprintf('Event # %d  (latency = %.4f s)', currentEvent, latency);
                            text(x,y,s,'color',orange);
                        end                  
                    end
                end
            end
           
           % check if markers exist in this window and draw...   
            if currentMarkerIndex > 1 && currentMarkerIndex < numMarkers+2      % marker menu count = "none" + numMarkers
                markerTrialNo = markerTrials{currentMarkerIndex-1};

                markerTimes = markerLatencies{currentMarkerIndex-1};
                markerName = char( markerNames{currentMarkerIndex} );

                % get valid markers, must be in time window and trial
                markers = find(markerTimes > epochStart & markerTimes < (epochStart+epochTime) );
                if ~isempty(markers)              
                    for k=1:length(markers)                       
                        % draw this marker latnecy if this is its trial number
                        if markerTrialNo(markers(k)) ~= trialNo
                            continue;
                        end

                        t = markerTimes(markers(k));
                        h = [t,t];
                        v = ylim;

                        line(h,v, 'color', 'blue');
                        x = t + epochTime * 0.001;
                        y =  v(2) - (v(2) - v(1))*0.05;
                        s = sprintf('%s', markerName);
                        text(x,y,s,'color','blue','interpreter','none');
                    end
                end            
            elseif currentMarkerIndex == numMarkers+2       % draw all markers
                for j=1:numMarkers
                    markerTimes = markerLatencies{j};
                    markerName = char( markerNames{j+1} );
                    markers = find(markerTimes > epochStart & markerTimes < (epochStart+epochTime));

                    if ~isempty(markers)
                        for k=1:length(markers)
                            t = markerTimes(markers(k));
                            h = [t,t];
                            v = ylim;

                            line(h,v, 'color', 'blue');
                            x = t + epochTime * 0.001;
                            y =  v(2) - (v(2) - v(1))*0.05;
                            s = sprintf('%s', markerName);
                            text(x,y,s,'color','blue','interpreter','none');
                        end
                    end                                    
                end
            end

            s = sprintf('%s',char( channelSets(channelMenuIndex)));
            if enableMarking
                tt = legend(channelName,'threshold', 'min. amplitude', 'max. amplitude');
                set(tt,'interpreter','none','Autoupdate','off');
            end

            plotAxes(col) = gca;
            ax = axis;
            cursorHandles(col)=line([cursorLatency cursorLatency], [ax(3) ax(4)],'color','black','linestyle','--');
                
        end % next column
               
        linkaxes(plotAxes,'x');  % autoscales x axis between plots
        
        hold off;
    end

    % get processed trial data...
    function [timebase, fd] = getTrial(channel, startTime)
              
        % check trial boundaries

        if startTime < header.epochMinTime 
            startTime = header.epochMinTime;
        end 
        
        if startTime + epochTime >= header.epochMaxTime
            startTime = header.epochMaxTime-epochTime;
        end
        
        epochStart = startTime;
        
        % get data - note data indices start at 1;
        
        startSample = round( (epochStart - header.epochMinTime) * header.sampleRate) + 1;
        if startSample < 1
            startSample = 1;
        end
        endSample = startSample + epochSamples;
        if endSample > header.numSamples
            endSample = header.numSamples;
        end
        
        fd = dataarray(channel, startSample:endSample);
        timebase = timeVec(startSample:endSample);

    end

    function clickedOnLabel(src,~)

        % do shift select
        modifier = get(gcf,'CurrentModifier');
        idx = get(src,'UserData');
                
        if strcmp(modifier,'shift')
            sel = find(selectedMask == 1);
            if ~isempty(selectedChannelList)
                if idx > sel(1)
                    selectedMask(sel:idx) = 1;
                else
                    selectedMask(idx:sel) = 1;
                end
            else
                selectedMask(idx) = 1;
            end                
        elseif strcmp(modifier,'command')
            selectedMask(idx) = ~selectedMask(idx);       
        elseif strcmp(modifier,'control')
            % if right contextual mouse click always is selected
            selectedMask(idx) = 1;
        else 
            val = selectedMask(idx);
            selectedMask(:) = 0;
            selectedMask(idx) = ~val;           
        end

        drawTrial;
    end

function  capturekeystroke(~,evt)          

    % capture keypresses here
    if contains(evt.EventName,'KeyPress')       
        if contains(evt.Key,'a')
            if contains(evt.Modifier,'command')
                selectedMask(:) = 1;
                drawTrial;
            end
        end
        if contains(evt.Key,'rightarrow')
            dwel = 1.0 / header.sampleRate;
            t = cursorLatency + dwel;
            if t <= epochStart + epochTime
                cursorLatency = t;
            end
            updateCursors;
        end
        if contains(evt.Key,'leftarrow')
            dwel = 1.0 / header.sampleRate;
            t = cursorLatency - dwel;
            if t >= epochStart
                cursorLatency = t;
            end
            updateCursors;
        end
    end    
end


    % version 4.0 - new cursor function
    
    function updateCursors           
        if ~isempty(cursorHandles)
            set(cursorHandles(:), 'XData', [cursorLatency cursorLatency]);      
        end 
        
        s = sprintf('Latency: %.4f s', cursorLatency);
        set(cursor_text, 'string', s);

        if ~isempty(amplitudeLabels) && showAmplitudes                      
                for k=1:length(selectedChannelList)
                    if ~overlayPlots && ~enableMarking
                       [~,fd] = getTrial(selectedChannelList(k),epochStart);
                        % get sample offset from beginning of trial
                        offset = cursorLatency - epochStart;
                        sample = round( offset * header.sampleRate) + 1;        
                        val = fd(sample);
                        units = char(amplitudeUnits(k));
                        s = sprintf('%.2g %s', val,units);
                        set(amplitudeLabels(k),'string',s)               
                    else
                        set(amplitudeLabels(k),'string','')
                    end
                end
        end
        updateMap;
    end
        
    function buttondown(~,~) 
        
        if isempty(cursorHandles) || isempty(header)
            return;
        end

        ax = gca;
        % get current latency in s (x coord)
        mousecoord = get(ax,'currentpoint');
        x = mousecoord(1,1);
        y = mousecoord(1,2);
        xbounds = xlim;
        ybounds = ylim;
        if x < xbounds(1) | x > xbounds(2) | y < ybounds(1) | y > ybounds(2)
            return;
        end
        
        % move to current location on click ...
        cursorLatency = mousecoord(1,1);
        % snap to nearest sample
        dwel = 1.0 / header.sampleRate;
        cursorLatency = dwel * floor(cursorLatency/dwel);
        updateCursors;        
        ax = gca;
        set(fh,'WindowButtonMotionFcn',{@dragCursor,ax}) % need to explicitly pass the axis handle to the motion callback   
    end

    % button down function - drag cursor
    function dragCursor(~,~, ax)
        mousecoord = get(ax,'currentpoint');
        x = mousecoord(1,1);
        y = mousecoord(1,2);
        xbounds = xlim;
        ybounds = ylim;
        if x < xbounds(1) | x > xbounds(2) | y < ybounds(1) | y > ybounds(2)
            return;
        end        
        cursorLatency = mousecoord(1,1);
         % snap to nearest sample
        dwel = 1.0 / header.sampleRate;
        cursorLatency = dwel * floor(cursorLatency/dwel);
       
        updateCursors;
    end

    % on button up event set motion event back to no callback 
    function stopdrag(~,~)
        set(fh,'WindowButtonMotionFcn','');
    end

    %%%%%%%%%%%%%%%%% search for events %%%%%%%%%%%%%%%

    function markData
      
        
        if ~isempty(eventList)
            hadEvents = true;
        else
            hadEvents = false;
        end
        
        eventList = [];
        numEvents = 0;
        
        excludedEvents = 0;
        
        if reverseScan
            fprintf('Searching for events in reverse direction...\n');
            sample = header.numSamples;
            inc = -1;
        else
            fprintf('searching for events...\n');
            sample = 1;
            inc = 1;
        end
        
        % now works on normalized data scale, avoids have to rescale when
        % switching channel types

        chanIdx = selectedChannelList(1);
        fd = dataarray(chanIdx,:);
        mx = max(fd);
        fd = fd ./ mx;

        while (1)          
            value = fd(sample);
            
            % scan until found suprathreshold value
            if value < threshold
                sample = sample + inc;
            else
                % mark event    
                eventSample = sample;
                latency = timeVec(eventSample);     
                includeEvent = true;
                            
                % scan data until value drops below threshold
                % this is the duration of the event. 
                while (value > threshold)
                    sample = sample + inc;
                    if sample >= header.numSamples || sample < 1
                        break;
                    end
                    value = fd(sample);
                end  
                
                if sample > header.numSamples || sample < 1
                    break;
                end
                
                if reverseScan
                    eventData = fd(sample:eventSample);
                else
                    eventData = fd(eventSample:sample);
                end
                
                % *** exclusion critera ***
                
                % check for mininum required duration of the event (time till drops below
                % threshold
                
                eventDuration = (length(eventData) - 1) / header.sampleRate;                
                if eventDuration < minDuration
                    includeEvent = false;
                end                 

                % NEW *** instead of only checking min amplitude (meaning the minimum 
                % amplitude event had to reach to be accepted .. now checks that peak is within
                %  min / max range. 
                % i.e., has to surpass threshold AND reach at least min amplitude but not exceed
                %  max ampltiude.
                
                peak = max( eventData );
                
                % event must not exceed max amplitude
                if peak > maxAmplitude
                    includeEvent = false;
                end
                
                % event must reach min amplitude
                 if peak < minAmplitude
                    includeEvent = false;
                end   
                                
                % check for min. separation.
                
                if length(eventList) > 1
                    previousLatency = eventList(end);
                    separation = abs(latency - previousLatency);
                    
                    if separation < minSeparation
                        includeEvent = false;
                    end             
                end
                
                
                
                % if meets criteria add this event
                if includeEvent
                    if isempty(eventList)
                        eventList = latency;
                    else
                        eventList = [eventList latency];
                    end
                else
                    excludedEvents = excludedEvents + 1;
                end
                
                sample = sample+inc;
                               
            end

            if sample > header.numSamples || sample < 1
                break;
            end
        
        end
  
        if isempty(eventList)
            fprintf('No events found...\n');
            numEvents = 0;
        else
        
            if reverseScan
                eventList = sort(eventList);
            end

            % go to first event found
            numEvents = length(eventList);
            if ~hadEvents
                currentEvent = 1;
                epochStart = eventList(currentEvent) - (epochTime / 2.0);
            end
            drawTrial;
            fprintf('Excluded %d event(s)...\n', excludedEvents);   
        end

        
        if numEvents > 0
            s = sprintf('Events: %d of %d', currentEvent, numEvents);
            set(numEventsTxt,'String',s);
            set(eventCtl(:),'enable','on');
        else
            s = sprintf('Events: %d of %d', 0, numEvents);
            set(numEventsTxt,'String',s);
            set(eventCtl(:),'enable','off');
        end
         
    end



   % need dialog to create conditional events ...
   
    function create_event_callback(~,~)

        latencies = bw_conditionalMarker(markerFileName);
        eventList = latencies;      
        numEvents = length(eventList);
      
        s = sprintf('Events: %d of %d', currentEvent, numEvents);
        set(numEventsTxt,'String',s);  
        currentEvent = 1;
        drawTrial;
    
    end


   % save latencies
    function load_events_callback(~,~)
        
        [loadlatname, loadlatpath, ~] = uigetfile( ...
            {'*.txt','Text File (*.txt)'}, ...
               'Select Event File', dsName);
        
        loadlatfull=[loadlatpath,loadlatname];
        if isequal(loadlatname,0) 
           return;
        end
  
        new_list = importdata(loadlatfull);
        
        new_list = new_list';
        
        if ~isempty(eventList)          
           response = questdlg('Replace or add to current events?','Event Marker','Replace','Add','Cancel','Replace');
           if strcmp(response,'Cancel')
               return;
           end
           if strcmp(response,'Replace')
               eventList = new_list;           
           else
               eventList = [eventList new_list];
               eventList = sort(eventList);
           end
        else        
            eventList = new_list;
        end
        
        numEvents = length(eventList);
        
        currentEvent = 1;
        drawTrial;
        
        set(eventCtl(:),'enable','on');
        s = sprintf('Events: %d of %d', currentEvent, numEvents);
        set(numEventsTxt,'String',s);
        
    end

    function import_fif_data_callback(~,~)

       % check if we can convert fiff files...
        if ismac || ispc
            warndlg('Linux OS required to run fiff2ctf conversion program');
            return;
        end
        
        fileList = uigetdir2(pwd,'Select dataset(s) to import...');  
        if isempty(fileList)
            return;
        end

        numFiles = size(fileList,2);

        s = sprintf('Convert %d datasets to CTF format?', numFiles);           
        response = questdlg(s,'BrainWave','Yes','No','Yes');
        if strcmp(response,'No')    
            return;
        end

        wbh = waitbar(0,'Converting datasets...');

        datafile = [];

        for j=1:numFiles
            datafile = char(fileList{j});
            [loadpath, name, ext] = bw_fileparts(datafile);
            
            fiffPath = sprintf('%s%s%s%s%s',BW_PATH,'external',filesep,'linux',filesep);

            s = sprintf('Converting Elekta-Neuromag dataset %d of %d', j, numFiles);
            waitbar(j/numFiles,wbh,s);

            dsName = strrep(strcat(name,ext),'.fif','.ds');           
            tempDir = strrep(datafile, '.fif','_tempDir');
            tempFile = sprintf('%s%s%s', tempDir,filesep,dsName);

            cmd = sprintf('%sfiff2ctf %s %s',fiffPath,datafile,tempDir);
            system(cmd);

            if ~exist(tempFile,'dir')
                fprintf('Cannot find <%s> Conversion may have failed...', tempFile) 
                return;
            end

            cmd = sprintf('mv %s %s', tempFile, loadpath);
            system(cmd);
            cmd = sprintf('rmdir %s', tempDir);
            system(cmd);

            datafile = strrep(datafile,'.fif','.ds');
            newDsList(j) = cellstr(datafile);

            % from Paul Ferrari
            % create MarkerFile here from the Neuromag STIM channels
            trig = bw_getNMTriggers(datafile);
            bw_write_MarkerFile(datafile,trig);                
        end

        delete(wbh);

        dsName = datafile; % load last file converted ..?
        initData;
        drawTrial;

    end

    function import_kit_data_callback(~,~)

        % v. 4.1 - use multiple file select dialog to avoid confusion due to
        % different KIT file naming conventions.

        [conFile, markerFile, evtFile] = import_KIT_files;

        if isempty(conFile) || isempty(markerFile)
            errordlg('Insufficient data files specified...');
            return;
        end

        wbh = waitbar(0,'Converting dataset...');
        s = sprintf('Converting Yokagawa-KIT dataset...');   
        waitbar(0.5,wbh, s);
        
        success = con2ctf(conFile, markerFile, evtFile);
        delete(wbh);
        if success == -1 
            errordlg('KIT conversion failed...');
            return;
        end

        dsName = strrep(conFile,'.con','.ds');      
        
        initData;
        drawTrial;
        
    end


    function [conFile, markerFile, eventFile] = import_KIT_files

        conFile = [];
        markerFile = [];
        eventFile = [];

        d = figure('Position',[500 800 1000 300],'Name','Import KIT Data', ...
            'numberTitle','off','menubar','none');

        if ispc
            movegui(d,'center')
        end

        uicontrol('Style','text',...
            'fontsize',12,...
            'HorizontalAlignment','left',...
            'units', 'normalized',...
            'Position',[0.02 0.77 0.4 0.1],...
            'String','Select KIT continuous data file (e.g., 123_3_09_2022_B1.con)');

        con_edit = uicontrol('Style','edit',...
            'fontsize',10,...
            'units', 'normalized',...
            'HorizontalAlignment','left',...
            'Position',[0.02 0.7 0.75 0.1],...
            'String','');

        uicontrol('Style','pushbutton',...
            'fontsize',12,...
            'units', 'normalized',...
            'Position',[0.8 0.7 0.15 0.1],...
            'String','Browse',...
            'Callback',@select_con_callback);  

        uicontrol('Style','text',...
            'fontsize',12,...
            'HorizontalAlignment','left',...
            'units', 'normalized',...
            'Position',[0.02 0.57 0.4 0.1],...
            'String','Select Marker file for co-registration (e.g., 123_3_09_2022_ini.mrk):');

        marker_edit = uicontrol('Style','edit',...
            'fontsize',10,...
            'units', 'normalized',...
            'HorizontalAlignment','left',...
            'Position',[0.02 0.5 0.75 0.1],...
            'String','');
        uicontrol('Style','pushbutton',...
            'fontsize',12,...
            'units', 'normalized',...
            'Position',[0.8 0.5 0.15 0.1],...
            'String','Browse',...
            'Callback',@select_marker_callback);              

        uicontrol('Style','text',...
            'fontsize',12,...
            'HorizontalAlignment','left',...
            'units', 'normalized',...
            'Position',[0.02 0.37 0.4 0.1],...
            'String','Optional: Select event file (.evt) containing trigger events:');

        event_edit = uicontrol('Style','edit',...
            'fontsize',10,...
            'units', 'normalized',...
            'HorizontalAlignment','left',...
            'Position',[0.02 0.3 0.75 0.1],...
            'String','');    

        uicontrol('Style','pushbutton',...
            'fontsize',12,...
            'units', 'normalized',...
            'Position',[0.8 0.3 0.15 0.1],...
            'String','Browse',...
            'Callback',@select_event_callback);          

        uicontrol('Style','pushbutton',...
            'fontsize',12,...
            'foregroundColor','blue',...
            'units', 'normalized',...
            'Position',[0.75 0.1 0.2 0.1],...
            'String','Import',...
            'Callback',@OK_callback);  

        uicontrol('Style','pushbutton',...
            'fontsize',12,...
            'units', 'normalized',...
            'Position',[0.5 0.1 0.2 0.1],...
            'String','Cancel',...
            'Callback','delete(gcf)') 

        function OK_callback(~,~)
            % get from text box in case user typed in
            conFile = get(con_edit,'string');
            markerFile = get(marker_edit,'string');        
            eventFile = get(event_edit,'string');

            delete(gcf)
        end

        function select_con_callback(~,~)
            s =uigetfile('*.con','Select KIT .con file ...');
            if isequal(s,0)
                return;
            end    

            set(con_edit,'string',s);      

         end

        function select_marker_callback(~,~)
            s =uigetfile('*.mrk','Select KIT .con file ...');
            if isequal(s,0)
                return;
            end    
            set(marker_edit,'string',s);
        end

        function select_event_callback(~,~)    
            s =uigetfile('*.evt','Select KIT .con file ...');
            if isequal(s,0)
                return;
            end 
            set(event_edit,'string',s);      
        end

        % make modal   


        uiwait(d);

    end

    function load_KIT_events_callback(~,~)
        
        [loadlatname, loadlatpath, ~] = uigetfile( ...
            {'*.evt','KIT Event File (*.evt)'}, ...
               'Select Event File', dsName);
        
        loadlatfull=[loadlatpath,loadlatname];
        if isequal(loadlatname,0) 
           return;
        end
  
        [new_list, ~] = bw_readMACCSEventFile(loadlatfull);
       
        new_list = new_list';
        
        if ~isempty(eventList)          
           response = questdlg('Replace or add to current events?','Event Marker','Replace','Add','Cancel','Replace');
           if strcmp(response','Cancel')
               return;
           end
           if strcmp(response','Replace')
               eventList = new_list;           
           else
               eventList = [eventList new_list];
               eventList = sort(eventList);
           end
        else        
            eventList = new_list;
        end
        
        numEvents = length(eventList);
        
        currentEvent = length(eventList);
        drawTrial;

        set(eventCtl(:),'enable','on');
        s = sprintf('Events: %d of %d', currentEvent, numEvents);
        set(numEventsTxt,'String',s);
        
    end

    function load_marker_events_callback(~,~)
                  
        if ~exist(markerFileName,'file')
            errordlg('No marker file exists yet. Create or import latencies then save events as markers.');
            return;
        end
        
        [new_list, ~] = bw_readCTFMarkers(markerFileName);       
        new_list = new_list';
        
        if ~isempty(eventList)          
           response = questdlg('Replace or add to current events?','Event Marker','Replace','Add','Cancel','Replace');
           if strcmp(response','Cancel')
               return;
           end
           if strcmp(response,'Replace')
               eventList = new_list;           
           else
               eventList = [eventList new_list];
               eventList = sort(eventList);
           end
        else        
            eventList = new_list;
        end
        
        numEvents = length(eventList);
        
        currentEvent = length(eventList);
        drawTrial;

        set(eventCtl(:),'enable','on')
        
        s = sprintf('Events: %d of %d', currentEvent, numEvents);
        set(numEventsTxt,'String',s);
        
        % first_event_callback;
        
    end

    % save latencies
    function save_events_callback(~,~)
        if isempty(eventList)
            errordlg('No events defined ...');
            return;
        end

        saveName = strcat(dsName, filesep, '*.txt');
        [name,path,~] = uiputfile('*.txt','Save Event latencies in File:',saveName);
        if isequal(name,0)
            return;
        end         

        eventFile = fullfile(path,name);
        fprintf('Saving event times to text file %s \n', eventFile);

        fid = fopen(eventFile, 'w');
        fprintf(fid,'%.5f\n', eventList');
        fclose(fid);
    end

    % save current events as a marker 
    function save_marker_callback(~,~)
        if isempty(eventList)
            errordlg('No events defined ...');
            return;
        end 
        
        newName = getMarkerName(markerNames);
        
        if isempty(newName)
            return;
        end
            
        if ~isempty( find( strcmp(newName, markerNames) == 1))
            beep;
            warndlg('A Marker with this name already exists for this dataset...');
            return;
        end
        
        % add current event as marker               
        if numMarkers > 0        
            for i=1:numMarkers
                trig(i).ch_name = char(markerNames(i+1));
                markerTimes = markerLatencies{i};                
                trig(i).times = markerTimes;
            end         
        end     
        
        % add to MarkerFile
        trig(numMarkers+1).ch_name = newName;
        trig(numMarkers+1).times = eventList;
        success = bw_writeCTFMarkerFile(dsName, trig);  
        
        % if not cancelled reload new file      
        if success
            loadMarkerFile(markerFileName);
        end
        
    end

    % save current marker as ascii event file
    function save_marker_as_excel_callback(~,~)
        
        saveName = strrep(markerFileName,'.mrk','.csv');
        [name,path,~] = uiputfile('*.csv','Export Markers to File:',saveName);
        if isequal(name,0)
            return;
        end         
        
        eventFile = fullfile(path,name);
        fprintf('Saving marker data to csv file %s \n', eventFile);
        
        bw_markerFile2Excel(markerFileName, saveName)
        
    end


    function updateMap
        
        if ~showMap     
            return;
        end
        
        currentAx = gca;
        subplot('Position', mapbox);

        if isempty(mapLocs)
            mapLocs = initMap;
        end

        % get sample offset from beginning of trial   
        sample = round(cursorLatency * header.sampleRate) + header.numPreTrig + 1; % add one sample for t=0
        map_data = dataarray(meg_idx,sample);  % read processed data! 

        tempLocs = mapLocs;
        
        % don't pass null data
        if mean(map_data) ~= 0.0
            topoplot(map_data', tempLocs, 'colormap',jet,'numcontour',8,'electrodes','on','shrink',0.15);
        end
        axes(currentAx);
%         subplot('Position', plotbox(1));

    end

    function map_locs = initMap

%         % create an EEGLAB chanlocs structure to avoid having to save .locs file
        map_locs = struct('labels',{},'theta', {}, 'radius', {});

        channelIndex = 1;
        for i=1:header.numChannels

            chan = header.channel(i);
            if ~chan.isSensor
                continue;
            end

            name = chan.name;
            % remove dashes in CTF names
            idx = strfind(name,'-');
            if ~isempty(idx)
                temp=name(1:idx-1);
                name = temp;
            end

            X = chan.xpos;
            Y = chan.ypos;
            Z = chan.zpos;
            

            [th, phi, ~] = cart2sph(X,Y,Z);

            decl = (pi/2) - phi;
            radius = decl / pi;
            theta = th * (180/pi);
            if (theta < 180)
                theta = -theta;
            else
                theta = 360 - theta;
            end

            map_locs(channelIndex).labels = name;
            map_locs(channelIndex).theta = theta;
            map_locs(channelIndex).radius = radius;  
            channelIndex = channelIndex + 1;
        end
    end


    % if ~exist('dsName','var')
    %     dsName = uigetdir('.ds', 'Select CTF dataset ...');
    %     if dsName == 0
    %         return;
    %     end    
    % end
    % 
    % % draw after controls defined.
    % initData;
    % drawTrial;
    % 
end

%%% helper functions...

% strip sensor version number from channel names for CTF
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

function [markerName] = getMarkerName(existingNames)
 
    markerName = 'newMarker';

    scrnsizes=get(0,'MonitorPosition');
    fg=figure('color','white','name','Choose Marker Name','numbertitle','off','menubar','none','position',[300 (scrnsizes(1,4)-300) 500 150]);
    

    markerNameEdit = uicontrol('style','edit','units','normalized','HorizontalAlignment','Left',...
        'position',[0.2 0.6 0.6 0.2],'String',markerName,'Backgroundcolor','white','fontsize',12);
                
    % check for duplicate name
    if ~isempty(existingNames)
        
    end      

    uicontrol('style','pushbutton','units','normalized','position',...
        [0.6 0.3 0.25 0.25],'string','OK','backgroundcolor','white',...
        'foregroundcolor','blue','callback',@ok_callback);
    
    function ok_callback(~,~)
        markerName = get(markerNameEdit,'string');
        uiresume(gcf);
    end

    uicontrol('style','pushbutton','units','normalized','position',...
        [0.2 0.3 0.25 0.25],'string','Cancel','backgroundcolor','white',...
        'foregroundcolor','black','callback',@cancel_callback);
    
    function cancel_callback(~,~)
        markerName = [];
        uiresume(gcf);
    end
    
    
    %%PAUSES MATLAB
    uiwait(gcf);
    %%CLOSES GUI
    close(fg);   
    
    
end



