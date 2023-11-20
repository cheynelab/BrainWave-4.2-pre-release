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
currentEvent = 1;
numEvents = 0;

threshold = 0;
maxAmplitude = 0.0;
minAmplitude = 0.0;
minSeparation = 0.0;


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

cursorHandle = 0;
cursorLatency = 0.0;

% conditional marking
showMarkerWindow = 0;
markerWindowStart = -0.1;
markerWindowEnd = 0.1;

overlayPlots = 0;

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

fh = figure('numbertitle','off','position',[200, 800, 2000, 1500],...
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
loadMenu = uimenu(markerMenu,'label','Load Events from ...');
uimenu(loadMenu,'label','CTF MarkerFile','callback',@load_marker_events_callback)
uimenu(loadMenu,'label','KIT Event File','callback',@load_KIT_events_callback)
uimenu(loadMenu,'label','Text File','callback',@load_events_callback)

uimenu(markerMenu,'label','Create Conditional Event...','callback',@create_event_callback)
uimenu(markerMenu,'label','Edit Marker File...','separator','on','callback',@edit_markers_callback)
uimenu(markerMenu,'label','Export MarkerFile to Excel...','callback',@save_marker_as_excel_callback)

% uimenu(markerMenu,'label','Export Markers to Text File...','callback',@save_events_callback)


% build channel menu once
for kk=1:numel(channelSets)
    s = char(channelSets(kk));
    channelMenuItems(kk) = uimenu(channelMenu,'label',s,'callback',@channel_menu_callback); 
end
channelMenuItems(end+1) = uimenu(channelMenu,'label','Edit Custom',...
        'separator','on','callback',@editChannelSet_callback); 


% +++++++++++++ set plot window +++++++++++++
subplot('position',[0.05 0.31 0.9 0.65]);
    

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
    answer=inputdlg({'Start time (s):','End Time (s)','Downsample Factor:',...
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
    
    err = bw_CTFNewDs(dsName, newDsName, filterFlag, bandPass, badChans, badTrials, sampleRange, ds, gradient);  
    if err ~= 0
        errordlg('bw_CTFNewDs returned error');
        return;
    end

    
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

    loadData;
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

    set(overlayPlotsCheck,'value',overlayPlots);

    updateChannelMenu; 
    
end

function saveDefaults(defaultsFile)
    params.channelMenuIndex = channelMenuIndex;
    params.selectedChannelList = selectedChannelList;
    params.overlayPlots = overlayPlots;
    params.epochTime = epochTime;
    
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
    currentEvent = 1;
    numEvents = 0;
    enableMarking = 0;

    threshold = 0;
    maxAmplitude = 0.0;
    minAmplitude = 0.0;
    minSeparation = 0.0;
    
    header = bw_CTFGetHeader(dsName);
    
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
    else
        set(trialIncButton,'enable','off');
        set(trialDecButton,'enable','off');
    end

    cursorLatency = header.epochMinTime + (epochTime/2);
    epochSamples = round(epochTime * header.sampleRate);

    % when loading new dataset reset bandpass and turn filter off
    bandPass(1) = 1;
    bandPass(2) = header.sampleRate / 2.0;
    set(filt_hi_pass,'string',bandPass(1),'enable','off')  
    set(filt_low_pass,'string',bandPass(2),'enable','off')
    filterOff = 1;
    set(filterCheckBox,'value',0);
    
    fprintf('Loading data...\n\n');
    
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
    
    channelMenuIndex = 1;

    % set default display on opening to first MEG sensor...
    idx = find(channelTypes == 5);
    if ~isempty(idx)
        selectedChannelList = idx(1);
    else
        selectedChannelList = 1;
    end
    
    customChannelList = selectedChannelList;
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

    loadData;   
    drawTrial;

    updateChannelMenu;
    updateMarkerControls;
    updateSlider;
    updateCursors;

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
        
    for k=1:numel(channelSets)
        % turn off default channel types that don't exist
        switch k 
            case 2
                if ~ismember(channelTypes,MEG_CHANNELS)
                    set(channelMenuItems(k),'enable','off');
                end
            case 3
                if ~ismember(channelTypes,ADC_CHANNELS)
                    set(channelMenuItems(k),'enable','off');
                end    
            case 4
                if ~ismember(channelTypes,TRIGGER_CHANNELS)
                    set(channelMenuItems(k),'enable','off');
                end                     
            case 5
                if ~ismember(channelTypes,DIGITAL_CHANNELS)
                    set(channelMenuItems(k),'enable','off');
                end         
            case 6
                if ~ismember(channelTypes,EEG_CHANNELS)
                    set(channelMenuItems(k),'enable','off');
                end
        end        
        if channelMenuIndex == k
            set(channelMenuItems(k),'checked','on')
        else
            set(channelMenuItems(k),'checked','off')
        end
    end

end

function channel_menu_callback(src,~)

    if isempty(channelNames)
        return;
    end
    
    channelMenuIndex = get(src,'position');
    menuItems = get(channelMenu,'Children');

    nchans = numel(channelNames);
    channelExcludeFlags = ones(1,nchans);             

    switch channelMenuIndex 
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
        
    set(get(channelMenu,'Children'),'Checked','off');
    if channelMenuIndex < 7
        set(src,'Checked','on')
    end
     
    updateMarkerControls;

    loadData; 

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
    % get handles to menu items - indices always in reverse order
    channelMenuItems = get(channelMenu,'Children');
    set(channelMenuItems(:),'Checked', 'off')
    set(channelMenuItems(end), 'Checked','on')  

    loadData;
    drawTrial;
    autoScale_callback;
    updateMarkerControls;

end

%+++++++++ event detection controls +++++++++  

annotation('rectangle',[0.6 0.02 0.35 0.18],'EdgeColor','blue');

uicontrol('style','text','units','normalized','position',[0.63 0.18 0.12 0.025],...
    'string','Mark Events (Single Channel)','backgroundcolor','white','foregroundcolor','blue','fontweight','bold',...
    'FontSize',11);
   

function updateMarkerControls
    if numel(selectedChannelList) == 1 && header.numTrials == 1
        set(enable_marking_check,'enable','on');
        setMarkingCtls('on');
    else
        setMarkingCtls('off');
        set(enable_marking_check,'enable','off');
        enableMarking = 0; 
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
    set(add_event_button,'enable',state)  
    set(delete_event_button,'enable',state)  
    set(delete_all_button,'enable',state)  
    set(save_events_button,'enable',state)  
    set(eventCtl(:),'enable',state)  
    set(numEventsTxt,'enable',state)    

end

threshold_text = uicontrol('style','text','units','normalized','position',[0.62 0.135 0.2 0.02],...
    'enable','off','string','Threshold:','fontsize',11,'backgroundcolor','white','horizontalalignment','left');
threshold_edit=uicontrol('style','edit','units','normalized','position',[0.69 0.14 0.035 0.02],...
    'enable','off','FontSize', 11, 'BackGroundColor','white','string',threshold,'callback',@threshold_callback);
    function threshold_callback(src,~)
        threshold =str2double(get(src,'string'));
        drawTrial;
    end
min_duration_text = uicontrol('style','text','units','normalized','position',[0.62 0.1 0.12 0.02],...
    'enable','off','string','Min. duration (s):','fontsize',11,'backgroundcolor','white','horizontalalignment','left');
min_duration_edit=uicontrol('style','edit','units','normalized','position',[0.69 0.105 0.035 0.02],...
    'enable','off','FontSize', 11, 'BackGroundColor','white','string',minDuration,...
    'callback',@min_duration_callback);

    function min_duration_callback(src,~)
        minDuration=str2double(get(src,'string'));
    end

min_amp_text = uicontrol('style','text','units','normalized','position',[0.62 0.06 0.05 0.03],...
    'enable','off','string','Amplitude Range:','fontsize',11,'backgroundcolor','white','horizontalalignment','left');

min_amplitude_edit=uicontrol('style','edit','units','normalized','position',[0.69 0.07 0.035 0.02],...
    'enable','off','FontSize', 11, 'BackGroundColor','white','string',minAmplitude,...
    'callback',@min_amplitude_callback);
    function min_amplitude_callback(src,~)
        minAmplitude=str2double(get(src,'string'));
        drawTrial;
    end

max_amplitude_edit=uicontrol('style','edit','units','normalized','position',[0.75 0.07 0.035 0.02],...
    'enable','off','FontSize', 11, 'BackGroundColor','white','string',maxAmplitude,...
    'callback',@max_amplitude_callback);
min_amp_text2 = uicontrol('style','text','units','normalized','position',[0.735 0.07 0.01 0.02],...
    'enable','off','string','to','fontsize',11,'backgroundcolor','white','horizontalalignment','left');
    function max_amplitude_callback(src,~)
        maxAmplitude=str2double(get(src,'string'));
        drawTrial;
    end

min_sep_text = uicontrol('style','text','units','normalized','position',[0.62 0.03 0.12 0.02],...
    'enable','off','string','Min. separation (s):','fontsize',11,'backgroundcolor','white','horizontalalignment','left');
min_sep_edit = uicontrol('style','edit','units','normalized','position',[0.69 0.0345 0.035 0.02],...
    'enable','off','FontSize', 11, 'BackGroundColor','white','string',minSeparation,...
    'callback',@min_separation_callback);

    function min_separation_callback(src,~)
        minSeparation=str2double(get(src,'string'));
        drawTrial;
    end

forward_scan_radio = uicontrol('style','radiobutton','units','normalized','position',[0.62 0.165 0.08 0.02],...
    'enable','off','string','rising edge','backgroundcolor','white','value',~reverseScan,'FontSize',11,'callback',@forward_check_callback);
reverse_scan_radio = uicontrol('style','radiobutton','units','normalized','position',[0.68 0.165 0.08 0.02],...
    'enable','off','string','falling edge','backgroundcolor','white','value',reverseScan,'FontSize',11,'callback',@reverse_check_callback);

find_events_button = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.86 0.16 0.06 0.025],...
    'enable','off','string','Find Events','Foregroundcolor','blue','backgroundcolor','white','callback',@find_events_callback);

add_event_button = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.78 0.16 0.06 0.025],...
    'enable','off','string','Insert Event','Foregroundcolor','blue','backgroundcolor','white','callback',@add_event_callback);

delete_event_button = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.78 0.12 0.06 0.025],...
    'enable','off','string','Delete Event','Foregroundcolor','blue','backgroundcolor','white','callback',@delete_event_callback);

delete_all_button = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.86 0.12 0.06 0.025],...
    'enable','off','string','Clear Events','Foregroundcolor','blue','backgroundcolor','white','callback',@delete_all_callback);

save_events_button = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.86 0.08 0.06 0.025],...
    'enable','off','string','Save Events','Foregroundcolor','blue','backgroundcolor','white','callback',@save_marker_callback);

enable_marking_check = uicontrol('style','checkbox','units','normalized','fontsize',11,'position',[0.75 0.03 0.05 0.025],...
    'enable','off','string','Enable','Foregroundcolor','blue','backgroundcolor','white','callback',@enable_marking_callback);


eventCtl(1) = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.81 0.03 0.025 0.02],...
    'enable','off','CData',leftarrow_im,'Foregroundcolor','black','backgroundcolor','white','callback',@event_dec_callback);
eventCtl(2) = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.84 0.03 0.025 0.02],...
    'enable','off', 'CData', rightarrow_im,'Foregroundcolor','black','backgroundcolor','white','callback',@event_inc_callback);

s = sprintf('# Events = %d', numEvents);
numEventsTxt = uicontrol('style','text','units','normalized','position',[0.88 0.03 0.06 0.02],...
    'string',s,'fontsize',12,'fontweight','bold','backgroundcolor','white', 'foregroundcolor', 'red','horizontalalignment','left');

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
        end
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

        s = sprintf('# events = %d', numEvents);
        set(numEventsTxt,'String',s);

        drawTrial;
    end

    function delete_event_callback(~,~)

       if numEvents < 1
           errordlg('No events to delete ...');
           return;
       end
       % make sure we have currentEvent in window
       if ( eventList(currentEvent) > epochStart & eventList(currentEvent) < (epochStart+epochTime) )        
           s = sprintf('Delete event #%d? (Cannot be undone)', currentEvent);
           response = questdlg(s,'Mark Events','Yes','No','No');
           if strcmp(response,'No')
               return;
           end
           eventList(currentEvent) = [];
           numEvents = length(eventList);         
           drawTrial;
           s = sprintf('# events = %d', numEvents);
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
        drawTrial;
        s = sprintf('# events = %d', numEvents);
        set(numEventsTxt,'string',s)
    end

    % now works on normalized scale
    function enable_marking_callback(~,~)
        
        enableMarking = ~enableMarking;
      
        if enableMarking
            % reset params for marking        

            maxAmplitude = 0.8;
            minAmplitude = 0.0;
            threshold = 0.1;

            s = sprintf('%.2g', maxAmplitude);
            set(max_amplitude_edit,'string',s);
            s = sprintf('%.2g', minAmplitude);
            set(min_amplitude_edit,'string',s);
            s = sprintf('%.2g', threshold);
            set(threshold_edit,'string',s);

            eventList = [];
            numEvents = 0;

            s = sprintf('# events = %d', numEvents);
            set(numEventsTxt,'string',s)
                       
            % make sure data is plotted showing normalized range and
            % disable scaling


            setMarkingCtls('on');
                
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
    'sliderStep', [1 100],'BackGroundColor',[0.8 0.8 0.8],'ForeGroundColor',...
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

uicontrol('style','popupmenu','units','normalized','fontsize',11,'position',[0.05 0.21 0.08 0.03],...
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

uicontrol('style','text','units','normalized','position',[0.14 0.22 0.1 0.02],...
    'string','Max:','fontsize',11,'backgroundcolor','white','horizontalalignment','left');
max_scale=uicontrol('style','edit','units','normalized','position',[0.16 0.22 0.04 0.025],...
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

uicontrol('style','text','units','normalized','position',[0.22 0.22 0.08 0.02],...
    'string','Min:','fontsize',11,'backgroundcolor','white','horizontalalignment','left');
min_scale=uicontrol('style','edit','units','normalized','position',[0.24 0.22 0.04 0.025],...
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

    function scaleUp_callback(~,~)
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


     
scaleUpArrow = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.3 0.22 0.03 0.025],...
    'CData',uparrow_im,'Foregroundcolor','black','backgroundcolor','white','callback',@scaleUp_callback);

scaleDownArrow = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.34 0.22 0.03 0.025],...
    'CData',downarrow_im,'Foregroundcolor','black','backgroundcolor','white','callback',@scaleDown_callback);


autoScaleButton = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.38 0.22 0.06 0.025],...
    'Foregroundcolor','black','string','Autoscale','backgroundcolor','white','callback',@autoScale_callback);
    function autoScale_callback(~,~)
        maxRange(currentScaleMenuIndex) = NaN;
        minRange(currentScaleMenuIndex) = NaN;
        drawTrial;
    end


trialNumTxt = uicontrol('style','text','fontsize',12,'units','normalized','horizontalalignment','left','position',...
     [0.46 0.21 0.08 0.03],'string','Trial: 1 of 1','BackgroundColor','white','foregroundcolor','black');

trialDecButton = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.515 0.22 0.025 0.025],...
    'CData',leftarrow_im,'Foregroundcolor','black','backgroundcolor','white','callback',@trial_dec_callback);

trialIncButton = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.545 0.22 0.025 0.025],...
    'CData',rightarrow_im,'Foregroundcolor','black','backgroundcolor','white','callback',@trial_inc_callback);


    function trial_inc_callback(~,~)  
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



function edit_markers_callback(~,~)    
    bw_editCTFMarkers(dsName);
    loadMarkerFile(markerFileName);
end

markerNames = {'none'};
uicontrol('style','text','fontsize',12,'units','normalized','horizontalalignment','left','position',...
     [0.6 0.21 0.08 0.03],'string','Markers:','BackgroundColor','white','foregroundcolor','black');
marker_Popup =uicontrol('style','popup','units','normalized','fontsize',12,'position',[0.64 0.19 0.1 0.05],...
    'string',markerNames,'value',currentMarkerIndex,'Foregroundcolor','black','backgroundcolor','white','callback',@marker_popup_callback);

    function marker_popup_callback(src,~)    
        currentMarkerIndex = get(src,'value');
        if currentMarkerIndex == 1
            set(markerIncButton,'enable','off')
            set(markerDecButton,'enable','off')
        else
            set(markerIncButton,'enable','on')
            set(markerDecButton,'enable','on')
        end
        drawTrial;
    end


markerDecButton = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.75 0.22 0.025 0.025],...
    'CData',leftarrow_im,'Foregroundcolor','black','backgroundcolor','white','callback',@marker_dec_callback);
markerIncButton = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.78 0.22 0.025 0.025],...
    'CData',rightarrow_im,'Foregroundcolor','black','backgroundcolor','white','callback',@marker_inc_callback);

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



uicontrol('style','text','units','normalized','position',[0.825 0.22 0.1 0.02],...
    'string','Window Duration (s):','fontsize',11,'backgroundcolor','white','horizontalalignment','left');
epochDurationEdit = uicontrol('style','edit','units','normalized','position',[0.9 0.225 0.05 0.02],...
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


    function updateSlider
        epochSamples = round(epochTime * header.sampleRate);
        dataRange = epochTime / header.trialDuration * 2;
        sliderScale =[(dataRange * 0.05) dataRange];
        % sliderScale = [(epochTime / header.trialDuration) * 0.02 (epochTime / header.trialDuration) * 0.08]
        set(latency_slider,'sliderStep',sliderScale);
        drawTrial;
    end

% +++++++++++++++++++++++++

% ++++++++++ plot settings 

annotation('rectangle',[0.05 0.02 0.53 0.18],'EdgeColor','blue');
uicontrol('style','text','fontsize',11,'units','normalized','position',...
     [0.08 0.18 0.1 0.025],'string','Plot Settings','BackgroundColor','white','foregroundcolor','blue','fontweight','b');

sampleRateTxt = uicontrol('style','text','units','normalized','position',[0.08 0.16 0.08 0.02],...
    'string','Sample Rate:','backgroundcolor','white','FontSize',11, 'HorizontalAlignment','Left');

totalSamplesTxt = uicontrol('style','text','units','normalized','position',[0.16 0.16 0.08 0.02],...
    'string','Total Samples:','backgroundcolor','white','FontSize',11, 'HorizontalAlignment','Left');

totalDurationTxt = uicontrol('style','text','units','normalized','position',[0.24 0.16 0.08 0.02],...
    'string','Trial Duration:','backgroundcolor','white','FontSize',11, 'HorizontalAlignment','Left');

numTrialsTxt = uicontrol('style','text','units','normalized','position',[0.32 0.16 0.08 0.02],...
    'string','# of Trials:','backgroundcolor','white','FontSize',11, 'HorizontalAlignment','Left');


numChannelsTxt = uicontrol('style','text','units','normalized','position',[0.08 0.14 0.08 0.02],...
    'string','Total Channels:','backgroundcolor','white','FontSize',11, 'HorizontalAlignment','Left');

numSensorsTxt = uicontrol('style','text','units','normalized','position',[0.16 0.14 0.08 0.02],...
    'string','MEG Sensors:','backgroundcolor','white','FontSize',11, 'HorizontalAlignment','Left');

numReferencesTxt = uicontrol('style','text','units','normalized','position',[0.24 0.14 0.08 0.02],...
    'string','MEG References:','backgroundcolor','white','FontSize',11, 'HorizontalAlignment','Left');

numAnalogTxt = uicontrol('style','text','units','normalized','position',[0.32 0.14 0.08 0.02],...
    'string','ADC Channels:','backgroundcolor','white','FontSize',11, 'HorizontalAlignment','Left');



filterCheckBox = uicontrol('style','checkbox','units','normalized','position',[0.1 0.1 0.04 0.02],...
    'string','Filter','backgroundcolor','white','value',~filterOff,'FontSize',11,'callback',@filter_check_callback);

% replace with 50 / 60 Hz checks - check code 
uicontrol('style','checkbox','units','normalized','position',[0.1 0.07 0.06 0.02],...
    'string','60 Hz Notch','backgroundcolor','white','value',notchFilter,'FontSize',11,'callback',@notch_check_callback);
uicontrol('style','checkbox','units','normalized','position',[0.18 0.07 0.06 0.02],...
    'string','50 Hz Notch','backgroundcolor','white','value',notchFilter2,'FontSize',11,'callback',@notch2_check_callback);

uicontrol('style','checkbox','units','normalized','position',[0.45 0.1 0.1 0.02],...
    'string','Remove Offset','backgroundcolor','white','value',removeOffset,'FontSize',11,'callback',@remove_offset_callback);

uicontrol('style','checkbox','units','normalized','position',[0.1 0.04 0.08 0.02],...
    'string','Invert','backgroundcolor','white','value',invertData,'FontSize',11,'callback',@invert_check_callback);

uicontrol('style','checkbox','units','normalized','position',[0.18 0.04 0.08 0.02],...
    'string','Rectify','backgroundcolor','white','value',rectify,'FontSize',11,'callback',@rectify_check_callback);

uicontrol('style','checkbox','units','normalized','position',[0.26 0.04 0.12 0.02],...
    'string','Differentiate','backgroundcolor','white','value',differentiate,'FontSize',11,'callback',@firstDiff_check_callback);

uicontrol('style','checkbox','units','normalized','position',[0.34 0.04 0.08 0.02],...
    'string','Envelope (hilbert)','backgroundcolor','white','value',envelope,'FontSize',11,'callback',@envelope_check_callback);

overlayPlotsCheck = uicontrol('style','checkbox','units','normalized','position',[0.45 0.14 0.1 0.02],...
    'string','Overlay Plots','backgroundcolor','white','value',overlayPlots,'FontSize',11,'callback',@overlay_plots_callback);


uicontrol('style','text','units','normalized','position',[0.15 0.1 0.06 0.02],...
    'string','High Pass (Hz):','fontsize',11,'backgroundcolor','white','horizontalalignment','left');
filt_hi_pass=uicontrol('style','edit','units','normalized','position',[0.22 0.1 0.04 0.02],...
    'FontSize', 11, 'BackGroundColor','white','string',bandPass(1),...
    'callback',@filter_hipass_callback);

    function filter_hipass_callback(src,~)
        bandPass(1)=str2double(get(src,'string'));
        loadData;
        drawTrial;
    end

uicontrol('style','text','units','normalized','position',[0.28 0.1 0.06 0.02],...
    'string','Low Pass (Hz):','fontsize',11,'backgroundcolor','white','horizontalalignment','left');
filt_low_pass=uicontrol('style','edit','units','normalized','position',[0.35 0.1 0.04 0.02],...
    'FontSize', 11, 'BackGroundColor','white','string',bandPass(2),...
    'callback',@filter_lowpass_callback);

    function filter_lowpass_callback(src,~)
        bandPass(2)=str2double(get(src,'string'));
        loadData;
        drawTrial;
    end

if filterOff
    set(filt_hi_pass, 'enable','off');
    set(filt_low_pass, 'enable','off');
else
    set(filt_hi_pass, 'enable','on');
    set(filt_low_pass, 'enable','on');
end



function rectify_check_callback(src,~)
    rectify=get(src,'value');
    loadData;
    drawTrial;
end

function envelope_check_callback(src,~)
    envelope=get(src,'value');
    loadData;
    drawTrial;
end

function invert_check_callback(src,~)
    invertData=get(src,'value');
    loadData;
    drawTrial;
end


function firstDiff_check_callback(src,~)
    differentiate=get(src,'value');
    loadData;
    drawTrial;
end

function notch_check_callback(src,~)
    notchFilter=get(src,'value');
    loadData;
    drawTrial;
end

function notch2_check_callback(src,~)
    notchFilter2=get(src,'value');
    loadData;
    drawTrial;
end

function remove_offset_callback(src,~)
    removeOffset=get(src,'value');
    loadData;
    drawTrial;
end


function filter_check_callback(src,~)
    filterOff=~get(src,'value');
    if filterOff
        bandPass = [0 header.sampleRate / 2.0];
        set(filt_hi_pass, 'enable','off');
        set(filt_low_pass, 'enable','off');
    else
        bandPass(1)=str2double(get(filt_hi_pass,'string'));
        bandPass(2)=str2double(get(filt_low_pass,'string'));
        set(filt_hi_pass, 'enable','on');
        set(filt_low_pass, 'enable','on');
    end
    
    loadData;
    drawTrial;
end

function overlay_plots_callback(src,~)
    overlayPlots=get(src,'value');
    drawTrial;
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

% load (all) data with current filter settings etc and adjust scale
    
function loadData
    
    % reload all data
    numChannelsToDisplay = numel(selectedChannelList);
    dataarray = zeros(numChannelsToDisplay, header.numSamples);

    % if displaying whole head 
    if numChannelsToDisplay > 64
        showProg = 1;
    else
        showProg = 0;
    end

    if showProg
        wbh = waitbar(0,'Processing data...');
    end

    % transpose since bw_getCTFData returns [nsamples x nchannels]
    readAllChannels = 1;

        % if header.numTrials > 1
        %     data = tmp_data(:,trialNo);
        % else
        %     data = tmp_data;
        % end

    all_data = bw_getCTFData(dsName, 0, header.numSamples, readAllChannels)';  
    timeVec = header.epochMinTime: 1/ header.sampleRate: header.epochMaxTime;
    timeVec = timeVec(1:header.numSamples);  % should be same ...

    for k=1:numChannelsToDisplay

        if showProg
            s = sprintf('Processing channel %d of %d', k, numChannelsToDisplay);
            waitbar(k/numChannelsToDisplay,wbh,s);
        end

        chanIdx = selectedChannelList(k);
        channelName = channelNames{chanIdx};                         
        chanType = channelTypes(chanIdx);   

        aTypes = [MEG_CHANNELS ADC_CHANNELS EEG_CHANNELS];
        isAnalog = ismember(chanType,aTypes);

        data(1:header.numSamples) = all_data(chanIdx,1:header.numSamples);                        

        % filter data
        if ~filterOff
            trial = data;
            y = bw_filter(trial, header.sampleRate, bandPass); 
            data = y;                   
        end

                
        if notchFilter && isAnalog
            nyquist = header.sampleRate/2.0;
            d = data';
            data = bw_filter(d, header.sampleRate, [58 62], 4, 1, 1)';
            if nyquist > 120 
                d = data';
                data = bw_filter(d, header.sampleRate, [115 125], 4, 1, 1)';
            end
            if nyquist > 180
                d = data';           
                data = bw_filter(d, header.sampleRate, [175 185], 4, 1, 1)';
            end
            data = detrend(data);
        end

        if notchFilter2 && isAnalog
            nyquist = header.sampleRate/2.0;
            d = data';
            data = bw_filter(d, header.sampleRate, [48 52], 4, 1, 1)';
            if nyquist > 100 
                d = data';
                data = bw_filter(d, header.sampleRate, [95 105], 4, 1, 1)';
            end
            if nyquist > 150
                d = data';           
                data = bw_filter(d, header.sampleRate, [145 155], 4, 1, 1)';
            end
            data = detrend(data);
        end

        if removeOffset
            offset = mean(data);
            data = data - offset;
        end

        if differentiate && isAnalog
            data = diff(data);
            data = [data; 0.0];  % keep num Samples the same!
        end
                
        if rectify
            data = abs(data);
        end
        
        if envelope && isAnalog
            data = abs(hilbert(data));
        end
        
        if invertData
            data = data * -1.0;
        end
                
        dataarray(k,1:header.numSamples) = data;
    end
            
    if showProg
        delete(wbh);    
    end

end

    function drawTrial
        

        %  ** here have to loop over channels and set colours and scales separately *** 

        % get segment of already processed full trial data and scaling
        % factors

        numChannelsToDisplay = numel(selectedChannelList);

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
                        
                        [timebase, fd] = getTrial(k, epochStart);                      
                        maxRange(j) = max(abs(fd));
                        minRange(j) = -maxRange(j);        
                     end
                        
                end
            end
            
            
        end

        % normalize plot data for plotting channels together.
        % normalize to 0.5 of full plot range for readability 
             
        plotMax = 2.0;              
        plotMin = -2.0;

        if ~overlayPlots
            plotMax = plotMax * numChannelsToDisplay;
            plotMin = plotMin * numChannelsToDisplay;
        end
        % add an offset to stack plots
        singleRange = (plotMax-plotMin) / numChannelsToDisplay;
        
        % plot channels
        amplitudeLabels = [];               


        for k=1:numChannelsToDisplay
            chanIdx = selectedChannelList(k);
            chanType = channelTypes(chanIdx);
               
            [timebase, fd] = getTrial(k, epochStart); 
            
            switch chanType
                case num2cell(MEG_CHANNELS)
                    plotColour = 'blue';
                    amplitudeUnits(k)  = {'Tesla'};
                    maxAmp = maxRange(1);
                case num2cell(ADC_CHANNELS)
                    plotColour = darkGreen;
                    amplitudeUnits(k) = {'Volts'};
                    maxAmp = maxRange(2);
                case num2cell(TRIGGER_CHANNELS) 
                    plotColour = orange;
                    amplitudeUnits(k) = {' '};
                    maxAmp = maxRange(3);
                case num2cell(DIGITAL_CHANNELS)
                    plotColour = 'black';
                    amplitudeUnits(k) = {'Bits'};
                    maxAmp = maxRange(4);
                case num2cell(EEG_CHANNELS)
                    plotColour = 'cyan';
                    amplitudeUnits(k) = {'Volts'};
                    maxAmp = maxRange(5);
                otherwise
                    plotColour = 'magenta';
                    amplitudeUnits(k) = {' '};
                    maxAmp = maxRange(6);
            end

            if badChannelMask(chanIdx) == 1 || badTrialMask(trialNo) == 1
                plotColour = [0.8 0.8 0.8];
            end
            
            if selectedMask(k) == 1
                plotColour = 'red';
            end


            pd = fd ./ maxAmp;

            chanIndex = selectedChannelList(k);
            channelName = char( channelNames(chanIndex) );
                         
            % add offset...
            offset = 0.0;
            if ~overlayPlots
                start = plotMin + (singleRange / 2.0);        
                offset =  start + ((numChannelsToDisplay-k) * singleRange);
            end
            pd = pd + offset;
            plot(timebase,pd,'Color',plotColour);
%              plot(timebase,pd,'Color',plotColour,'UserData',k,'ButtonDownFcn',@clickedOnLabel,'ContextMenu',labelContextMenu);

            if k==1
                ylim([plotMin plotMax])
                xlim([timebase(1) timebase(end)])
            end
            % plot baseline 
            if numChannelsToDisplay > 1 
                v = [offset offset];
                h = xlim;            
                line(h,v, 'color', 'black');
            end        

            % plot channel Name and amplitude
            if  ~overlayPlots
                x = epochStart - (epochTime * 0.035);
                y = offset;
                if ~enableMarking
                    s = sprintf('%s', channelName);
                    text(x,y,s,'color',plotColour,'interpreter','none','UserData',k,'ButtonDownFcn',@clickedOnLabel);
                end
           
                sample = round( (cursorLatency - header.epochMinTime - epochStart) * header.sampleRate) + 1;
                if sample > length(fd), sample = length(fd); end
                if sample < 1, sample = 1; end
                amplitude = fd(sample);
                s = sprintf('%.2g %s', amplitude, char(amplitudeUnits(k)));
                x = epochStart + epochTime + (epochTime * 0.007);
                amplitudeLabels(k) = text(x,y,s,'color',plotColour,'interpreter','none');

            end

            hold on;        

        end            
             
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
        

        nsamples = length(timebase);
        if enableMarking
  
            th = ones(nsamples,1) * threshold';

            plot(timebase, th', 'r:', 'lineWidth',1.5);

            th = ones(nsamples,1) * minAmplitude';
            plot(timebase, th', 'g:', 'lineWidth',1.5);

            th = ones(nsamples,1) * maxAmplitude';
            plot(timebase, th', 'c:', 'lineWidth',1.5);
        end

        % check if events exist in this window and draw...      
        if ~isempty(eventList) && numChannelsToDisplay == 1
            events = find(eventList > epochStart & eventList < (epochStart+epochTime));
            
            if ~isempty(events)
                for i=1:length(events)
                    thisEvent = events(i);
                    t = eventList(thisEvent);
                    h = [t,t];
                    v = ylim;
                    cursor = line(h,v, 'color', orange,'linewidth',1);
                    if thisEvent == currentEvent
                        pos = h;
                        latency = eventList(thisEvent);
                        sample = round( (latency - header.epochMinTime) * header.sampleRate) + 1;
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
                    if showMarkerWindow
                        t1 = t + markerWindowStart;
                        t2 = t + markerWindowEnd;
                        xpoints = [t1, t1, t2, t2];
                        ypoints = [v(1), v(2), v(2), v(1)];
                        a = fill(xpoints,ypoints,'blue','linestyle','none');
                        a.FaceAlpha = 0.08;
                    end
                      
                    line(h,v, 'color', 'blue');
                    x = t + epochTime * 0.001;
                    y =  v(2) - (v(2) - v(1))*0.05;
                    s = sprintf('%s', markerName);
                    text(x,y,s,'color','blue','interpreter','none');
                end
            end                         
        elseif currentMarkerIndex == numMarkers+2
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
        
        if k == numChannelsToDisplay
            ax=axis;
            cursorHandle=line([cursorLatency cursorLatency], [ax(3) ax(4)],'color','black');
        end
        
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

% contextmenu on right mouse click on object code only works for later versions of Matlab ...
%     function cm = labelContextMenu
%         cm = uicontextmenu;
%         uimenu(cm,"Text","Set Good","MenuSelectedFcn",@setGood);
%         uimenu(cm,"Text","Set Bad","MenuSelectedFcn",@setBad);       
%     end
% 
%         function setGood(~,~)
%             h = gco;
%             idx = find(selectedMask == 1);
%             chanIdx = selectedChannelList(idx);  
%             badChannelMask(chanIdx) = 0;        % set this channel good
%             selectedMask(:) = 0;
%             drawTrial;
%         end
% 
%         function setBad(~,~)
%             h = gco;
%             idx = find(selectedMask == 1);
%             chanIdx = selectedChannelList(idx); % set this channel bad
%             badChannelMask(chanIdx) = 1;
%             selectedMask(:) = 0;
%             drawTrial;
%         end

    % version 4.0 - new cursor function
    
    function updateCursors           
        if ~isempty(cursorHandle)
            set(cursorHandle, 'XData', [cursorLatency cursorLatency]);      
        end 
        
        s = sprintf('Latency: %.4f s', cursorLatency);
        set(cursor_text, 'string', s);

        if ~isempty(amplitudeLabels)                        
                for k=1:length(selectedChannelList)
                    if ~overlayPlots && ~enableMarking
                       [~,fd] = getTrial(k,epochStart);
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
    end
        
    function buttondown(~,~) 
        
        if isempty(cursorHandle) || isempty(header)
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

 
        fd = dataarray(1,:);
        mx = max(fd);
        fd = fd ./ mx;

        while (true)          
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

        set(eventCtl(:),'enable','on');
        s = sprintf('# events = %d', numEvents);
        set(numEventsTxt,'String',s);
         
    end



   % need dialog to create conditional events ...
   
    function create_event_callback(~,~)

        latencies = bw_conditionalMarker(markerFileName);
        eventList = latencies;      
        numEvents = length(eventList);
      
        s = sprintf('# events = %d', numEvents);
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
  
        new_list =importdata(loadlatfull);
        
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
        
        currentEvent = 1;
        drawTrial;
        
        set(eventCtl(:),'enable','on');
        s = sprintf('# events = %d', numEvents);
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
        s = sprintf('# events = %d', numEvents);
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
        
        s = sprintf('# events = %d', numEvents);
        set(numEventsTxt,'String',s);
        
        % first_event_callback;
        
    end

    % % save latencies
    % function save_events_callback(~,~)
    %     if isempty(eventList)
    %         errordlg('No events defined ...');
    %         return;
    %     end
    % 
    %     saveName = strcat(dsName, filesep, '*.txt');
    %     [name,path,~] = uiputfile('*.txt','Save Event latencies in File:',saveName);
    %     if isequal(name,0)
    %         return;
    %     end         
    % 
    %     eventFile = fullfile(path,name);
    %     fprintf('Saving event times to text file %s \n', eventFile);
    % 
    %     fid = fopen(eventFile, 'w');
    %     fprintf(fid,'%.5f\n', eventList');
    %     fclose(fid);
    % end

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



