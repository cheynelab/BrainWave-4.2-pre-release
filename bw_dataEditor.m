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
data = [];
dataarray = [];
markerFileName = [];

numMarkers = 0;
currentMarkerIndex = 1;
markerNames = {};
markerLatencies = {};
currentMarker = 1;

channelName = [];
channelMenuItems = [];
maxChannels = 16;       % max channels to plot at once 

eventList = [];  
currentEvent = 1;
numEvents = 0;

threshold = 0;
maxAmplitude = 0.0;
minAmplitude = 0.0;
minSeparation = 0.0;

header = [];
channelNames = [];
channelTypes = [];
selectedChannelList = [];
selectedChannelNames = [];
channelMenuIndex = 1;
customChannelList = 1;

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

bandPass = [1 50];
minDuration = 0.01;
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
channelSets = {'Custom';'MEG Sensors';'ADC Channels';'Digital Channels';'Trigger Channel'};
channelSetIndex = [-1, 5, 18, 20, 19];

darkGreen = [0.1 0.7 0.1];
orange = [0.9 0.6 0.1];
gray = [0.3 0.3 0.3];

% set defaults
% Draw arrows - calls:  uparrow.m and downarrow.m - %%%% ADDED BY CECILIA %%%%
uparrow_im=draw_uparrow;
downarrow_im=draw_downarrow;
rightarrow_im=draw_rightarrow;
leftarrow_im=draw_leftarrow;


fh = figure('numbertitle','off','position',[200, 800, 1600, 1200],...
    'Name','Data Editor', 'Color','white','menubar','none','WindowButtonUpFcn',@stopdrag,'WindowButtonDownFcn',@buttondown);

if ispc
    movegui(fh, 'center');
end

filemenu=uimenu('label','File');
uimenu(filemenu,'label','Open Dataset','accelerator','O','callback',@openFile_callback)
uimenu(filemenu,'label','Close','accelerator','W','separator','on','callback',@quit_filemenu_callback)

channelMenu=uimenu('label','ChannelSets');

markerMenu=uimenu('label','Edit Markers');
uimenu(markerMenu,'label','Import Markers from Text File...','callback',@load_events_callback)
uimenu(markerMenu,'label','Import Markers from KIT Event File...','callback',@load_KIT_events_callback)
uimenu(markerMenu,'label','Export Markers to Excel...','separator','on','callback',@save_marker_as_excel_callback)
uimenu(markerMenu,'label','Export Markers to Text File...','callback',@save_events_callback)

% +++++++++++++ set plot window +++++++++++++
ph = subplot('position',[0.05 0.3 0.9 0.68]);
    

function openFile_callback(~,~)
    dsName = uigetdir('.ds', 'Select CTF dataset ...');
    if dsName == 0
        return;
    end
    
    initData;
    drawTrial;

end

function quit_filemenu_callback(~,~)
    close(fh);
end

function initData
        
    if ~exist(dsName,'file')
        fprintf('Could not find file %s...\n', dsName);
        return;      
    end

    eventList = [];  
    currentEvent = 1;
    numEvents = 0;
    
    threshold = 0;
    maxAmplitude = 0.0;
    minAmplitude = 0.0;
    minSeparation = 0.0;
    
    header = bw_CTFGetHeader(dsName);
    
    if header.numTrials > 1
        errordlg('eventMarker can only be used with single-trial (raw) data...');
        return;
    end
    
    [~,n,e] = fileparts(dsName);
    tStr = sprintf('Data Editor: %s', [n e]);
    set(fh,'Name', tStr);

    t = round( (cursorLatency - header.epochMinTime) * header.sampleRate) + 1;
    s = sprintf('Cursor = %.4f s (%d)', cursorLatency, t);
    set(cursor_text,'string',s);

    sliderScale = [(epochTime / header.trialDuration) * 0.02 (epochTime / header.trialDuration) * 0.08];
    set(latency_slider, 'sliderStep',sliderScale);
    
    timeVec = [];
    data = [];
    dataarray = [];
    
    numMarkers = 0;
    currentMarkerIndex = 1;
    currentMarker = 1;
    set(marker_Popup,'value',1)
    markerLatencies = {};
    set(markerIncButton,'enable','off')
    set(markerDecButton,'enable','off')

    
    minScale = NaN; % forces autoscale;
    maxScale = NaN;
    
    epochStart = header.epochMinTime;
    epochTime = 10;
    if epochTime > header.epochMaxTime
        epochTime = header.epochMaxTime;
    end
    s = sprintf('%.4f', epochTime);
    set(epochDurationEdit,'string', s)

    cursorLatency = epochTime / 2.0;

    epochSamples = round(epochTime * header.sampleRate);

    bandpass(1) = 1;
    bandpass(2) = header.sampleRate / 2.0;
    set(filt_low_pass,'string',bandpass(2))
    
    fprintf('Loading data...\n\n');
    
    markerFileName = strcat(dsName,filesep,'MarkerFile.mrk');
    
    markerNames = {'none'};
    if exist(markerFileName,'file')
        fprintf('found marker file %s\n', markerFileName); 
        [names, markerData] = bw_readCTFMarkerFile(markerFileName);
       
        % drop trial numbers for now
        numMarkers = size(names,1);
        fprintf('dataset has %d markers ...\n', numMarkers); 
        if (numMarkers > 0)
            for j = 1:numMarkers
                x = markerData{j}; 
                markerLatencies{j} = x(:,2);
                markerNames{j+1} = names{j};
            end
            markerNames{numMarkers+2} = 'All Markers';
        end
        set(marker_Popup,'string',markerNames)
    else
       fprintf('no marker file found...\n'); 
       numMarkers = 0;
    end
    

    % ** need to rebuild channel set menu with valid channel types *** 
    

    s = sprintf('Sample Rate: %4.1f S/s',header.sampleRate);
    set(sampleRateTxt,'string',s);

    s = sprintf('MEG Sensors: %d',header.numSensors);
    set(numSensorsTxt,'string',s);

    s = sprintf('Total Channels: %d',header.numChannels);
    set(numChannelsTxt,'string',s);

    longnames = {header.channel.name};   
    channelNames = cleanChannelNames(longnames);
    numChannels = length(channelNames);

    channelTypes = [header.channel.sensorType];
    updateChannelMenu;

    channelMenuIndex = 1;
    customChannelList = 1;
    selectedChannelList = 1;
    selectedChannelNames = channelNames(selectedChannelList);


    loadData;   

end

function updateChannelMenu
        
    if ~isempty(channelMenuItems)
        delete(channelMenuItems(:)) 
    end
    for k=1:numel(channelSets)
        s = char(channelSets(k));
        channelMenuItems(k) = uimenu(channelMenu,'label',s,'callback',@channel_menu_callback); 
        menuIndex = channelSetIndex(k);

        % turn off channel types that don't exist
        if menuIndex > 1
            idx = find(channelTypes == menuIndex);
            if isempty(idx)
                set(channelMenuItems(k),'enable','off');
            end
        end
        
        if channelMenuIndex == k
            set(channelMenuItems(k),'checked','on')
        else
            set(channelMenuItems(k),'checked','off')
        end
    end
    channelMenuItems(end+1) = uimenu(channelMenu,'label','Select Channels',...
        'separator','on','callback',@editChannelSet_callback); 
    
end
    
      
function channel_menu_callback(src,~)

    if isempty(channelNames)
        return;
    end
    
    channelMenuIndex = get(src,'position'); 
    menuItems = get(channelMenu,'Children');

    if channelMenuIndex == 1
        selectedChannelList = customChannelList;
    else

        % set by channel type - need to flag too many channels?
        nchans = numel(channelNames);
        channelExcludeFlags = ones(1,nchans);
             
        chanType = channelSetIndex(channelMenuIndex);
        for i=1:nchans    
            if (channelTypes(i) == chanType)
                channelExcludeFlags(i) = 0;             
            end
        end
        
        selectedChannelList = find(channelExcludeFlags == 0);
        if numel(selectedChannelList) > maxChannels
            s = sprintf('Do you want to plot %d channels at once?',numel(selectedChannelList));        
            r = questdlg(s,'Data Editor','Yes','No','No');
            if strcmp(r,'No')
                channelMenuIndex = 1;
                selectedChannelList = customChannelList;
                set(menuItems(:),'Checked','off');
                set(menuItems(end),'Checked','on'); % first item is last! 
                return;
            end
        end
    end

    selectedChannelNames = channelNames(selectedChannelList);
        
    set(get(channelMenu,'Children'),'Checked','off');
    set(src,'Checked','on');

    if numel(selectedChannelNames) == 1
        setMarkingCtls('on');
        set(markerEventsCheck,'enable','on');
        reset_events_button;
    else
        setMarkingCtls('off');
        set(markerEventsCheck,'enable','off');
    end

    loadData; 
    drawTrial;
    autoScale_callback;

      
end

function editChannelSet_callback(~,~)  
       
    [selected, menuIndex] = bw_channelSelector(header, channelSets, selectedChannelList, channelMenuIndex);
    if isempty(selected)
        return;
    end   
    selectedChannelList = selected;
    channelMenuIndex = menuIndex;
    if menuIndex == 1
        customChannelList = selectedChannelList;
    end
    selectedChannelNames = channelNames(selectedChannelList);

    loadData;
    drawTrial;
    autoScale_callback;

    fprintf('Number of Selected Channels = %d\n', length(selectedChannelList));
end

%+++++++++ event detection controls +++++++++  

annotation('rectangle',[0.6 0.02 0.35 0.18],'EdgeColor','blue');

markerEventsCheck = uicontrol('style','checkbox','units','normalized','position',[0.63 0.185 0.08 0.025],...
    'string','Mark Events','backgroundcolor','white','foregroundcolor','blue','fontweight','bold',...
    'FontSize',11, 'value',enableMarking,'callback', @enable_marking_callback);

function enable_marking_callback(src,~)
    enableMarking = get(src,'value');
    if enableMarking
        setMarkingCtls('on');
    else
        setMarkingCtls('off');
    end
    drawTrial;
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
    set(reset_events_button,'enable',state)  
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

reset_events_button = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.74 0.03 0.05 0.025],...
    'enable','off','string','Update','Foregroundcolor','blue','backgroundcolor','white','callback',@reset_events_callback);


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
            drawTrial;       
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
            drawTrial;       
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

    function reset_events_callback(~,~)
        % reset params for marking 
        peakAmp = maxScale * 0.5;
        maxAmplitude = peakAmp * 1.2;
        minAmplitude = 0.0;
        threshold = peakAmp * 0.2;

        s = sprintf('%.2g', maxAmplitude);
        set(max_amplitude_edit,'string',s);
        s = sprintf('%.2g', minAmplitude);
        set(min_amplitude_edit,'string',s);
        s = sprintf('%.2g', threshold);
        set(threshold_edit,'string',s);
        drawTrial;
        
    end



% +++++++++ amplitude and time controls +++++++++

cursor_text = uicontrol('style','text','fontsize',12,'units','normalized','position',[0.78 0.305 0.15 0.02],...
     'string','Cursor =','BackgroundColor','white','foregroundColor',[0.7,0.41,0.1]);

latency_slider = uicontrol('style','slider','units', 'normalized',...
    'position',[0.05 0.22 0.9 0.05],'min',0,'max',1,'Value',0,...
    'sliderStep', [1 100],'BackGroundColor',[0.8 0.8 0.8],'ForeGroundColor',...
    'white'); 

% this callback is called everytime slider value is changed. 
% replaces slider callback function

addlistener(latency_slider,'Value','PostSet',@slider_moved_callback);

    function slider_moved_callback(~,~)   
       val = get(latency_slider,'Value');
       epochStart = (val * header.trialDuration) - header.epochMinTime;
       if (epochStart + epochTime > header.epochMaxTime)
            epochStart = header.epochMaxTime - epochTime;
       end
       
       if epochStart < header.epochMinTime
           epochStart = header.epochMinTime;
       end
       drawTrial;
    end


% add popup for channel types for scaling...






uicontrol('style','text','units','normalized','position',[0.14 0.22 0.1 0.02],...
    'string','Max:','fontsize',11,'backgroundcolor','white','horizontalalignment','left');
max_scale=uicontrol('style','edit','units','normalized','position',[0.16 0.22 0.04 0.025],...
    'FontSize', 11, 'BackGroundColor','white','string',maxScale,...
    'callback',@max_scale_callback);

    function max_scale_callback(src,~)
        val =str2double(get(src,'string'));
        if val < minScale
            s = sprintf('%.3g', maxcale);
            set(max_scale,'string',s);
            return;
        end
        maxScale = val;
        drawTrial;
    end

uicontrol('style','text','units','normalized','position',[0.22 0.22 0.08 0.02],...
    'string','Min:','fontsize',11,'backgroundcolor','white','horizontalalignment','left');
min_scale=uicontrol('style','edit','units','normalized','position',[0.24 0.22 0.04 0.025],...
    'FontSize', 11, 'BackGroundColor','white','string',minScale,...
    'callback',@min_scale_callback);

    function min_scale_callback(src,~)
        val =str2double(get(src,'string'));
        if val > maxScale
            s = sprintf('%.3g', minScale);
            set(min_scale,'string',s);
            return;
        end
        minScale = val;
        drawTrial;
    end

    function scaleUp_callback(~,~)
        inc = maxScale * 0.1;
        maxScale = maxScale  - inc;
        minScale = minScale + inc;
        s = sprintf('%.3g', maxScale);
        set(max_scale,'String', s);
        s = sprintf('%.3g', minScale);
        set(min_scale,'String', s);
        drawTrial;
    end

    function scaleDown_callback(~,~)
        inc = maxScale * 0.1;
        maxScale = maxScale  + inc;
        minScale = minScale - inc;
        s = sprintf('%.3g', maxScale);
        set(max_scale,'String', s);
        s = sprintf('%.3g', minScale);
        set(min_scale,'String', s);
        drawTrial;
    end


     
uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.3 0.22 0.03 0.025],...
    'CData',uparrow_im,'Foregroundcolor','black','backgroundcolor','white','callback',@scaleUp_callback);

uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.34 0.22 0.03 0.025],...
    'CData',downarrow_im,'Foregroundcolor','black','backgroundcolor','white','callback',@scaleDown_callback);


uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.38 0.22 0.06 0.025],...
    'Foregroundcolor','black','string','Autoscale','backgroundcolor','white','callback',@autoScale_callback);
    function autoScale_callback(~,~)
        maxScale = NaN;
        drawTrial;
    end



markerNames = {'none'};
uicontrol('style','text','fontsize',14,'units','normalized','horizontalalignment','left','position',...
     [0.48 0.21 0.08 0.03],'string','Marker:','BackgroundColor','white','foregroundcolor','black');
marker_Popup =uicontrol('style','popup','units','normalized','fontsize',12,'position',[0.52 0.19 0.1 0.05],...
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

markerDecButton = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.62 0.22 0.03 0.025],...
    'CData',leftarrow_im,'Foregroundcolor','black','backgroundcolor','white','callback',@marker_dec_callback);
markerIncButton = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.656 0.22 0.03 0.025],...
    'CData',rightarrow_im,'Foregroundcolor','black','backgroundcolor','white','callback',@marker_inc_callback);

    function marker_inc_callback(~,~)  
        if currentMarkerIndex == 1
            return;
        end
        markerTimes = markerLatencies{currentMarkerIndex-1};
        if currentMarker < numel(markerTimes)
            currentMarker = currentMarker + 1;
            epochStart = markerTimes(currentMarker) - (epochTime / 2);
            drawTrial;       
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
            drawTrial;       
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

        if isnan(t) || epochStart+t > header.epochMaxTime
            t = header.epochMaxTime;
            set(src,'string',t);
        end    
        if t < (1.0 / header.sampleRate) * 2.0
            t = (1.0 / header.sampleRate) * 2.0;
            set(src,'string',t);
        end
        epochTime = t;
       
        epochSamples = round(epochTime * header.sampleRate);
        sliderScale = [(epochTime / header.trialDuration) * 0.02 (epochTime / header.trialDuration) * 0.08];
        set(latency_slider,'sliderStep',sliderScale);
        drawTrial;
    end


% +++++++++++++++++++++++++

% ++++++++++ plot settings 

annotation('rectangle',[0.05 0.02 0.53 0.18],'EdgeColor','blue');
uicontrol('style','text','fontsize',11,'units','normalized','position',...
     [0.08 0.18 0.1 0.025],'string','Plot Settings','BackgroundColor','white','foregroundcolor','blue','fontweight','b');

sampleRateTxt = uicontrol('style','text','units','normalized','position',[0.1 0.135 0.08 0.02],...
    'string','Sample Rate:','backgroundcolor','white','FontSize',11,'callback',@filter_check_callback);

numSensorsTxt = uicontrol('style','text','units','normalized','position',[0.18 0.135 0.08 0.02],...
    'string','MEG Sensors:','backgroundcolor','white','FontSize',11,'callback',@filter_check_callback);

numChannelsTxt = uicontrol('style','text','units','normalized','position',[0.26 0.135 0.08 0.02],...
    'string','Total Channels:','backgroundcolor','white','FontSize',11,'callback',@filter_check_callback);

uicontrol('style','checkbox','units','normalized','position',[0.1 0.1 0.04 0.02],...
    'string','Filter','backgroundcolor','white','value',~filterOff,'FontSize',11,'callback',@filter_check_callback);

% replace with 50 / 60 Hz checks - check code 
uicontrol('style','checkbox','units','normalized','position',[0.1 0.07 0.06 0.02],...
    'string','60 Hz Notch','backgroundcolor','white','value',notchFilter,'FontSize',11,'callback',@notch_check_callback);
uicontrol('style','checkbox','units','normalized','position',[0.18 0.07 0.06 0.02],...
    'string','50 Hz Notch','backgroundcolor','white','value',notchFilter2,'FontSize',11,'callback',@notch2_check_callback);

uicontrol('style','checkbox','units','normalized','position',[0.45 0.1 0.1 0.02],...
    'string','Remove Offset','backgroundcolor','white','value',overlayPlots,'FontSize',11,'callback',@remove_offset_callback);

uicontrol('style','checkbox','units','normalized','position',[0.1 0.04 0.08 0.02],...
    'string','Invert','backgroundcolor','white','value',invertData,'FontSize',11,'callback',@invert_check_callback);

uicontrol('style','checkbox','units','normalized','position',[0.18 0.04 0.08 0.02],...
    'string','Rectify','backgroundcolor','white','value',rectify,'FontSize',11,'callback',@rectify_check_callback);

uicontrol('style','checkbox','units','normalized','position',[0.26 0.04 0.12 0.02],...
    'string','Differentiate','backgroundcolor','white','value',differentiate,'FontSize',11,'callback',@firstDiff_check_callback);

uicontrol('style','checkbox','units','normalized','position',[0.34 0.04 0.08 0.02],...
    'string','Envelope (hilbert)','backgroundcolor','white','value',envelope,'FontSize',11,'callback',@envelope_check_callback);


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
        set(filt_hi_pass, 'enable','off');
        set(filt_low_pass, 'enable','off');
    else
        set(filt_hi_pass, 'enable','on');
        set(filt_low_pass, 'enable','on');
    end
    
    loadData;
    drawTrial;
end

uicontrol('style','checkbox','units','normalized','position',[0.45 0.14 0.1 0.02],...
    'string','Overlay Plots','backgroundcolor','white','value',overlayPlots,'FontSize',11,'callback',@overlay_plots_callback);

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
        numChannelsToDisplay = numel(selectedChannelNames);
        dataarray = zeros(numChannelsToDisplay, header.numSamples);


        for k=1:numChannelsToDisplay
            channelName = selectedChannelNames{k};
            
            if filterOff
                [timeVec, data] = bw_CTFGetChannelData(dsName, channelName);
            else
                [timeVec, data] = bw_CTFGetChannelData(dsName, channelName, bandPass);
            end
            
            nyquist = header.sampleRate/2.0;
                    
            if (notchFilter)
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
    
            if (notchFilter)
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
    
            if differentiate
                data = diff(data);
                data = [data; 0.0];  % keep num Samples the same!
            end
                    
            if rectify
                data = abs(data);
            end
            
            if envelope
                data = abs(hilbert(data));
            end
            
            if invertData
                data = data * -1.0;
            end
                    
            dataarray(k,1:header.numSamples) = data;
        end

    end

    function drawTrial
        

        %  ** here have to loop over channels and set colours and scales separately *** 

        % get segment of already processed full trial data and scaling
        % factors

        numChannelsToDisplay = numel(selectedChannelNames);
        maxAmps = [];
        plotData = [];

        % get global scale and timebase for this 
        for k=1:numChannelsToDisplay
            [timebase, fd] = getTrial(k, epochStart); 
            maxAmps(k) = max(fd);
            plotData(k,:) = fd;
        end

        % autoscale if scale empty
        if isnan(minScale) || isnan(maxScale)
            maxScale = max(maxAmps) * 2;
            if maxScale == 0
                maxScale = 1;
            end
            minScale = -maxScale;
            s = sprintf('%.3g', maxScale);
            set(max_scale,'String', s);
            s = sprintf('%.3g', minScale);
            set(min_scale,'String', s);
        end
            
        plotMax = maxScale; 
        plotMin = minScale;
        if ~overlayPlots
            plotMax = plotMax * numChannelsToDisplay;
            plotMin = plotMin * numChannelsToDisplay;
        end

        for k=1:numChannelsToDisplay
            chanIndex = selectedChannelList(k);
            chanType = channelTypes(chanIndex);
            channelName = char(selectedChannelNames{k});

            % ** check channel type 
            if chanType == 5
                plotColour = 'blue';
            elseif chanType == 18
                plotColour = darkGreen;
            elseif chanType == 20
                plotColour = 'black';
            elseif chanType == 19
                plotColour = gray;
            else
                plotColour = 'black';
            end
    
            % add offset...
            offset = 0.0;

            % add an offset to stack plots
            singleRange = (plotMax-plotMin) / numChannelsToDisplay;

            if ~overlayPlots
                start = plotMin + (singleRange / 2.0);        
                offset =  start + ((numChannelsToDisplay-k) * singleRange);
            end
            fd = plotData(k,:) + offset;

            plot(timebase,fd,'Color',plotColour);
     
            % draw baseline and channel name
            s = sprintf('%s', channelName);
            x = epochStart + epochTime * 0.008;

            % plot baseline
            if numChannelsToDisplay > 1
                v = [offset offset];
                h = xlim;            
                line(h,v, 'color', 'black');
            end

            % plot channel Name
            if numChannelsToDisplay == 1 || ~overlayPlots
                y = offset + singleRange * 0.2;
                text(x,y,s,'color',plotColour,'interpreter','none');
            end
            hold on;        

        end            

                
        ylim([plotMin plotMax]);
        set(gca,'ytick',[])



        % plot xscale and other markers
   
        xlim([timebase(1) timebase(end)])
        xlabel('Time (sec)', 'fontsize', 12);
        nsamples = length(timebase);
        if enableMarking && numChannelsToDisplay == 1
  
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
            markerTimes = markerLatencies{currentMarkerIndex-1};
            markerName = char( markerNames{currentMarkerIndex} );
            markers = find(markerTimes > epochStart & markerTimes < (epochStart+epochTime));
            
            if ~isempty(markers)
                for k=1:length(markers)
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
        if enableMarking && numChannelsToDisplay == 1
            tt = legend(s, 'threshold', 'min. amplitude', 'max. amplitude');
        else
            tt = legend(s);
        end
        set(tt,'interpreter','none','Autoupdate','off');
        
      
        ax=axis;
        cursorHandle=line([cursorLatency cursorLatency], [ax(3) ax(4)],'color',[0.8,0.4,0.1]);

        
        hold off;
        
        
    end

    % get processed trial data...
    function [timebase, fd] = getTrial(channel, startTime)
              
        % check trial boundaries
       
        if startTime < header.epochMinTime 
            startTime = header.epochMinTime;
        end 
        
        if startTime + epochTime > header.epochMaxTime
            startTime = header.epochMaxTime-epochTime;
        end
        
        epochStart = startTime;
        
        % get data - note data indices start at 1;
        
        startSample = round( (epochStart - header.epochMinTime) * header.sampleRate) + 1;
        endSample = startSample + epochSamples;
       
        fd = dataarray(channel, startSample:endSample);
        timebase = timeVec(startSample:endSample);
               
    end

    % version 4.0 - new cursor function
    
    function updateCursors                 
        if ~isempty(cursorHandle)
            set(cursorHandle, 'XData', [cursorLatency cursorLatency]);      
        end 
        sample = round( (cursorLatency - header.epochMinTime) * header.sampleRate) + 1;
       
        s = sprintf('Cursor = %.4f s (%d)', cursorLatency, sample);
        set(cursor_text, 'string', s);

    end
        
    function buttondown(~,~) 
        
        if isempty(cursorHandle)
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
        
        while (true)          
            value = data(sample);
            
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
                    value = data(sample);
                end  
                
                if sample > header.numSamples || sample < 1
                    break;
                end
                
                if reverseScan
                    eventData = data(sample:eventSample);
                else
                    eventData = data(eventSample:sample);
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
        success = write_MarkerFile(dsName, trig);  
        
        % if not cancelled save file 
        if success
            numMarkers = numMarkers + 1;
            % add to list
            markerNames{numMarkers+1} = newName;
            markerLatencies{numMarkers} = eventList;
            markerNames{numMarkers+2} = 'All Markers';
            set(marker_Popup,'string',markerNames);
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

    if ~exist('dsName','var')
        dsName = uigetdir('.ds', 'Select CTF dataset ...');
        if dsName == 0
            return;
        end    
    end

    % draw after controls defined.
    initData;
    drawTrial;

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function bw_write_MarkerFile(dsName,trig)
%   Writes a markerFile.mrk for the given input 'trig'
%   Returns on error if a markerFile.mrk exists
%   INPUT:
%         dsName = ctf dataset name
%         trig = 1 X N(num of trigs) structure
%                trig(1:N).ch_name 
%                         .onset_idx %sample index of trigger onsets                 
%                         .times     %dataset times of triggers onsets
%
%  pferrari@meadowlandshospital.org, Aug2012
%
% this is a modified version for eventMarker - apply changes to bw_write_MarkerFile?
%
function result = write_MarkerFile(dsName,trig)

result = 0;

no_trigs=numel(trig);

filepath=strcat(dsName, filesep, 'MarkerFile.mrk');
          
fprintf('writing marker file %s\n', filepath);

fid = fopen(filepath,'w','n');
fprintf(fid,'PATH OF DATASET:\n');
fprintf(fid,'%s\n\n\n',dsName);
fprintf(fid,'NUMBER OF MARKERS:\n');
fprintf(fid,'%g\n\n\n',no_trigs);

for i = 1:no_trigs
    
    fprintf(fid,'CLASSGROUPID:\n');
    fprintf(fid,'3\n');
    fprintf(fid,'NAME:\n');
    fprintf(fid,'%s\n',trig(i).ch_name);
    fprintf(fid,'COMMENT:\n\n');
    fprintf(fid,'COLOR:\n');
    fprintf(fid,'blue\n');
    fprintf(fid,'EDITABLE:\n');
    fprintf(fid,'Yes\n');
    fprintf(fid,'CLASSID:\n');
    fprintf(fid,'%g\n',i);
    fprintf(fid,'NUMBER OF SAMPLES:\n');
    fprintf(fid,'%g\n',length(trig(i).times));
    fprintf(fid,'LIST OF SAMPLES:\n');
    fprintf(fid,'TRIAL NUMBER\t\tTIME FROM SYNC POINT (in seconds)\n');
    for t = 1:length(trig(i).times)-1
        fprintf(fid,'                  %+g\t\t\t\t               %+0.6f\n',0,trig(i).times(t));
    end
    fprintf(fid,'                  %+g\t\t\t\t               %+0.6f\n\n\n',0,trig(i).times(end));
end

fclose(fid);

result = 1;

end



% ++++++++++++ old event controls 
% 
% evt_ctrl(1) = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.6 0.2 0.06 0.03],...
%     'string','First Event','Foregroundcolor','blue','backgroundcolor','white','callback',@first_event_callback);
% 
% evt_ctrl(2) = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.68 0.2 0.06 0.03],...
%     'string','Last Event','Foregroundcolor','blue','backgroundcolor','white','callback',@last_event_callback);
% 
% evt_ctrl(3) = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.78 0.2 0.09 0.03],...
%     'string','Previous Event','Foregroundcolor','blue','backgroundcolor','white','callback',@event_dec_callback);
% 
% evt_ctrl(4) = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.88 0.2 0.07 0.03],...
%     'string','Next Event','Foregroundcolor','blue','backgroundcolor','white','callback',@event_inc_callback);
% 
% evt_ctrl(5) = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.05 0.2 0.07 0.03],...
%     'string','Add Event','Foregroundcolor', [0.3 0.6 0],'backgroundcolor','white','callback',@add_callback);
% 
% evt_ctrl(6) = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.13 0.2 0.07 0.03],...
%     'string','Delete Event','Foregroundcolor','red','backgroundcolor','white','callback',@delete_callback);
% 
% evt_ctrl(7) = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.22 0.2 0.07 0.03],...
%     'string','Clear Events','Foregroundcolor','black','backgroundcolor','white','callback',@delete_all_callback);
% 
% evt_ctrl(8) = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.33 0.2 0.1 0.03],...
%     'string','Load Marker Events...','Foregroundcolor','blue','backgroundcolor','white','callback',@load_marker_events_callback);
% 
% evt_ctrl(9) = uicontrol('style','pushbutton','units','normalized','fontsize',11,'position',[0.45 0.2 0.1 0.03],...
%     'string','Include/Exclude Events','Foregroundcolor','blue','backgroundcolor','white','callback',@filter_events_callback);
% 
% set(evt_ctrl(:),'enable','off')

   %  function first_event_callback(~,~)   
   %      if numEvents < 1
   %          return;
   %      end
   %      currentEvent = 1;
   %      epochStart = eventList(currentEvent) - (epochTime / 2);
   %      drawTrial;   
   %      % adjust slider position
   %      val = (epochStart + header.epochMinTime) / header.trialDuration;
   %      if val < 0, val = 0.0; end
   %      if val > 1.0, val = 1.0; end
   %      set(latency_slider, 'value', val);
   %  end
   % 
   %  function last_event_callback(~,~)   
   %      if numEvents < 1
   %          return;
   %      end
   %      currentEvent = numEvents;
   %      epochStart = eventList(currentEvent) - (epochTime / 2);
   %      drawTrial;   
   %      % adjust slider position
   %      val = (epochStart + header.epochMinTime) / header.trialDuration;
   %      if val < 0, val = 0.0; end
   %      if val > 1.0, val = 1.0; end
   %      set(latency_slider, 'value', val);
   %  end

   % 


   %  function filter_events_callback(~,~)
   %      if numEvents < 1
   %         errordlg('No events defined ...');
   %         return;
   %      end
   % 
   %      markerTimes = markerLatencies{currentMarkerIndex-1};
   % 
   %      wStart = markerWindowStart;
   %      wEnd = markerWindowEnd;
   % 
   %      windowEvents = [];
   %      for k=1:numEvents  
   %          latency = eventList(k);
   %          for j=1:numel(markerTimes)
   %              wStart = markerTimes(j) + markerWindowStart;               
   %              wEnd = markerTimes(j) + markerWindowEnd;
   %              if latency > wStart && latency < wEnd
   %                  windowEvents(end+1) = k;
   %              end
   %          end
   %      end     
   % 
   %     if isempty(windowEvents)
   %         errordlg('No Events found within Marker Window ...');
   %         return;
   %     end
   % 
   %     response = questdlg('If event is within Marker window','Event Marker','Include Event','Exclude Event','Cancel','Cancel');
   % 
   %     if strcmp(response,'Include Event')
   %          s = sprintf('Including %d of %d events...Continue? (Cannot be undone)', numel(windowEvents), numEvents);
   %          response = questdlg(s,'Event Marker','Yes','No','No');
   %          if strcmp(response,'Yes')
   %              eventList = eventList(windowEvents);
   %          end  
   %     elseif strcmp(response,'Exclude Event')
   %          s = sprintf('Excluding %d of %d events...Continue? (Cannot be undone)', numel(windowEvents), numEvents);
   %          response = questdlg(s,'Event Marker','Yes','No','No');
   %          if strcmp(response,'Yes')
   %              eventList(windowEvents) = [];
   %          end  
   %     end       
   %     numEvents = numel(eventList);
   % 
   %     s = sprintf('No. events = %d', numEvents);
   %      set(numEventsTxt,'String',s);
   % 
   %      % since event list has changed goto first event
   %      first_event_callback;
   % 
   %      drawTrial;
   %  end
   % 
   % 
    % markerW(1) = uicontrol('style','checkbox','units','normalized','fontsize',12,'position',[0.47 0.955 0.13 0.04],...
    %     'string','Exclude/Include Window','value',0,'Foregroundcolor','black','backgroundcolor','white','callback',@show_window_callback);
    % 
    %     function show_window_callback(src,~)    
    %         showMarkerWindow = get(src,'value');
    %         drawTrial;
    %     end
    % 
    % markerW(2) = uicontrol('style','edit','units','normalized','position',[0.61 0.96 0.04 0.035],...
    %     'FontSize', 11, 'BackGroundColor','white','string',markerWindowStart,'callback',@windowStart_edit_callback);
    %     function windowStart_edit_callback(src,~)
    %         markerWindowStart =str2double(get(src,'string'));
    %         drawTrial;
    %     end
    % markerW(3) = uicontrol('style','text','fontsize',11,'units','normalized','position',[0.66 0.96 0.02 0.03],...
    %      'string','to','BackgroundColor','white');
    % markerW(4)  = uicontrol('style','edit','units','normalized','position',[0.68 0.96 0.04 0.035],...
    %     'FontSize', 11, 'BackGroundColor','white','string',markerWindowEnd,'callback',@windowEnd_edit_callback);
    %     function windowEnd_edit_callback(src,~)
    %         markerWindowEnd =str2double(get(src,'string'));
    %         drawTrial;
    %     end

    % 
    % markerW(5) = uicontrol('style','text','fontsize',11,'units','normalized','position',[0.72 0.96 0.05 0.03],...
    %      'string','seconds','BackgroundColor','white');
    
    % set(markerW(:),'enable','off');




