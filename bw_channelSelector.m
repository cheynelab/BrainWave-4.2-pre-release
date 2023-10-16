function [selectedChannels] = bw_channelSelector(header,oldSelected, badChannelMask)

% D. Cheyne, Sept, 2023 - new channel set editor.
% - returns one of a number of predefined channel sets or a custom channel set. 

defaultSets = {'Custom';'MEG Sensors';'ADC Channels';'Trigger Channel';'Digital Channels';'None'};

maxChannels = 32;       % max channels to plot at once 


channelTypes = [header.channel.sensorType];
nchans = length(channelTypes);
if ~exist('badChannelMask','var')
    badChannelMask = zeros(1,nchans);
end

longnames = {header.channel.name};   
channelNames = cleanChannelNames(longnames); 
for k=1:nchans
    if badChannelMask(k) == 1
        channelNames(k) = strcat('*',channelNames(k),'*');
    end
end

scrnsizes=get(0,'MonitorPosition');

fh = figure('color','white','name','CTF Channel Selector','MenuBar','none',...
    'numbertitle','off', 'Position', [scrnsizes(1,4)/2 scrnsizes(1,4)/2  700 900], 'closeRequestFcn', @cancel_button_callBack);
if ispc
    movegui(fh,'center');
end

displaylistbox=uicontrol('Style','Listbox','FontSize',10,'Units','Normalized',...
    'Position',[0.05 0.06 0.4 0.37],'HorizontalAlignment',...
    'Center','BackgroundColor','White','max',10000,'Callback',@displaylistbox_callback);

hidelistbox=uicontrol('Style','Listbox','FontSize',10,'Units','Normalized',...
    'Position',[0.55 0.06 0.4 0.37],'HorizontalAlignment',...
    'Center','BackgroundColor','White','max',10000,'Callback',@hidelistbox_callback);

uicontrol('style','text','fontsize',12,'units','normalized',...
    'position',[0.06 0.51 0.2 0.05],'string','Channel Sets:','HorizontalAlignment',...
    'left','backgroundcolor','white','FontWeight','bold');

includetext=uicontrol('style','text','fontsize',12,'units','normalized',...
    'position',[0.05 0.43 0.25 0.03],'string','Included Channels:','HorizontalAlignment',...
    'left','backgroundcolor','white','FontWeight','bold');
excludetext=uicontrol('style','text','fontsize',12,'units','normalized',...
    'position',[0.55 0.43 0.25 0.03],'string','Excluded Channels:','HorizontalAlignment',...
    'left','backgroundcolor','white','FontWeight','bold');


uicontrol('style','text','fontsize',12,'units','normalized',...
    'position',[0.05 0.03 0.3 0.03],'string','Bad Channels indicated by: * *','HorizontalAlignment',...
    'left','backgroundcolor','white');

%Apply button
uicontrol('Style','PushButton','FontSize',13,'Units','Normalized','Position',...
    [0.57 0.5 0.15 0.06],'String','Apply','HorizontalAlignment','Center',...
    'BackgroundColor',[0.99,0.64,0.3],'ForegroundColor','white','Callback',@apply_button_callback);

    function apply_button_callback(~,~)
        selectedChannels=find(channelExcludeFlags == 0);

        if numel(selectedChannels) > maxChannels
            s = sprintf('Setting display to more than %d channels will slow down plotting. Proceed?',numel(selectedChannels));        
            r = questdlg(s,'Data Editor','Yes','No','No');
            if strcmp(r,'No')
                return;
            end
        end

        delete(fh);
    end

%Cancel button

uicontrol('Style','PushButton','FontSize',13,'units','normalized','Position',...
    [0.78 0.5 0.15 0.06],'String','Cancel',...
    'BackgroundColor','white','FontSize',13,'ForegroundColor','black','callback',@cancel_button_callBack);
              
    function cancel_button_callBack(~,~)
        selectedChannels = [];
        delete(fh);
    end

%title
uicontrol('style','text','units','normalized','position',[0.1 0.95 0.8 0.04],...
        'String','Channel Selector','FontSize',20,'ForegroundColor',[0.93,0.6,0.2], 'HorizontalAlignment','center','BackGroundColor', 'white');


%%%%%%%%%%%%
% init flags

goodChans = {};
badChans = {};

channelExcludeFlags = ones(numel(channelNames),1);
channelExcludeFlags(oldSelected) = 0;               % flag previous selected channels 


function displaylistbox_callback(src,~)  
    % get selected rows    
    xx = get(src,'value');
end

function hidelistbox_callback(~,~)
end

right_arrow=draw_rightarrow;
uicontrol('Style','pushbutton','FontSize',10,'Units','Normalized',...
    'Position',[0.46 0.3 0.08 0.05],'CData',right_arrow,'HorizontalAlignment',...
    'Center','BackgroundColor','White','Callback',@tohidearrow_callback);
left_arrow=draw_leftarrow;
uicontrol('Style','pushbutton','FontSize',10,'Units','Normalized',...
    'Position',[0.46 0.2 0.08 0.05],'CData',left_arrow,'HorizontalAlignment',...
    'Center','BackgroundColor','White','Callback',@todisplayarrow_callback);

    function tohidearrow_callback(~,~)
        idx=get(displaylistbox,'value');
        list = get(displaylistbox,'String');
        if isempty(list)
            return;
        end
        selected = list(idx,:);
        for i=1:size(selected,1)
            a = selected(i); 
            idx = find(strcmp(a,channelNames));
            channelExcludeFlags(idx) = 1;
        end
        updateChannelLists;

    end

    function todisplayarrow_callback(~,~)
        idx=get(hidelistbox,'value');
        list = get(hidelistbox,'String');
        if isempty(list)
            return;
        end
        selected = list(idx,:);
        for i=1:size(selected,1)
            a = selected(i);
            idx = find(strcmp(deblank(a),channelNames));
            channelExcludeFlags(idx) = 0;
        end
        updateChannelLists;
                
    end



function updateChannelLists
    goodChans = {};
    badChans = {};
    badChanCount = 0;
    goodChanCount = 0;
    for i=1:size(channelExcludeFlags,1)
        if channelExcludeFlags(i) == 1
            badChanCount = badChanCount + 1;
            badChans(badChanCount) = channelNames(i);
        else
            goodChanCount = goodChanCount + 1;
            goodChans(goodChanCount) = channelNames(i);
        end                
    end
    
    
    % make sure we are setting list beyond range.
    
    set(displaylistbox,'String',goodChans);
    set(hidelistbox,'String',badChans);   
     
    if ~isempty(goodChans)
        idx = get(displaylistbox,'value');
        if idx(end) > size(goodChans,2) && size(goodChans,2) > 0
            set(displaylistbox,'value',size(goodChans,2));
        end
    end
    
    if ~isempty(badChans)     
        idx = get(hidelistbox,'value');
        if idx(end) > size(badChans,2) && size(badChans,2) > 0
            set(hidelistbox,'value',size(badChans,2));
        end     
    end
        
    s = sprintf('Included channels (%d):',goodChanCount);
    set(includetext,'string',s);

    s = sprintf('Excluded channels (%d):',badChanCount);
    set(excludetext,'string',s);
    
end

% shortcut to default channnels...
uicontrol('style','popup','units','normalized',...
    'position',[0.05 0.45 0.35 0.06],'String',defaultSets, 'Backgroundcolor','white','fontsize',12,...
    'value',1,'callback',@channel_popup_callback);

    function channel_popup_callback(src,~)
        menuSelect=get(src,'value');
        updateMenuSelection(menuSelect);
    end

function updateMenuSelection(idx)

    % set exclude flag for all
    nchans = numel(channelNames);
    for i=1:nchans 
        channelExcludeFlags(i) = 1;
    end

    switch idx
        case 1  % All channels
            channelExcludeFlags(:) = 1;   
            channelExcludeFlags(oldSelected) = 0;
        case 2  % MEG
            for i=1:nchans 
                if (channelTypes(i) == 5); channelExcludeFlags(i) = 0;
                end
            end
        case 3  % ADC channels
            for i=1:nchans  
                if (channelTypes(i) == 18); channelExcludeFlags(i) = 0; 
                end                
            end    
        case 4  % digital (PPT) channels (CTF only)
            for i=1:nchans  
                if (channelTypes(i) == 20); channelExcludeFlags(i) = 0; 
                end                
            end         
        case 5  % trigger channel (CTF only)
            for i=1:nchans  
                if (channelTypes(i) == 19); channelExcludeFlags(i) = 0; 
                end                
            end              
        case 6  % trigger channel (CTF only)
            for i=1:nchans  
                channelExcludeFlags(i) = 1;               
            end


    end

   %update listbox
   updateChannelLists;     
   
end

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

% end of initialization
updateChannelLists;
    
% PAUSES MATLAB
uiwait(gcf);
end
