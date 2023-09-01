function [selectedChannels, menuSelect] = bw_channelSelector(channelNames, channelTypes, oldSelected, oldMenuSelect)


defaultPrefsFile = sprintf('channelsets.mat');

channelLists = [];
numDefaultChanList = 5;

scrnsizes=get(0,'MonitorPosition');

fh = figure('color','white','name','CTF Channel Selector','MenuBar','none',...
    'numbertitle','off', 'Position', [scrnsizes(1,4)/2 scrnsizes(1,4)/2  700 700], 'closeRequestFcn', @cancel_button_callBack);
if ispc
    movegui(fh,'center');
end

displaylistbox=uicontrol('Style','Listbox','FontSize',10,'Units','Normalized',...
    'Position',[0.05 0.05 0.4 0.35],'HorizontalAlignment',...
    'Center','BackgroundColor','White','max',10000,'Callback',@displaylistbox_callback);

hidelistbox=uicontrol('Style','Listbox','FontSize',10,'Units','Normalized',...
    'Position',[0.55 0.05 0.4 0.35],'HorizontalAlignment',...
    'Center','BackgroundColor','White','max',10000,'Callback',@hidelistbox_callback);

uicontrol('style','text','fontsize',12,'units','normalized',...
    'position',[0.06 0.5 0.2 0.05],'string','Channel List:','HorizontalAlignment',...
    'left','backgroundcolor','white');

includetext=uicontrol('style','text','fontsize',12,'units','normalized',...
    'position',[0.05 0.42 0.25 0.02],'string','Included Channels:','HorizontalAlignment',...
    'left','backgroundcolor','white','FontWeight','bold');
excludetext=uicontrol('style','text','fontsize',12,'units','normalized',...
    'position',[0.55 0.42 0.25 0.02],'string','Excluded Channels:','HorizontalAlignment',...
    'left','backgroundcolor','white','FontWeight','bold');

%%%%%%%%%%%%
% init flags

goodChans = {};
badChans = {};

channelExcludeFlags = ones(numel(channelNames),1);
channelExcludeFlags(oldSelected) = 0;
menuSelect = oldMenuSelect;


% not used...
function displaylistbox_callback(src,~)       

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

uicontrol('units','normalized','Position',[0.55 0.465 0.19 0.06],'String','Save ChannelSet',...
              'BackgroundColor','white','FontSize',13,'ForegroundColor','black','callback',@save_button_callBack);
              
function save_button_callBack(~,~)
    
    input = inputdlg({'Name for Channelset'},'Save ChannelSet',[1 50],{'myList'});
    if isempty(input)
        return;
    end   
    listName = input{1}; 
    
    % add to current channelLists structure and save in prefs
    listNo = size(channelLists,2) + 1;
    channelLists(listNo).name = listName;
    channelLists(listNo).list = get(displaylistbox,'String');  % save list of channel names           
    
    prefs.channelLists = channelLists;
    
    fprintf('saving current settings to file %s\n', defaultPrefsFile);
    save(defaultPrefsFile, '-struct', 'prefs')    
        
    % add to menu
    buildChannelMenu;
    newlist = get(channel_popup,'String');
    menuSelect = numel(newlist);
    
    set(channel_popup,'value',menuSelect);
    updateMenuSelection(menuSelect);      

end

delete_button = uicontrol('units','normalized','Position',[0.76 0.465 0.19 0.06],'String','Delete ChannelSet','enable','off',...
              'BackgroundColor','white','FontSize',13,'ForegroundColor','black','callback',@delete_button_callBack);
              
function delete_button_callBack(~,~)

    idx = get(channel_popup,'value');  

    listNo = idx - numDefaultChanList;
    s = sprintf('Delete channel set [%s]?',char(channelLists(listNo).name));
    response = questdlg(s,'Channel Selector','Yes','No','No');
    if strcmp(response,'No')
        return;
    end
    
    % delete list
    if size(channelLists,2) == 1
        channelLists = [];
    else    
        channelLists(listNo) = [];
    end

    prefs.channelLists = channelLists;
    
    fprintf('saving current settings to file %s\n', defaultPrefsFile)
    save(defaultPrefsFile, '-struct', 'prefs')    
    
    buildChannelMenu;
    
    % switch to default
    menuSelect = 1;
    set(channel_popup,'value',menuSelect);
    updateMenuSelection(menuSelect);
           
end

%Apply button
uicontrol('Style','PushButton','FontSize',13,'Units','Normalized','Position',...
    [0.57 0.55 0.15 0.06],'String','Apply','HorizontalAlignment','Center',...
    'BackgroundColor',[0.99,0.64,0.3],'ForegroundColor','white','Callback',@apply_button_callback);

    function apply_button_callback(~,~)
        selectedChannels=find(channelExcludeFlags == 0);
        menuSelect = get(channel_popup,'value')
        delete(fh);
    end

%Cancel button

uicontrol('units','normalized','Position',[0.78 0.55 0.15 0.06],'String','Cancel',...
              'BackgroundColor','white','FontSize',13,'ForegroundColor','black','callback',@cancel_button_callBack);
              
    function cancel_button_callBack(~,~)
        selectedChannels = [];
        menuSelect = [];
        delete(fh);
    end

%title
uicontrol('style','text','units','normalized','position',[0.1 0.95 0.8 0.04],...
        'String','Channel Selector','FontSize',20,'ForegroundColor',[0.93,0.6,0.2], 'HorizontalAlignment','center','BackGroundColor', 'white');

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

channel_popup = uicontrol('style','popup','units','normalized',...
    'position',[0.05 0.45 0.35 0.06],'String',{}, 'Backgroundcolor','white','fontsize',12,...
    'value',1,'callback',@channel_popup_callback);

    function channel_popup_callback(src,~)
        menuSelect=get(src,'value');
        updateMenuSelection(menuSelect);
    end

function buildChannelMenu
    
    menuList = {'All Channels';'MEG Sensors';'ADC Channels';'Digital Channels';'Trigger Channel'};
    
    if ~isempty(channelLists)
        for j=1:size(channelLists,2)
            menuList{j+numDefaultChanList} = channelLists(j).name;
        end
    end
    
    if menuSelect > numel(menuList)
        menuSelect = numel(menuList);
    end
    
    set(channel_popup,'String',menuList);    
    set(channel_popup,'value',menuSelect);

end

function updateMenuSelection(idx)

    % set exclude flag for all
    nchans = numel(channelNames);
    for i=1:nchans 
        channelExcludeFlags(i) = 1;
    end

    switch idx
        case 1  % All channels
            channelExcludeFlags(:) = 0;   
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

    end

   %update listbox
   updateChannelLists;     
           
    if menuSelect > numDefaultChanList
       set(delete_button,'enable','on');
    else
       set(delete_button,'enable','off');
    end      
   
end

% end of initialization

buildChannelMenu;
updateChannelLists;
    
% PAUSES MATLAB
uiwait(gcf);
end
