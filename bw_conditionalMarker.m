%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function [latencies] = bw_readCTFMarkers( markerFileName )
% GUI to select a marker from CTF Marker.mrk file
%
% input:   name of a CTF MarkerFile (e.g., dsName/MarkerFile.mrk)
%
% returns: latencies and label for selected marker
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [latencies] = bw_conditionalMarker( markerFileName )
 

    latencies = [];
    selectedMarker = 1;
                
            
    include_latencies = [];       
    includeWindowStart = -1.0;
    includeWindowEnd = 1.0;

    exclude_latencies = [];       
    excludeWindowStart = -1.0;
    excludeWindowEnd = 1.0;



    if ~exist(markerFileName,'file')
        errordlg('No marker file exists yet. Create or import latencies then save events as markers.');
        return;
    end
    
    scrnsizes=get(0,'MonitorPosition');

    fg=figure('color','white','name','Create Conditional Event','numbertitle','off','menubar','none','position',[300 (scrnsizes(1,4)-300) 650 450]);
    
    uicontrol('style','text','units','normalized','HorizontalAlignment','Left','position',[0.1 0.88 0.9 0.05],...
    'string','Select Sync Marker (t = 0.0):','Backgroundcolor','white','fontsize',14,'FontWeight','bold');

    sync_popup = uicontrol('style','popup','units','normalized',...
    'position',[0.09 0.75 0.3 0.1],'String','No events','Backgroundcolor','white','fontsize',12,'value',1,'callback',@sync_popup_callback);
        
        function sync_popup_callback(src,~)
            selectedMarker=get(src,'value');

            t = trials{selectedMarker};
            latencies = t(:,2);
            str = sprintf('Number of Events = %d\n',length(latencies));   
            set(eventCount,'String',str);
        end

    str = sprintf('Number of Markers = %d\n',length(latencies));   
    eventCount = uicontrol('style','text','units','normalized','HorizontalAlignment','Left',...
    'position',[0.5 0.8 0.6 0.05],'String',str,'Backgroundcolor','white','fontsize',12);
   

    % get marker data
    [names, trials] = bw_readCTFMarkerFile( markerFileName );

    if isempty(names)
        errordlg('No Markers defined');
        uiresume(gcf);
    end

    % initialize list to first marker
    set(sync_popup,'String',names,'value',1);  
    t = trials{1};
    latencies = t(:,2);
    str = sprintf('Number of Events = %d\n',length(latencies));   
    set(eventCount,'String',str);
          

   %%%%% inclusion controls

    uicontrol('style','text','units','normalized','HorizontalAlignment','Left','position',[0.1 0.66 0.9 0.05],...
    'string','Select Inclusion Mask (only include sync events within this marker time window)','Backgroundcolor','white','fontsize',14,'FontWeight','bold');

    uicontrol('style','text','units','normalized','HorizontalAlignment','Left','position',[0.45 0.6 0.15 0.05],...
    'string','Window Start (s):','Backgroundcolor','white','fontsize',12);

    uicontrol('style','text','units','normalized','HorizontalAlignment','Left','position',[0.7 0.6 0.15 0.05],...
    'string','Window End (s):','Backgroundcolor','white','fontsize',12);

    menuText = ['None';names];
    uicontrol('style','popup','units','normalized',...
    'position',[0.09 0.5 0.3 0.1],'String',menuText,'Backgroundcolor','white','fontsize',12,'value',1,'callback',@include_popup_callback);

        function include_popup_callback(src,~)
            menu_select=get(src,'value');
            if menu_select > 1
                t = trials{menu_select-1};
                include_latencies = t(:,2);
            else
                include_latencies = [];
            end
        end

    
    s = sprintf('%.4g', includeWindowStart);
    includeEditStart = uicontrol('style','edit','units','normalized','HorizontalAlignment','center','position',[0.45 0.55 0.15 0.05],...
        'string',s,'Backgroundcolor','white','fontsize',12, 'callback',@includeStart_callback);

    function includeStart_callback(src,~)
        s = get(src,'string');
        includeWindowStart = str2double(s);
    end

    s = sprintf('%.4g', includeWindowEnd);
    includeEditEnd = uicontrol('style','edit','units','normalized','HorizontalAlignment','center','position',[0.7 0.55 0.15 0.05],...
    'string',s,'Backgroundcolor','white','fontsize',12, 'callback',@includeEnd_callback);

    function includeEnd_callback(src,~)
        s = get(src,'string');
        includeWindowEnd = str2double(s);    
    end

    %%%%% exclusion controls

    uicontrol('style','text','units','normalized','HorizontalAlignment','Left','position',[0.1 0.4 0.9 0.05],...
    'string','Select Exclusion Mask (exclude sync events within this marker time window)','Backgroundcolor','white','fontsize',14,'FontWeight','bold');

    uicontrol('style','text','units','normalized','HorizontalAlignment','Left','position',[0.45 0.34 0.15 0.05],...
    'string','Window Start (s):','Backgroundcolor','white','fontsize',12);

    uicontrol('style','text','units','normalized','HorizontalAlignment','Left','position',[0.7 0.34 0.15 0.05],...
    'string','Window End (s):','Backgroundcolor','white','fontsize',12);
    
    uicontrol('style','popup','units','normalized',...
    'position',[0.09 0.24 0.3 0.1],'String',menuText,'Backgroundcolor','white','fontsize',12,'value',1,'callback',@exclude_popup_callback);

        function exclude_popup_callback(src,~)
            menu_select=get(src,'value');
            if menu_select > 1
                t = trials{menu_select-1};
                exclude_latencies = t(:,2);
            else
                exclude_latencies = [];
            end
        end

    s = sprintf('%.4g', excludeWindowStart);
    excludeEditStart = uicontrol('style','edit','units','normalized','HorizontalAlignment','center','position',[0.45 0.29 0.15 0.05],...
        'string',s,'Backgroundcolor','white','fontsize',12, 'callback',@excludeStart_callback);

    function excludeStart_callback(src,~)
        s = get(src,'string');
        excludeWindowStart = str2double(s);
    end

    s = sprintf('%.4g', excludeWindowEnd);
    excludeEditEnd = uicontrol('style','edit','units','normalized','HorizontalAlignment','center','position',[0.7 0.29 0.15 0.05],...
    'string',s,'Backgroundcolor','white','fontsize',12, 'callback',@excludeEnd_callback);

    function excludeEnd_callback(src,~)
        s = get(src,'string');
        excludeWindowEnd = str2double(s);    
    end


    function updateEvents

        % reset to original count
        t = trials{selectedMarker};
        latencies = t(:,2);
        totalEvents = length(latencies);
        numEvents = length(latencies);

        if numEvents > 0 && ~isempty(include_latencies)        
            idx = [];
            for k=1:numEvents  
                latency = latencies(k);
                for j=1:numel(include_latencies)
                    wStart = include_latencies(j) + includeWindowStart;               
                    wEnd = include_latencies(j) + includeWindowEnd;
                    if latency > wStart && latency < wEnd
                        idx(end+1) = k;
                    end
                end
            end    
            latencies = latencies(idx);
        end

        % apply exclude after include?

        numEvents = length(latencies);
        if numEvents > 0 && ~isempty(exclude_latencies)        
            mask = ones(1,numEvents);
            for k=1:numEvents  
                latency = latencies(k);
                for j=1:numel(exclude_latencies)
                    wStart = exclude_latencies(j) + excludeWindowStart;               
                    wEnd = exclude_latencies(j) + excludeWindowEnd;
                    if latency > wStart && latency < wEnd
                        mask(k) = 0;
                    end
                end
            end    
            idx = find(mask == 1);
            latencies = latencies(idx);
        end

        s = sprintf('Including %d of %d events', length(latencies), totalEvents);
        msgbox(s)

        str = sprintf('Number of Events = %d\n',length(latencies));     
        set(eventCount,'String',str);

    end

    %%%% 
    
    uicontrol('style','pushbutton','units','normalized','position',...
        [0.15 0.1 0.2 0.1],'string','Create Event','backgroundcolor','white',...
        'foregroundcolor','blue','callback',@ok_callback);
    
    function ok_callback(~,~)
        uiresume(gcf);
    end
    
    uicontrol('style','pushbutton','units','normalized','position',...
        [0.43 0.1 0.2 0.1],'string','Update Events','backgroundcolor','white',...
        'callback',@update_callback);
    
    function update_callback(~,~)
        s = get(includeEditStart,'string');
        includeWindowStart = str2double(s);
        s = get(includeEditEnd,'string');
        includeWindowEnd = str2double(s);

        s = get(excludeEditStart,'string');
        excludeWindowStart = str2double(s);
        s = get(excludeEditEnd,'string');
        excludeWindowEnd = str2double(s);
        
        updateEvents;

    end

    uicontrol('style','pushbutton','units','normalized','position',...
        [0.7 0.1 0.2 0.1],'string','Cancel','backgroundcolor','white',...
        'foregroundcolor','black','callback',@cancel_callback);
    
    function cancel_callback(~,~)
        latencies = [];
        uiresume(gcf);
    end
    
    
    %%PAUSES MATLAB
    uiwait(gcf);
    %%CLOSES GUI
    close(fg);   
    
    
end
