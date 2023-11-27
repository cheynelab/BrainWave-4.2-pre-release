%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function bw_editCTFMarkers( markerFileName )
% GUI to select a marker from CTF Marker.mrk file
%
% input:   name of a CTF MarkerFile (e.g., dsName/MarkerFile.mrk)
%
% returns: latencies and label for selected marker
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bw_editCTFMarkers( dsName )
 
    markerFileName = strcat(dsName, filesep, 'MarkerFile.mrk');
    
    if ~exist(markerFileName,'file')
        errordlg('No marker file exists for this dataset.');
        return;
    end
    
    scrnsizes=get(0,'MonitorPosition');

    fg=figure('color','white','name','Edit Markers','numbertitle','off','menubar','none','position',[300 (scrnsizes(1,4)-300) 500 250]);   
    
    Event_popup = uicontrol('style','popup','units','normalized',...
    'position',[0.05 0.55 0.35 0.1],'String','No events','Backgroundcolor','white','fontsize',12,'value',1,'callback',@event_popup_callback);

        function event_popup_callback(src,~)
            menu_select=get(src,'value');
            markerTimes = markerData{menu_select};
            latencies = markerTimes(:,2);
            str = sprintf('Number of Events = %d\n',length(latencies));   
            set(MarkerCountText,'String',str);

        end
   
    uicontrol('style','text','units','normalized','HorizontalAlignment','Left',...
    'position',[0.08 0.7 0.6 0.1],'String','Marker Name:','Backgroundcolor','white','fontsize',12);

    MarkerCountText = uicontrol('style','text','units','normalized','HorizontalAlignment','Left',...
    'position',[0.08 0.25 0.6 0.1],'String','Number of Events =','Backgroundcolor','white','fontsize',12);
   
    uicontrol('style','pushbutton','units','normalized','position',...
        [0.6 0.4 0.25 0.2],'string','Delete Marker','backgroundcolor','white',...
        'foregroundcolor','black','callback',@delete_callback);
    
    uicontrol('style','pushbutton','units','normalized','position',...
        [0.6 0.7 0.25 0.2],'string','Rename Marker','backgroundcolor','white',...
        'foregroundcolor','black','callback',@rename_callback);

    [names, markerData] = bw_readCTFMarkerFile(markerFileName);     
    set(Event_popup,'string',names);

    % initialize list to first marker
    set(Event_popup,'String',names,'value',1);  
    markerTimes = markerData{1};
    latencies = markerTimes(:,2);
    str = sprintf('Number of Events = %d\n',numel(latencies));   
    set(MarkerCountText,'String',str);
          
    function rename_callback(~,~)
        idx = get(Event_popup,'val'); 
        s = names(idx);
        output = inputdlg('Enter Marker Name','Edit Markers', [1 50],s);
        newName = output{1};
              
        if ~isempty( find( strcmp(newName, names) == 1))
            warndlg('A Marker with this name already exists ...');
            return;
        end
        names(idx) = {newName};
        set(Event_popup,'string',names);
    end

       
    function delete_callback(~,~)          
        idx = get(Event_popup,'val'); 
        names(idx) = [];
        markerData(idx) = [];
        set(Event_popup,'string',names);  
        if idx > numel(names)
            idx = idx-1;
        end      
        set(Event_popup,'value',idx-1);
    end

    uicontrol('style','pushbutton','units','normalized','position',...
        [0.75 0.1 0.2 0.15],'string','OK','backgroundcolor','white',...
        'foregroundcolor','blue','callback',@ok_callback);
    
    function ok_callback(~,~)
        r = questdlg('Save Changes (these cannot be undone)','Edit Markers', 'Yes', 'No','No');
        if strcmp(r,'Yes')
            for k=1:numel(names)
                % note structure passed to bw_writeNewMarkerFile is slightly different that       
                % that returned from readMarkerData including shift of trial numbering ...
                newMarkerData(k).ch_name = char(names(k)); 
                markerTimes = markerData{k};
                newMarkerData(k).trials = markerTimes(:,1) - 1;
                newMarkerData(k).latencies = markerTimes(:,2);
            end
            bw_writeNewMarkerFile(dsName, newMarkerData);   
        end

        uiresume(gcf);
    end

    uicontrol('style','pushbutton','units','normalized','position',...
        [0.5 0.1 0.2 0.15],'string','Cancel','backgroundcolor','white',...
        'callback',@cancel_callback);

    function cancel_callback(~,~)
        uiresume(gcf);
    end

        % % note bw_readCTFMarkerFile adds 1 to the trial numbers ...
        % [markerNames, markerData] = bw_readCTFMarkerFile(markerFileName);     
        % 
        % % for each marker correct trial number and/or latency
        % for k=1:numel(markerData)            
        %     newMarkerData(k).ch_name = char(markerNames{k}); 
        % 
        %     markerTimes = markerData{k};
        %     trials = markerTimes(:,1);      % numbering starts at 1
        %     latencies = markerTimes(:,2);
        % 
        %     % remove deleted trials and renumber
        %     if ~isempty(badTrialIdx)     
        %         fprintf('correcting marker trial numbers...');
        %         trials(badTrialIdx) = [];       % compress lists to valid trials only
        %         latencies(badTrialIdx) = [];
        % 
        %         % renumber the trials - tricky!               
        %         for j=1:length(trials)
        %             n = trials(j);
        %             idx = find(badTrialIdx < n);   % which trials preceding this one deleted?
        %             shift = length(idx);           % how many? shift = 0 if idx = []
        %             trials(j) = n - shift;                
        %         end
        %     end
        % 
        %     % correct latencies if time zero was shifted forward
        %     if t1 > 0.0 
        %         fprintf('correcting marker latencies (subtracting %.4f seconds)\n', t1);
        %         latencies = latencies - t1;
        %     end
        % 
        %     newMarkerData(k).trials = trials - 1;      % correct back to base zero before writing
        %     newMarkerData(k).latencies = latencies;          
        % end
        % 
        % % write corrected markerFile to the new dataset (overwrite
        % bw_writeNewMarkerFile(dsName, newMarkerData);     
        

    
    
    %%PAUSES MATLAB
    uiwait(gcf);
    %%CLOSES GUI
    close(fg);   
    
    
end
