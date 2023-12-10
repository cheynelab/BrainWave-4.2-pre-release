% strip sensor version number from channel names for CTF
% replaces bw_truncateChannelNames - returns cellstr
% D. Cheyne, Dec, 2023
function channelNames = bw_cleanChannelNames(names) 
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