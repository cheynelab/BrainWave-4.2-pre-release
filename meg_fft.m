function meg_fft(dsName,channelNames)

    averageSpectra = true;
    usePowerOf2 = false;
    plotData = false;

    header = bw_CTFGetHeader(dsName);
    if isempty(header)
        return;
    end
    fs = header.sampleRate;
    numTrials = header.numTrials;
    
    nsamples = header.numSamples;
     
    % get power of two samples
    
    if usePowerOf2
        N = 2^nextpow2(nsamples);
        if N > header.numSamples
            N = N/2;      
        end       
    else
        N = nsamples;
    end
       
    numChannelsToPlot = numel(channelNames);
              
    % compute fft with window
    % CTF uses 50% cosine window?
    %         win = hann(N);   
    win = tukeywin(N,0.5);
    norm = sqrt(1.0/(N*fs));  


    % create window 
    pos1 = get(gcf,'Position');
    [~,n,e] = fileparts(dsName);
    s = sprintf('FFT: %s',[n e]);
    fh = figure('numbertitle','off','Name',s);
    pos = [pos1(1)+pos1(3) pos1(2) 600 500];
    set(fh,'Position',pos);

    uicontrol('Style','radiobutton',...
        'fontsize',12,...
        'units', 'normalized',...
        'Position',[0.15 0.93 0.2 0.05],...
        'String','Plot Average',...
        'Value',averageSpectra,...
        'Callback',@averageCallback);          

    draw;


    function draw

        for j=1:numChannelsToPlot
            
            chanName = char(channelNames{j});
            % returns all trials for one channel
            [tvec, data] = bw_CTFGetChannelData(dsName, chanName);
          
            timeVec = tvec(1:N);        
        
            % plot fft 
            for k=1:numTrials
                  
                d = data(1:N,k) * 1e15;
        
                % remove offset
                offset = mean(d);
                d = d - offset;
        
                % plot data segment
                if k==1 && plotData
                  figure(98)
                  plot(timeVec,d);
                  hold on;
                end
        
                y = fft(d.* win);               
                % scale to femtTesla / sqrtHz
                amp(:,k) = 2.0 * abs(y) .* norm;
                       
            end
           
            freq = 0:fs/N:fs/2;
               
            if averageSpectra
                amp = mean(amp,2);
            end
        
            loglog(freq, amp(1:length(freq),:));
            hold on;
    
        end
    
        ylim([0.001 1000]);
        xlim([0 1000]);
        ax = gca;
        
        ax.YAxis.TickLabels = compose('%g', ax.YAxis.TickValues);
        ylabel('Magnitude (fT / sqrt(Hz) )');
        
        ax.XAxis.TickLabels = compose('%g', ax.XAxis.TickValues);
        xlabel('Frequency (Hz)');
    
        grid on;
        ax. GridColor = 'black';
        ax. GridAlpha = 0.4;
    
        legend(channelNames);
    
        hold off
    
    end
    
    function averageCallback(~,~)
        averageSpectra = ~averageSpectra;
        draw;
    end
 
end
