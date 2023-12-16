function meg_fft(dsName,channelName)

    averageSpectra = true;
    usePowerOf2 = false;
   
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
    
    % only reads trial 1 of N
    [tvec, data] = bw_CTFGetChannelData(dsName, channelName);
    
    timeVec = tvec(1:N);
    
         
    figure(99)

    % plot fft 
    for k=1:numTrials
          
        d = data(1:N,k) * 1e15;

        % remove offset
        offset = mean(d);
        d = d - offset;

        % plot data segment
        plot(timeVec,d,'blue');

        % compute fft with hanning window
%         win = hann(N);   
        win = tukeywin(N,0.5);
        
        y = fft(d.* win);   
        
        norm = sqrt(1.0/(N*fs));
        amp(:,k) = 2.0 * abs(y) .* norm;
%          pow =  1/(N*fs) * 2 * abs(y).^2;            % units are T^2 / Hz        
%          amp(:,k) = sqrt(pow);                        % units are T / sqrt(Hz) 
               
        hold on;
    end
   
    freq = 0:fs/N:fs/2;
        
    figure(98);
 
    if averageSpectra
        amp = mean(amp,2);
    end

    loglog(freq, amp(1:length(freq),:), 'blue');

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

    hold off
    

    

 
end
