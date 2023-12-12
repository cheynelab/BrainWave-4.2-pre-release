function meg_fft(dsName,channelName)

    stackTrials = 1;
    
    header = bw_CTFGetHeader(dsName);
    if isempty(header)
        return;
    end
    fs = header.sampleRate;
    numTrials = header.numTrials;
    
    nsamples = header.numSamples;
    
        
    % get power of two samples
    
%     N = 2^nextpow2(nsamples);
%     if N > header.numSamples
%         N = N/2;      
%     end       

    N = nsamples;
    
    % only reads trial 1 of N
    [tvec, data] = bw_CTFGetChannelData(dsName, channelName);

    timeVec = tvec(1:N);
    
    % plot data segment
    figure(1)
    plot(timeVec,data(1:N,:),'blue');
       
    % plot fft 
    for k=1:numTrials
          
        d = data(1:N,k) * 1e15;
        y = fft(d);
        pow =  2/(N*fs) * abs(y).^2;        % units are T^2 / Hz
        yt(1:N,k) = sqrt(pow);                  % units are T / sqrt(Hz) 
        
    end
    yy = mean(yt,2);
    
    freq = 0:fs/N:fs/2;
        
    figure(2);
    loglog(freq, yy(1:length(freq)), 'black');
    ylim([0.01 1000]);
    ax = gca;
    
    ax.YAxis.TickLabels = compose('%g', ax.YAxis.TickValues);
    ylabel('Magnitude (fT / sqrt(Hz) )');
    
    ax.XAxis.TickLabels = compose('%g', ax.XAxis.TickValues);
    xlabel('Frequency (Hz)');

    grid on;
    ax. GridColor = 'blue';
    ax. GridAlpha = 0.4;
    

    

 
end
