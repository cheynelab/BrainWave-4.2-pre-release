function brainwave
    % start brain wave
    
    % version 4.2 - this function is called from parent folder (e.g., BrainWave_Toolbox)
    
    s = which('brainwave');
    if isempty(s)
        fprintf('brainwave not found. Make sure CheyneLab_Toolbox is installed and in your Matlab path\n');
        return;
    end
    
    [bw_path,~,~] = fileparts(s);
    addpath(bw_path);
    bw_main_menu      % this adds all other paths needed
    
end