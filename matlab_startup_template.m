% add paths for Brainwave
disp('configuring BrainWave ...');
addpath('<path_to_BrainWave_Toolbox_folder>/BrainWave_Toolbox');
addpath('<path_to_matlab_folder>/spm12');

% setup FSL (MacOS and Linux only)
setenv('FSLDIR','<path_to_fsl_folder>/fsl');
setenv('FSLOUTPUTTYPE','NIFTI_GZ');
fsldir = getenv('FSLDIR');