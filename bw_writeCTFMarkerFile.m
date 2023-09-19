%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function bw_writeCTFMarkerFile(dsName,trig)
%   Writes a markerFile.mrk for the given input 'trig'
%   Returns on error if a markerFile.mrk exists
%   INPUT:
%         dsName = ctf dataset name
%         trig = 1 X N(num of trigs) structure
%                trig(1:N).ch_name 
%                         .onset_idx %sample index of trigger onsets                 
%                         .times     %dataset times of triggers onsets
%
%  pferrari@meadowlandshospital.org, Aug2012
%
% this is a modified version for eventMarker that overwrites the existing markerFile
% use bw_write_MarkerFile to flag overwrite.
%
function result = bw_writeCTFMarkerFile(dsName,trig)

result = 0;

no_trigs=numel(trig);

filepath=strcat(dsName, filesep, 'MarkerFile.mrk');
          
fprintf('writing marker file %s\n', filepath);

fid = fopen(filepath,'w','n');
fprintf(fid,'PATH OF DATASET:\n');
fprintf(fid,'%s\n\n\n',dsName);
fprintf(fid,'NUMBER OF MARKERS:\n');
fprintf(fid,'%g\n\n\n',no_trigs);

for i = 1:no_trigs
    
    fprintf(fid,'CLASSGROUPID:\n');
    fprintf(fid,'3\n');
    fprintf(fid,'NAME:\n');
    fprintf(fid,'%s\n',trig(i).ch_name);
    fprintf(fid,'COMMENT:\n\n');
    fprintf(fid,'COLOR:\n');
    fprintf(fid,'blue\n');
    fprintf(fid,'EDITABLE:\n');
    fprintf(fid,'Yes\n');
    fprintf(fid,'CLASSID:\n');
    fprintf(fid,'%g\n',i);
    fprintf(fid,'NUMBER OF SAMPLES:\n');
    fprintf(fid,'%g\n',length(trig(i).times));
    fprintf(fid,'LIST OF SAMPLES:\n');
    fprintf(fid,'TRIAL NUMBER\t\tTIME FROM SYNC POINT (in seconds)\n');
    for t = 1:length(trig(i).times)-1
        fprintf(fid,'                  %+g\t\t\t\t               %+0.6f\n',0,trig(i).times(t));
    end
    fprintf(fid,'                  %+g\t\t\t\t               %+0.6f\n\n\n',0,trig(i).times(end));
end

fclose(fid);

result = 1;

end






