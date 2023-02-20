%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write vertices and faces to .vtk format
% D. Cheyne Nov, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bw_save_surface_vtk(vertices, faces, filename )

    fid=fopen(filename,'w');
    if fid == 0
        return;
    end

    faces = faces - 1;  % faces are base 0 in vtk files.
    npts = size(vertices,1);

    fprintf('Saving %d points to VTK file %s\n', npts, filename);
    fprintf(fid,'# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET POLYDATA\n');
    % write vertices
    fprintf(fid,'POINTS %d float\n',npts);

    for k=1:npts
        fprintf(fid,'%3.7f %3.7f %3.7f\n', vertices(k,1), vertices(k,2), vertices(k,3));
    end
    % write faces
    nfaces=size(faces,1);
    fprintf(fid,'POLYGONS %d %d\n',nfaces,4*nfaces);
    for k=1:nfaces
       fprintf(fid,'3 %d %d %d\n',faces(k,1), faces(k,2), faces(k,3));
    end
    fclose(fid);        


end