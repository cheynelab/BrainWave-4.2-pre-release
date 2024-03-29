%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% October, 2019 Update (Merron Woodbury)
%% function meshFile = bw_importMeshFiles(mesh_dir, mri_dir, [mesh_path], [downSampleFactor))
%%
%% October, 2021 Update (D. Cheyne) - changed name to bw_importMeshFiles
%%  separate version for diffusionTools that combines all meshes into one .mat file
%%  also now saves the MEG to RAS transform matrix and thickness values
%%  
%% Feb, 2023 - moved new version back to BrainWave_Toolbox folder...
%%
%% Input:
%% mri_dir              - Subject's MRI directory
%% mesh_dir             - Subject's CIVET or Freesurfer directory
%% mesh_path            - (optional) specify path of saved mesh files
%% downSampleFactor     - (optional) downsample by factor 0.0 to 1.0 (e.g., 0.5 = reduce to 50% of original vertices)
%% 
%% Returns:       - filename of saved mesh file
%% 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function meshFile = bw_importMeshFiles(mesh_dir, mri_dir, mesh_path, downSampleFactor)

tpath=which('bw_importMeshFiles');
IMPORT_PATH=tpath(1:end-20);

% load rotation matrix from .mat file generated by mri2nii
%
if ~exist('downSampleFactor','var')
    downSampleFactor = 1.0;
end

if ~exist('mesh_path','var')
    mesh_path = [];
end

[~, n, ~] = fileparts(mri_dir);
postfix = regexp(n, '_MRI');
subjID = n(1:(postfix(end)-1));
matfile = fullfile(mri_dir, strcat(subjID,'.mat'));
tmat = load(matfile);       
M = tmat.M;             % voxel-to-head transformation matrix

% ** new save a copy of the head to RAS voxel transform that includes scaling from mm to cm
RAS_to_MEG = M * diag([0.1,0.1,0.1,1]);
MEG_to_RAS = inv(RAS_to_MEG);
RAS_to_MNI = [];

isFSaverage = 0; % for importing fsaverage surfaces ...

if isequal(mesh_dir, 0)
    return;
end
if isequal(mri_dir, 0)
    return;
end

% use folder structure to do initial assesment of CIVET vs Freesurfer
if exist(fullfile(mesh_dir, 'surfaces'), 'dir') && exist(fullfile(mesh_dir, 'transforms'), 'dir')
    isCIVET = true;
    importCIVETmeshes(mesh_dir, mri_dir);
elseif exist(fullfile(mesh_dir, 'surf'), 'dir') && exist(fullfile(mesh_dir, 'mri', 'transforms'), 'dir')
    isCIVET = false;
    importFSmeshes(mesh_dir, mri_dir);
end

% finish up by creating convex hull around gray surface mesh
if ~isempty(mesh_path)
    fprintf('Generating shape file from cortical surface ...\n');
    mesh = load(mesh_path);
    if isCIVET
        gray_mesh = mesh.gray_rsl;
    else
        gray_mesh = mesh.pial;
    end
    clear mesh;
    
    % 1) save surface as convex hull to shape file in MEG coordinates....
    headpts = gray_mesh.meg_vertices;   
    fprintf('read total of %d vertices  ...\n', size(headpts,1));

    % create convex hull from the current surface_points
    % which are mesh points in MEG coordinates if mesh has been loaded
    % restore full mesh after saving.

    % get vertices (in MEG coordinates) that make convex hull
    fprintf('Creating convex hull...\n');
    points = double(headpts);
    idx = convhulln(points);
    headpts = points(idx,:);    
   
    shapefile = fullfile(mri_dir, strcat(subjID, '_brainHull.shape'));

    % get updated array size!
    npts = size(headpts,1);
    fprintf('writing %d points to shape file (%s)\n', npts, shapefile);
    fid = fopen(shapefile,'w');
    fprintf(fid, '%d\n', npts);
    for k=1:npts
        fprintf(fid,'%6.2f %6.2f  %6.2f\n', headpts(k,1), headpts(k,2), headpts(k,3));
    end
    fclose(fid);  
end



% imports all CIVET meshes in chosen folder
 function importCIVETmeshes(mesh_dir, mri_dir)
        
        wbh = waitbar(0,'Reading CIVET meshfiles...');
        obj_files = dir(fullfile(mesh_dir, 'surfaces', '*.obj'));
        if isempty(obj_files)
            return;
        end
       
        % subject name
        names = cell(size(obj_files,1),1);
        for h=1:size(obj_files,1)
            names{h} = char(obj_files(h).name);
        end
        subj_name = names{1}(all(~diff(char(names(:)))));
        
        % different surfaces available
         for h=1:size(obj_files,1)
            surfs{h} = strrep(names(h), subj_name, '');
            surfs{h} = strrep(char(surfs{h}), '.obj', '');
        end
        
        surfaces{1, 1} = 'gray';
        surfaces{2, 1} = 'gray_rsl';
        surfaces{3, 1} = 'mid';
        surfaces{4, 1} = 'mid_rsl';
        surfaces{5, 1} = 'white';
        surfaces{6, 1} = 'white_rsl';
         
        % sort files into surface groups
        for h=1:size(surfs,2)
            name = char(surfs{h});
            
            % skip any files that are not followed by a hemisphere
            % specifier
            if strcmp(name((end-7):end), '_surface')
                continue;
            end

            if strncmp(name, 'gray_surface_left_', 18)
                surfaces{1, 2} = fullfile(mesh_dir, 'surfaces', char(obj_files(h).name));
            elseif strncmp(name, 'gray_surface_right_', 19)
                surfaces{1, 3} = fullfile(mesh_dir, 'surfaces', char(obj_files(h).name));
            elseif strncmp(name, 'gray_surface_rsl_left_', 22)
                surfaces{2, 2} = fullfile(mesh_dir, 'surfaces', char(obj_files(h).name));
            elseif strncmp(name, 'gray_surface_rsl_right_', 23)
                surfaces{2, 3} = fullfile(mesh_dir, 'surfaces', char(obj_files(h).name));
            elseif strncmp(name, 'mid_surface_left_', 17)
                surfaces{3, 2} = fullfile(mesh_dir, 'surfaces', char(obj_files(h).name));
            elseif strncmp(name, 'mid_surface_right_', 18)   
                surfaces{3, 3} = fullfile(mesh_dir, 'surfaces', char(obj_files(h).name));
            elseif strncmp(name, 'mid_surface_rsl_left_', 21)
                surfaces{4, 2} = fullfile(mesh_dir, 'surfaces', char(obj_files(h).name));
            elseif strncmp(name, 'mid_surface_rsl_right_', 22)
                surfaces{4, 3} = fullfile(mesh_dir, 'surfaces', char(obj_files(h).name));  
            elseif strncmp(name, 'white_surface_left_', 19)
                surfaces{5, 2} = fullfile(mesh_dir, 'surfaces', char(obj_files(h).name));
            elseif strncmp(name, 'white_surface_right_', 20)   
                surfaces{5, 3} = fullfile(mesh_dir, 'surfaces', char(obj_files(h).name));
            elseif strncmp(name, 'white_surface_rsl_left_', 23)
                surfaces{6, 2} = fullfile(mesh_dir, 'surfaces', char(obj_files(h).name));
            elseif strncmp(name, 'white_surface_rsl_right_', 24)
                surfaces{6, 3} = fullfile(mesh_dir, 'surfaces', char(obj_files(h).name));
            end

        end
        

        % create each surface
        num_surf = size(surfaces,1);
        for h=1:size(surfaces,1)
            
            leftMesh = surfaces{h,2};
            rightMesh = surfaces{h,3};
            
            % skip surfaces without the requisite files
            if isempty(leftMesh) || isempty(rightMesh)
                continue;
            end
            
            s = sprintf('Creating %s surface...\n', surfaces{h, 1});
            st = strrep(s, '_', ' ');
            waitbar((h/num_surf), wbh, st);
            
            % read one or more mesh files and combine all vertices and faces
            mesh = [];
            meshtype = [];
            idx = [];
            vertices = [];
            faces = [];
            normalized_vertices = [];
            RAS_to_MNI = [];
            thickness = [];

            for k=1:2

                if k == 1 
                    if isempty(leftMesh)
                        mesh.numLeftVertices = 0;
                        mesh.numLeftFaces = 0;
                        continue;
                    else
                        file = leftMesh;
                    end
                else
                    
                    if isempty(rightMesh)
                        mesh.numRightVertices = 0;
                        mesh.numRightFaces = 0;
                        continue;
                    else
                        file = rightMesh;
                    end
                end

                [meshtype, meshdata] = bw_readMeshFile(file);

                if isempty(meshdata)
                    fprintf('failed to read mesh file %s - check format...\n', file);
                    return;
                end
                
                % for CIVET meshes in non-native (MNI) coordinate space, we 
                % need transformation matrix from .xfm file to convert back to RAS coordinates

                % left hemisphere
                if k==1
                    % for CIVET mesh first save a copy of the original 
                    % vertices which are already in MNI coords                 
                    normalized_vertices = meshdata.vertices;

                % right hemisphere
                else
                    normalized_vertices = [normalized_vertices; meshdata.vertices];        
                end

                % get transform and thickness files - only one for both
                % hemispheres

                if (k==1)
                    % look for .xfm file in default location
                    transformPath = sprintf('%s%stransforms%slinear%s',mesh_dir,filesep,filesep,filesep);

                    s = sprintf('%st1_tal.xfm',subj_name);
                    transformFile = fullfile(transformPath,s);
                    fprintf('looking for transformation file (%s)\n', transformFile);

                    if exist(transformFile,'file') ~= 2
                        [tname, tpath, ext] = uigetfile(...
                           {'*.xfm','CIVET TRANSFORM (*.xfm)'},'Select the NIfTI to MNI transform file',mesh_dir);
                        if iequal(tname, 0) || isequal(tpath,0)
                            delete(wbh);
                            return;
                        end
                        transformFile = fullfile(tpath, tname);
                    end   

                    t = importdata(transformFile);
                    transform = [t.data; 0 0 0 1];
                    fprintf('transforming mesh from MNI to original NIfTI coordinates using transformation:\n');

                    mni_to_ras = inv(transform)';
                    RAS_to_MNI = transform';  % for saving


                    % look for whole head cortical thickness
                    % rsl_tlink20 in default location
                    txt_files = dir(fullfile(mesh_dir, 'thickness', '*.txt'));
                    tfile = 'rsl_tlink_20mm.txt';                    
                    if ~isempty(txt_files)
                        for j=1:size(txt_files,1)
                            s = txt_files(j).name;
                            if ~isempty(strfind(s,tfile))
                                thicknessFile = sprintf('%s%s%s%s',mesh_dir,filesep,'thickness',filesep,s);
                                fprintf('reading cortical thickness data from (%s)\n', thicknessFile);
                                fid1 = fopen(thicknessFile,'r');
                                C = textscan(fid1,'%f');
                                thickness = cell2mat(C);
                                fclose(fid1);                           
                            end
                        end
                    end                
                end

                t_vertices = [meshdata.vertices, ones(size(meshdata.vertices,1), 1) ];
                t_vertices = t_vertices * mni_to_ras; 
                t_vertices(:,4) = [];
                meshdata.vertices = t_vertices;

                % scale mesh back to RAS voxels (native) space
                fprintf('rescaling mesh from mm to voxels (scale = %g mm/voxel)\n', tmat.mmPerVoxel);
                meshdata.vertices = meshdata.vertices ./tmat.mmPerVoxel;

                % concatenate left and right meg relative meshes and faces
                if k == 1
                    mesh.numLeftVertices = size(meshdata.vertices,1);
                    mesh.numLeftFaces = size(meshdata.faces,1);
                    vertices = meshdata.vertices;
                    faces = meshdata.faces;
                else
                    mesh.numRightVertices = size(meshdata.vertices,1);
                    mesh.numRightFaces = size(meshdata.faces,1);

                    offset = size(vertices,1);
                    vertices = [vertices; meshdata.vertices];

                    % have to add offset to face indices
                    meshdata.faces = meshdata.faces + offset;
                    faces = [faces; meshdata.faces];
                end        

                clear meshdata;

            end
            
            % should now have mesh vertices (including inflated) in native (RAS) voxel space 
            % based on original .nii volume used to generate meshes

            % *** keep a copy of vertices in MRI (native voxel space for viewing on MRI)  ***
            mesh.mri_vertices = single(vertices);

            % rotate vertices into head coordinates 
            fprintf('rotating and scaling vertices from voxels into MEG (head) coordinates ...\n');

            vertices = [vertices, ones(size(vertices,1), 1) ];
            vertices = vertices * RAS_to_MEG;  % transform to head coordinates and scale to cm
            vertices(:,4) = [];
         
            Nvertices = size(vertices,1);
            Nfaces = size(faces,1);

            % get normal vectors...

            waitbar((h/num_surf), wbh, 'computing surface normals... ',k);

            nvoxels = Nvertices;
            voxels = vertices;

            % for each vertex find its adjoining faces and take the mean of each face
            % normal as the normal for that vertex
            fprintf('computing vertex normals ...\n');

            if exist('bw_computeFaceNormals','file') == 3
                % use mex function       
                normals = bw_computeFaceNormals(double(vertices'), double(faces')); 
                normals = normals';
            else  
                % matlab code to get vertices
                % get normal for each vertex.  This is the mean of the surrounding face normals
                fprintf('mex function bw_computeFaceNormals not found. Computing normals using m-file...\n');
                normals = zeros(nvoxels,3);

                for i=1:size(vertices,1)  
                    voxels(i,1:3) = vertices(i,1:3);
                    faceIdx = i-1;  % since face numbers start at zero
                    [idx, ~] = find(faces == faceIdx);
                    meanVertex = zeros(1,3);
                    numFaces = length(idx);
                    for j=1:numFaces
                        faceNo = idx(j);
                        v1 = vertices(faces(faceNo,1)+1,:);   % add 1 to vertex number since matlab array indices start at 1 not zero.
                        v2 = vertices(faces(faceNo,2)+1,:);
                        v3 = vertices(faces(faceNo,3)+1,:);
                        V = cross((v1-v2),(v2-v3));
                        V = V / norm(V); 
                        meanVertex = meanVertex + V;       
                    end
                    meanVertex = meanVertex / numFaces;
                    meanVertex = meanVertex / norm(meanVertex);  % rescale to unit vector 
                    normals(i,1:3) = meanVertex(1:3);         
                end
            end

           

            for k=80:5:100
                waitbar((h/num_surf), wbh,'finishing surface...');
            end

            %******************************************************************************** 
            % vertices are stored in 3 different coordinate systems
            % 1. mri_vertices:
            %    These are "native" MRI (RAS) voxel coordinates (i.e., original .nii volume).
            % 2. normalized_vertices:
            %    MNI coordinates in mm (native mesh coords for CIVET, have to be computed for Freesurfer).
            % 3. meg_vertices
            %    These are vertices in MEG "head" coordinates (in cm) based on
            %    fiducials used to save the SUBJ_ID.mat file
            %
            %    Face indices are valid for all
            %    *** normals are in MEG relative coordinates !!! 
            % 
            % mesh data saved in single precision to reduce memory / file size etc.
            %********************************************************************************

            mesh.meshType = meshtype;
            mesh.isCIVET = 1;
            mesh.hasFaceNormals = 0;
            mesh.RAS_to_MNI = RAS_to_MNI;
            mesh.MEG_to_RAS = MEG_to_RAS;
          
            % for convenience save direct transformation
            mesh.MNI_to_MEG = inv(RAS_to_MNI) * inv(MEG_to_RAS);
           
            mesh.faces = single(faces);
            
            mesh.thickness = thickness;

            mesh.numVertices = mesh.numLeftVertices + mesh.numRightVertices;
            mesh.numFaces = mesh.numLeftFaces + mesh.numRightFaces;

            if ~isempty(normalized_vertices)
                mesh.normalized_vertices = single(normalized_vertices);  
            else
                mesh.normalized_vertices = [];
            end

            mesh.vertices= single(vertices);        % drawing vertices in MEG coords
            mesh.meg_vertices = single(vertices);   % *new meg_vertices are saved separatelyz
            mesh.normals = single(normals);         % normals in MEG coords
            
            fieldname = surfaces{h, 1};
            mesh.meshName = fieldname;
            
            meshes.(fieldname) = mesh;

        end
        
        if isempty(mesh_path)
            meshFile = fullfile(mri_dir, sprintf('CIVET_SURFACES.mat'));
        else
            meshFile = mesh_path;
        end
        fprintf('Saving mesh data in file %s\n', meshFile);
        save(meshFile ,'-struct', 'meshes'); 
     
        delete (wbh);
        fprintf('all done...\n');

 end

% imports all FSL meshes in chosen folder
 function importFSmeshes(mesh_dir, mri_dir)
        
        wbh = waitbar(0,'Reading meshfiles...');
        
        surfaces{1,1} = 'pial';
        surfaces{2,1} = 'mid';
        surfaces{3,1} = 'white';
        surfaces{4,1} = 'inflated';
        
        % pial surface files
        sf = dir(fullfile(mesh_dir, 'surf', 'lh.pial'));
        if size(sf,1) == 1
            surfaces{1,2} = fullfile(mesh_dir, 'surf', char(sf.name));
        end
        sf = dir(fullfile(mesh_dir, 'surf', 'rh.pial'));
        if size(sf,1) == 1
            surfaces{1,3} = fullfile(mesh_dir, 'surf', char(sf.name));
        end
        
        % white surface files
        sf = dir(fullfile(mesh_dir, 'surf', 'lh.white'));
        if size(sf,1) == 1
            surfaces{3,2} = fullfile(mesh_dir, 'surf', char(sf.name));
        end
        sf = dir(fullfile(mesh_dir, 'surf', 'rh.white'));
        if size(sf,1) == 1
            surfaces{3,3} = fullfile(mesh_dir, 'surf', char(sf.name));
        end
        
        % inflated surface files
        sf = dir(fullfile(mesh_dir, 'surf', 'lh.inflated'));
        if size(sf,1) == 1
            surfaces{4,2} = fullfile(mesh_dir, 'surf', char(sf.name));
        end
        sf = dir(fullfile(mesh_dir, 'surf', 'rh.inflated'));
        if size(sf,1) == 1
            surfaces{4,3} = fullfile(mesh_dir, 'surf', char(sf.name));
        end
        
        % create each surface
        num_surf = size(surfaces,1);
        for h=1:size(surfaces,1)

            s = sprintf('Creating %s surfaces...\n', surfaces{h, 1});
            waitbar((h/num_surf), wbh, s);
            
            mesh =[];
            meshdata = [];
            meshtype = [];
            for k=1:2                
                if h == 1
                    % pial - save copy for mid surface
                    if k == 1 
                        file = surfaces{1,2};               
                        [meshtype, meshdata] = bw_readMeshFile(file);
                        pial_left = meshdata;  
                    else
                        file = surfaces{1,3};                 
                        [meshtype, meshdata] = bw_readMeshFile(file);
                        pial_right = meshdata;
                    end
                end
                if h == 2
                    % create mid-surface from pial and white, assumes they
                    % have same number of vertices (i.e., not downsampled)
                    if k == 1 
                        file = surfaces{3,2};        % read white...       
                        [~, meshdata] = bw_readMeshFile(file);
                        meshdata.vertices = (pial_left.vertices + meshdata.vertices) / 2.0;
                    else
                        file = surfaces{3,3};                 
                        [~, meshdata] = bw_readMeshFile(file);
                        meshdata.vertices = (pial_right.vertices + meshdata.vertices) / 2.0;
                    end
                end
                if h == 3
                   % white
                   if k == 1 
                        file = surfaces{3,2};               
                    else
                        file = surfaces{3,3};                 
                    end
                    [~, meshdata] = bw_readMeshFile(file);                    
                end
                
                if h == 4
                    % inflated
                   if k == 1 
                       file = surfaces{4,2};               
                    else
                        file = surfaces{4,3};                 
                    end
                    [~, meshdata] = bw_readMeshFile(file);     
                end 
                
                if isempty(meshdata)
                    fprintf('failed to read mesh file %s ...\n', file);
                    return;
                end
                
                % ** for freesurfer meshes we have to scale mesh from mm to RAS voxels (native) space
                % ** have to scale then translate origin which is center of image...

                fprintf('rescaling mesh from mm to voxels (scale = %g mm/voxel)\n', tmat.mmPerVoxel);
                meshdata.vertices = meshdata.vertices ./tmat.mmPerVoxel;

                % translate origin to correspond to original RAS volume
                % adding 129 instead of 128 seems to make MNI coords line up on midline better
                % both on .nii and for Talairach coordinates. Also corresponds to
                % conversion shown in surfaceRAS to Talairach conversion documentation
                
                % ****** don't translate if converting fsaverage! *******
                meshdata.vertices = meshdata.vertices + 129;
                
                % meshdata.vertices should now be in native (RAS) voxel space 
                % based on original .nii volume used to generate meshes               
                                           
                % inflated surface appears to be shifted, need to recenter mesh              
                if h == 4                                  
                    mx = max(meshdata.vertices(:,1));
                    mn = min(meshdata.vertices(:,1));
                    if k==1
                        offset = mx - 128;                
                        meshdata.vertices(:,1) = meshdata.vertices(:,1) - offset;
                    else
                        offset = 128 - mn;                
                        meshdata.vertices(:,1) = meshdata.vertices(:,1) + offset;
                    end                                       
                end                  
                
                % save copy of non-downsampled pial vertices in RAS for
                % inflated mesh lookup table
                if h==1
                    if k == 1
                        pial_vertices = meshdata.vertices;
                    else
                        pial_vertices = [pial_vertices; meshdata.vertices];   
                    end
                end
                
                vIndices = [];  
                % for downsampling need to know vertex indices of original mesh
                % can't downsample inflated mesh or it cannot reference
                % back to other meshes in real space
                if downSampleFactor < 1.0 && h ~= 4
                    fprintf('downsampling mesh by factor of %.1f (%.1f %%)\n', downSampleFactor, downSampleFactor * 100.0)
                    v = meshdata.vertices;
                    f = meshdata.faces + 1; 
                    fprintf('Original resolution ( %d vertices, %d faces)\n', size(v,1),size(f,1) );
                    ph = patch('Vertices',v,'Faces',f);
                    [nf, nv] = reducepatch(ph, downSampleFactor);
                    fprintf('Downsampled resolution ( %d vertices, %d faces)\n', size(nv,1),size(nf,1) );
                    meshdata.vertices = nv;
                    meshdata.faces = nf - 1;
                    
                    % get lookup table to indices of old vertices  
                    % find vertices that were preserved by reducePatch
                    [ind,loc] = ismember(v,nv,'rows');
                    
                    % now some tricky code from to map to corresponding indices in original list 
                    locb = loc(ind);                   
                    subind = find(ind);
                    [~,revsor] = sort(locb);
                    vIndices = subind(revsor);
                    
                    % sanity check - should not happen  
                    if size(vIndices) ~= size(nv,1)
                        fprintf('Something went wrong. reducePatch vertices do not match original\n');
                        return;
                    end
                    clear v f nv nf;      
                    
                else
                    vIndices = 1:size(meshdata.vertices,1);
                end                              
                
                % concatenate left and right meshes and faces in mesh struct
                if k == 1
                    mesh.numLeftVertices = size(meshdata.vertices,1);
                    mesh.numLeftFaces = size(meshdata.faces,1);
                    mesh.vertices = meshdata.vertices;
                    mesh.faces = meshdata.faces;
                else
                    mesh.numRightVertices = size(meshdata.vertices,1);
                    mesh.numRightFaces = size(meshdata.faces,1);                   
                    % have to add offset to right hemisphere face indices
                    face_offset = mesh.numLeftVertices;
                    meshdata.faces = meshdata.faces + face_offset;
                    mesh.vertices = [mesh.vertices; meshdata.vertices];
                    mesh.faces = [mesh.faces; meshdata.faces];                
                end
                                
                % get curvature data for this hemisphere
                % - preferably use pial for pial surface if exists
                if k == 1
                    hemiPath = sprintf('%s%ssurf%slh.',mesh_dir, filesep, filesep);
                else
                    hemiPath = sprintf('%s%ssurf%srh.',mesh_dir, filesep, filesep);
                end
                
                % different curv files for pial and white (use white for mid)
                if h == 1
                    curvFile = sprintf('%scurv.pial',hemiPath);
                else
                    curvFile = sprintf('%scurv',hemiPath);
                end
                if exist(curvFile, 'file')
                    fprintf('Reading curvature file %s\n',curvFile);            
                    [meshtype, tdata] = bw_readMeshFile(curvFile);                              
                    tdata.curv = tdata.curv(vIndices);
                    if k==1
                        mesh.curv = tdata.curv;
                    else
                        mesh.curv = [mesh.curv; tdata.curv];
                    end
                else
                    fprintf('curvature file %s not found...\n',curvFile);            
                end
                
                % get thickness for this hemisphere (same for all meshes)          
                thicknessFile = sprintf('%sthickness',hemiPath);
                if exist(thicknessFile, 'file')
                    fprintf('Reading thickness file %s\n',thicknessFile);
                    [meshtype, tdata] = bw_readMeshFile(thicknessFile); 
                    tdata.thickness = tdata.thickness(vIndices);
                    if k==1
                        mesh.thickness = tdata.thickness;
                    else
                        mesh.thickness = [mesh.thickness; tdata.thickness];
                    end
                else
                    fprintf('thickness file %s not found...\n',thicknessFile);            
                end     
                              
            end % for k = 1:2 
            
            clear meshdata            
            
            % now have combined left and right vertices in RAS in
            % mesh.vertices - apply tranformations
           
            % 1. keep a copy of vertices in RAS native voxel space
            mesh.mri_vertices = single(mesh.vertices);
           
            % 2. Compute MNI coordinates, have to get transformation matrix to scale FS meshes
            % from native to MNI (opposite to CIVET above). This can be done
            % using the transformation in talairach.xfm.  However, this
            % transform seems to assume that vertices are still in mm with
            % origin in the center of the image.

            % look for .xfm file in default location
            transformPath = sprintf('%s%smri%stransforms%s',mesh_dir, filesep,filesep,filesep);

            transformFile = sprintf('%stalairach.xfm',transformPath);
            fprintf('looking for transformation file (%s)\n', transformFile);

            if exist(transformFile,'file') ~= 2
                [tname, tpath, ~] = uigetfile(...
                   {'*.xfm','Talairach transform (*.xfm)'},'Select the talairach transform file',mesh_dir);    
                if isequal(tname, 0) || isequal(tpath, 0)
                    delete(wbh);
                    return;
                end
               transformFile = fullfile(tpath, tname);
            end   

            t = importdata(transformFile);
            transform = [t.data; 0 0 0 1];
                  
            if (isFSaverage)
                % for FSaverage meshes they are already in MNI coords so
                % transform is identity matrix with origin [0 0 0] so need to
                % undo the above translation
                transform(:,4) = [-129 -129 -129 1];  
            else
                transform = [t.data; 0 0 0 1];
            end
            
            fprintf('transforming mesh to MNI from original NIfTI coordinates using transformation:\n');
            RAS_to_MNI = transform';        

            % compute correct MNI lookup table for inflated surface
            % if the surface has been downsampled we have to use the
            % vIndices to select correct vertices from the non-downsampled
            % pial surface
            if h == 4
                t_vertices = pial_vertices;
            else
                t_vertices = mesh.vertices;
            end

            % need to scale from voxels to mm for talairach.xfm to work
            fs_vertices = t_vertices * tmat.mmPerVoxel;

            t_vertices = [fs_vertices, ones(size(fs_vertices,1), 1) ];
            t_vertices = t_vertices * RAS_to_MNI; 
            t_vertices(:,4) = [];            
            mesh.normalized_vertices = t_vertices;

            % 3. compute vertices in MEG (head) coordinates RAS_to_MEG
                
            % since inflated surface is not in standard space 
            % swap in pial vertices for MEG coords.  Note we still have to  
            % apply transform to inflated vertices for viewing 
            % pial vertices have already been downsampled 
            fprintf('rotating and scaling vertices from voxels into MEG (head) coordinates ...\n');
            
            t_vertices = [mesh.vertices, ones(size(mesh.vertices,1), 1) ];
            t_vertices = t_vertices * RAS_to_MEG;  % transform to head coordinates and scale to cm
            t_vertices(:,4) = [];
            meg_vertices = t_vertices;
    
            % set drawing vertices to MEG coordinates
            mesh.vertices = single(meg_vertices);
            
            % compute correct MEG lookup table for inflated surface   
            % pial vertices have already been downsampled 
            if h == 4
                t_vertices = pial_vertices;
                t_vertices = [t_vertices, ones(size(t_vertices,1), 1) ];
                t_vertices = t_vertices * RAS_to_MEG;  % transform to head coordinates and scale to cm
                t_vertices(:,4) = [];
                mesh.meg_vertices = single(t_vertices);
            else
                mesh.meg_vertices = single(meg_vertices);
            end
            
            % get normal vectors in MEG coordinates

            waitbar((h/num_surf),wbh, 'computing surface normals... ',k);

            % for each vertex find its adjoining faces and take the mean of each face
            % normal as the normal for that vertex
            fprintf('computing vertex normals ...\n');

            if exist('bw_computeFaceNormals','file') == 3
                % use mex function       
                normals = bw_computeFaceNormals(double(meg_vertices'), double(mesh.faces')); 
                normals = normals';
            else  
                % matlab code to get vertices
                % get normal for each vertex.  This is the mean of the surrounding face normals
                fprintf('mex function bw_computeFaceNormals not found. Computing normals using m-file...\n');
                normals = zeros(nvoxels,3);

                for i=1:size(meg_vertices,1)  
                    voxels(i,1:3) = meg_vertices(i,1:3);
                    faceIdx = i-1;  % since face numbers start at zero
                    [idx, ~] = find(mesh.faces == faceIdx);
                    meanVertex = zeros(1,3);
                    numFaces = length(idx);
                    for j=1:numFaces
                        faceNo = idx(j);
                        v1 = meg_vertices(mesh.faces(faceNo,1)+1,:);   % add 1 to vertex number since matlab array indices start at 1 not zero.
                        v2 = meg_vertices(mesh.faces(faceNo,2)+1,:);
                        v3 = meg_vertices(mesh.faces(faceNo,3)+1,:);
                        V = cross((v1-v2),(v2-v3));
                        V = V / norm(V); 
                        meanVertex = meanVertex + V;       
                    end
                    meanVertex = meanVertex / numFaces;
                    meanVertex = meanVertex / norm(meanVertex);  % rescale to unit vector 
                    normals(i,1:3) = meanVertex(1:3);         
                end
            end
         
            for k=80:5:100
                waitbar((h/num_surf), wbh,'finishing surface...');
            end

            %******************************************************************************** 
            % ** new 
            %  vertices now stored in 4 different coordinate systems
            %  since inflated is now a separate mesh we cannot assume vertices 
            %  (mesh.vertices) are always in real space
            %  drawing vertices (mesh.vertices) now reflect non-inflated and
            %  inflated vertices converted to MEG frame of reference but 
            %  true MEG coordinates are saved in separate array (mesh.meg_vertices).
            %  
            % 1. mri_vertices:
            %    These are "native" MRI (RAS) voxel coordinates (i.e., original .nii volume).
            % 2. normalized_vertices:
            %    MNI coordinates in mm (native mesh coords for CIVET, have to be computed for Freesurfer).
            % 3. vertices
            %    These are image drawing vertices - in MEG "head" frame of reference (in cm) 
            % 4. meg_vertices: 
            %    These are always true MEG coordinates in cm. Same as vertices except for
            %    inflated meshes which are not in true MEG space.
            %
            %    Face indices are same for all
            %    *** normals are in MEG relative coordinates !!! 
            % 
            % mesh data saved in single precision to reduce memory / file size etc.
            %********************************************************************************

            mesh.meshType = meshtype;
            mesh.isCIVET = 0;
            mesh.hasFaceNormals = 0;
            mesh.RAS_to_MNI = RAS_to_MNI;
            mesh.MEG_to_RAS = MEG_to_RAS;            
            % for convenience save direct transformation
            mesh.MNI_to_MEG = inv(RAS_to_MNI) * inv(MEG_to_RAS);
            
            mesh.faces = single(mesh.faces);
            mesh.numVertices = mesh.numLeftVertices + mesh.numRightVertices;
            mesh.numFaces = mesh.numLeftFaces + mesh.numRightFaces;
            mesh.normals = single(normals);                             % normals in MEG coord!                      
            
            fieldname = surfaces{h, 1};
            mesh.meshName = fieldname;         
            meshes.(fieldname) = mesh;   
        end      
        

        if isempty(mesh_path)
            meshFile = fullfile(mri_dir, sprintf('FS_SURFACES.mat'));
        else
            meshFile = mesh_path;
        end       
        
        fprintf('Saving mesh data in file %s\n', meshFile);
        save(meshFile ,'-struct', 'meshes'); 
        
        
        delete (wbh);
        fprintf('all done...\n');
    end


end

