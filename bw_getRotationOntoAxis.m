function orientVector = bw_getRotationOntoAxis( axis, angle )
%

% build CTF coordinate system
% origin is midpoint between ears
origin = (left_preauricular_pos + right_preauricular_pos) /2;

% x axis is vector from this origin to Nasion
z_axis= axis/norm(axis);

% y axis is origin to left ear vector
y_axis= left_preauricular_pos - origin;
y_axis=y_axis/norm(y_axis);

% not necessarely perpendicular  (if axis != xaxis)
orthog=cross(z_axis,[1 0 0]);
orthog=orthog/norm(orthog);
x_axis = 
y_axis=cross(z_axis,x_axis);
y_axis=y_axis/norm(y_axis);

% now build 4 x 4 affine transformation matrix

% rotation matrix is constructed from principal axes as unit vectors
% note transpose for correct direction of rotation 
rmat = [ [x_axis 0]; [y_axis 0]; [z_axis 0]; [0 0 0 1] ]';

% scaling matrix from mm to voxels
smat = diag([mmPerVoxel mmPerVoxel mmPerVoxel 1]);

% translation matrix - subtract origin
tmat = diag([1 1 1 1]);
tmat(4,:) = [-origin, 1];

% affine transformation matrix for voxels to CTF is concatenation of these
% three transformations. Order of first two operations is important. Since
% the origin is in units of voxels we must subtract it BEFORE scaling. Also
% since translation vector is in original coords must be also be rotated in
% order to rotate and translate with one matrix operation

M = tmat * smat * rmat;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




