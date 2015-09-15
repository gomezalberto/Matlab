function [ output_image ] = transform_rigid_ND(  input_image, c, varargin )
%TRANSFORM_RIGID_ND Transforms the input ND image using the rigid parameters in c
%   This function transforms the input image (ImageType) using the input rigid parameters specified in c. Angles are in radians!
% 
% Author: Alberto Gomez, Biomedical Engineering, KCL, 2015

transform_initialization = [];
dbg = false;
invert=false;
interpolation='linear';
centreOfRotation = [];
for i=1:size(varargin,2)
    if (strcmp(varargin{i},'init'))
        transform_initialization=varargin{i+1};
        i=i+1;
    elseif (strcmp(varargin{i},'interpolation'))
        interpolation=varargin{i+1};
        i=i+1;
    elseif (strcmp(varargin{i},'centreOfRotation'))
        centreOfRotation=varargin{i+1};
        i=i+1;
    elseif (strcmp(varargin{i},'debug'))
        dbg = true;
        i=i+1;
    elseif (strcmp(varargin{i},'invert'))
        invert = true;
        i=i+1;
    end
    
end

n = numel(input_image.spacing);

matrix_use = rigidMatrixFromParameters(c); % this should be 4x4 ( or 3x3 for 2D which is not currently implemented)

N_SPACE_COORDINATES = min(n,3);

imageMatrix = eye(N_SPACE_COORDINATES+1);
imageMatrix(1:N_SPACE_COORDINATES,1:N_SPACE_COORDINATES) = input_image.orientation(1:N_SPACE_COORDINATES,1:N_SPACE_COORDINATES);
imageMatrix(1:N_SPACE_COORDINATES,N_SPACE_COORDINATES+1) = input_image.origin(1:N_SPACE_COORDINATES);

CORMatrix = eye(N_SPACE_COORDINATES+1);

if numel(centreOfRotation)
    CORMatrix(1:N_SPACE_COORDINATES, N_SPACE_COORDINATES+1) = -centreOfRotation; 
end

% Move the image to be centred about the centreofrotation and then after
% the rotation undo.
matrix_use = CORMatrix \  matrix_use  * CORMatrix; 
if invert
    matrix_use = eye(N_SPACE_COORDINATES+1) / matrix_use;
end
orientation_matrix =   matrix_use \ imageMatrix ; 

newOrigin = matrix_use\[input_image.origin(1:N_SPACE_COORDINATES) ; 1];
no = input_image.origin;
no(1:N_SPACE_COORDINATES)=newOrigin(1:N_SPACE_COORDINATES);

om = input_image.orientation;
om(1:N_SPACE_COORDINATES,1:N_SPACE_COORDINATES) = orientation_matrix(1:N_SPACE_COORDINATES,1:N_SPACE_COORDINATES);

out = ImageType(input_image.size,no,input_image.spacing,om);
out.data = input_image.data;
% then resample to the initial grid

output_image = resampleImage(out, input_image,'interpolation', interpolation);

end



