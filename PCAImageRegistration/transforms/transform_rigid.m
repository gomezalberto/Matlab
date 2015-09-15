function [ output_image ] = transform_rigid(  input_image, c, varargin )
%TRANSFORM_RIGID Transforms the input image using the rigid parameters in c
%   This function transforms the input image (ImageType) using the input rigid parameters specified in c. Angles are in radians!
% 
% Author: Alberto Gomez, Biomedical Engineering, KCL, 2013

transform_initialization = [];
dbg = false;
interpolation='linear';
for i=1:size(varargin,2)
    if (strcmp(varargin{i},'init'))
        transform_initialization=varargin{i+1};
        i=i+1;
    elseif (strcmp(varargin{i},'interpolation'))
        interpolation=varargin{i+1};
        i=i+1;
    elseif (strcmp(varargin{i},'debug'))
        dbg = true;
        i=i+1;
    end
    
end


% first apply orientation
n = numel(input_image.spacing);
if n == 4
    n = 3;
end

matrix_use = rigidMatrixFromParameters(c);

imageMatrix = eye(n+1);
imageMatrix(1:n,1:n) = input_image.orientation(1:n,1:n);
imageMatrix(1:n,n+1) = input_image.origin(1:n);
% imageMatrix(1:n,n+1) = -1*((input_image.size(1:n)/2).*input_image.spacing(1:n)+input_image.origin(1:n));

orientation_matrix = matrix_use \ imageMatrix;

newOrigin = matrix_use\[input_image.origin(1:n) ; 1];
no = input_image.origin;
no(1:n)=newOrigin(1:n);

om = input_image.orientation;
om(1:n,1:n) = orientation_matrix(1:n,1:n);

out = ImageType(input_image.size,no,input_image.spacing,om);
out.data = input_image.data;
% then resample to the initial grid

output_image = resampleImage(out, input_image,'interpolation', interpolation);

end



