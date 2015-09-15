function [value, value_gradient] = similarityMetric_PCA( source, target, transform_object, transform_params, varargin)
%SIMILARITYMETRIC_PCA [value, value_gradient] = similarityMetric_PCA( source, target, transform_function, transform_params)
%   Calculated the value of the PCA based similarity measure
%   between the warped source and target (both of type ImageType)
%
%   input arguments:
%       source - moving image to be warped (of type imageType)
%       target - fixed image(of type imageType)
%       transform_object - structure containing the following elements
%           transform_object.transform_function
%           transform_object.transform_gradient_function (optional)
%       transform_params - column vector with the transform parameters
%
% Author: Devis Peressutti, Biomedical Engineering, KCL, 2015



transform_initialization = [];
% Weight of similarity terms
gamma_ = 0.1;
targetmask = [];
targetpad = [];
sourcemask=[];
sourcepad = [];
dbg = false;
reorder = true;


reg_functions = {};
% Parse input arguments
for i=1:size(varargin,2)
% Use initial estimate
    if (strcmp(varargin{i},'transformInit'))
        transform_initialization=varargin{i+1};
        i=i+1;
% Set value of gamma
    elseif (strcmp(varargin{i},'gamma'))
        gamma_=varargin{i+1};
        i=i+1;
% Read in target mask. Only voxels in mask are considered
    elseif (strcmp(varargin{i},'targetMask'))
        targetmask=varargin{i+1};
        i=i+1;
% Padding of target image. Consider voxels with intensities > targetPad
    elseif (strcmp(varargin{i},'targetPad'))
        targetpad=varargin{i+1};
        i=i+1;
% Read in source mask. Only voxels in mask are considered
    elseif (strcmp(varargin{i},'sourceMask'))
        sourcemask=varargin{i+1};
        i=i+1;
% Padding of source image. Consider voxels with intensities > targetPad
    elseif (strcmp(varargin{i},'sourcePad'))
        sourcepad=varargin{i+1};
        i=i+1;
% Debug
    elseif (strcmp(varargin{i},'debug'))
        dbg = true;
        i=i+1;
    end
end


%% calculate metric value

T = double(target.data);
[nx, ny, nz, nt] = size(T);

% Warp source and source mask with current transformation
warped_source = transform_object.transform_function(source, transform_params);
% Source sequence
wsourcemask = ImageType( warped_source );
wsourcemask.data = zeros( warped_source.size' );
wsourcemask.data( warped_source.data(:,:,:,1)>0 ) = 1; 
wsourcemask.data = repmat( wsourcemask.data(:,:,:,1), [1,1,1,nt]);
% wsourcemask = transform_object.transform_function( sourcemask, transform_params);
% Resample warped_source to target
warped_source = resampleImage( warped_source, target, 'interpolation', 'linear');
wsourcemask = resampleImage( wsourcemask, target, 'interpolation', 'NN');
% Get data 
S = double(warped_source.data);

% Compute perimeter of polygon in the PCA space of the target sequence.
% This term constitutes the second term of the similarity measure.
TX = reshape( T(:), [numel(T)/nt, nt] );
[ eVectT, eValT, ZT, meanT, stdDevT, UT, TXn ] = sPCA( TX, eye(nt), nt, 'type', 'dual' );
perimeter = 0;
for n=1:nt-1
    perimeter = perimeter + sqrt( sum( abs( ZT(n+1,:)-ZT(n,:) ).^2) );
end
clear eVectT eValT ZT meanT stdDevT UT TXn TX;

idx_toignore=[];
if numel(targetmask)
    idx_ = find(targetmask.data==0);
    idx_toignore = union(idx_toignore, idx_);
end

if numel(targetpad)
    idx_ = find(target.data<targetpad);
    idx_toignore = union(idx_toignore, idx_);
end

% This ensures that the number of voxels is a multiple of nt
% wsourcemask.data = repmat(wsourcemask.data(:,:,:,1), [1,1,1,nt]);
if numel(wsourcemask)
    idx_ = find(wsourcemask.data==0);
    idx_toignore = union(idx_toignore, idx_);
end

if numel(sourcepad)
    idx_ = find(warped_source.data<sourcepad);
    idx_toignore = union(idx_toignore, idx_);
end

idx_tokeep = setdiff(1:numel(T),idx_toignore);
idx_tokeep = idx_tokeep(:);

% Compute PCA on target sequence considering overlapping volume 
TX = reshape( T(idx_tokeep), [length(idx_tokeep)/nt, nt] );
[ eVectT, eValT, ZT, meanT, stdDevT, UT, TXn ] = sPCA( TX, eye(nt), nt, 'type', 'dual' );

% Compute PCA on source sequence considering overlapping volume (just used 
% to get a demeaned and normalised (if used) matrix)
SX = reshape( S(idx_tokeep), [length(idx_tokeep)/nt, nt] );
[ eVectS, eValS, ZS, meanS, stdDevS, US, SXn ] = sPCA( SX, eye(nt), nt, 'type', 'dual' );

% Project the source sequence on the target subspace
ZST = real( UT'* SXn )';

% Compute perimeter of the projected sequence in target subspace
st_perimeter = 0;
for n=1:nt-1
    st_perimeter = st_perimeter + norm(ZST(n+1,:)-ZST(n,:));
end

% First term of similarity measure
term_1 = sum( sqrt( sum( abs( ZT-ZST ).^2,2 ) ) );
term_2 = perimeter-st_perimeter;

value = gamma_*term_1+(1-gamma_)*term_2;

% disp(value)
% disp(transform_params')

clear eVectT eValT ZT meanT stdDevT UT TXn;
clear eVectS eValS ZS meanS stdDevS US SXn;

%% calculate gradient values (not implemented)
value_gradient = [];

end

