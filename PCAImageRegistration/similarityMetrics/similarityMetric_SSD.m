function [value, value_gradient] = similarityMetric_SSD( source, target, transform_object, transform_params, varargin)
%SIMILARITYMETRIC_SSD [value, value_gradient] = similarityMetric_SSD( source, target, transform_function, transform_params)
%   Calculated the value and gradient of the sum of squared differences
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
% Author: Alberto Gomez, Biomedical Engineering, KCL, 2013



transform_initialization = [];
targetmask = [];
targetpad = [];
sourcemask=[];
sourcepad = [];
dbg = false;
reorder = true;


reg_functions = {};

for i=1:size(varargin,2)
    if (strcmp(varargin{i},'transformInit'))
        transform_initialization=varargin{i+1};
        i=i+1;
    elseif (strcmp(varargin{i},'targetMask'))
        targetmask=varargin{i+1};
        i=i+1;
    elseif (strcmp(varargin{i},'targetPad'))
        targetpad=varargin{i+1};
        i=i+1;
    elseif (strcmp(varargin{i},'sourceMask'))
        sourcemask=varargin{i+1};
        i=i+1;
    elseif (strcmp(varargin{i},'sourcePad'))
        sourcepad=varargin{i+1};
        i=i+1;
    elseif (strcmp(varargin{i},'debug'))
        dbg = true;
        i=i+1;
    elseif (strcmp(varargin{i},'add_regularization'))
        n = numel(reg_functions);
        reg_functions{n+1}.reg_fun=varargin{i+1};
        reg_functions{n+1}.lambda=varargin{i+2};
        i=i+2;
    end
end


%% calculate metric value

T = double(target.data);

warped_source = transform_object.transform_function(source, transform_params);

% resample warped_source to target
warped_source = resampleImage(warped_source, target,'interpolation','linear');

S = double(warped_source.data);

idx_toignore=[];
if numel(targetmask)
    idx_ = find(targetmask.data==0);
    idx_toignore = union(idx_toignore, idx_);
end

if numel(targetpad)
    idx_ = find(target.data<targetpad);
    idx_toignore = union(idx_toignore, idx_);
end
% DP missing sourcemask check
if numel(sourcepad)
    idx_ = find(warped_source.data<sourcepad);
    idx_toignore = union(idx_toignore, idx_);
end

idx_tokeep = setdiff(1:numel(T),idx_toignore);
idx_tokeep = idx_tokeep(:);


value =  sum((T(idx_tokeep) - S(idx_tokeep)).^2)/numel(idx_tokeep);

%% calculate gradient values
N = numel(idx_tokeep);
n = warped_source.ndimensions;

% nabla S
out_ = gradientImage(warped_source);

out = reshape(out_,[prod(warped_source.size) n]);
out = out(idx_tokeep,:);

% This is for the reordering
totalrows = repmat(1:N,1,n); % This assumes that the vectorization of out goes starting by x
totalcols = 1:n*N;
totalcols = reshape(totalcols, n,[])';
totalcols = totalcols(:)';
grad_S = sparse(double(totalrows),double(totalcols) ,out(:)' ,N,n*N) ;

% nabla h
x = target.GetPosition(idx_tokeep')';
grad_h = transform_object.transform_gradient_function(x,transform_params);

% total gradient
value_gradient = 2*(S(idx_tokeep)-T(idx_tokeep))'*grad_S * grad_h/numel(idx_tokeep);

for i=1:numel(reg_functions)
    [value_r, grad_r ] = reg_functions{i}.reg_fun(transform_params);
    value = value + reg_functions{i}.lambda*value_r;
    value_gradient = value_gradient + reg_functions{i}.lambda*grad_r;
end


end

