function transform_initialization = transform_rigid_initialize(  input_image, transform_params, varargin )
%TRANSFORM_RIGID_INITIALISE Initialise rigid transformation
% 
% Author: Alberto Gomez, Biomedical Engineering, KCL, 2013

% transform_params = [];

dbg = false;
for i=1:size(varargin,2)
    if (strcmp(varargin{i},'debug'))
        dbg = true;
        i=i+1;
    end
end

NPARAMS = 6;

transform_params.bounds = input_image.GetBounds();
transform_initialization.output_dimensionality = numel(transform_params.bounds)/2;
transform_initialization.params0 = zeros(NPARAMS,1); 

end

