function grad_h = transform_rigid_gradient(  x, params, varargin )
%TRANSFORM_RIGID_GRADIENT Compute gradient of rigid transformation
% 
% Author: Alberto Gomez, Biomedical Engineering, KCL, 2013

transform_initialization = [];
dbg = false;
reorder = true;
for i=1:size(varargin,2)
    if (strcmp(varargin{i},'init'))
        transform_initialization=varargin{i+1};
        i=i+1;
    elseif (strcmp(varargin{i},'no_reorder'))
        reorder = false;
        i=i+1;
    elseif (strcmp(varargin{i},'debug'))
        dbg = true;
        i=i+1;
    end
    
end


NPARAMS = numel(params);

tx = params(1);
ty = params(2);
tz = params(3);
a = params(4);
b = params(5);
c = params(6);

NDIMSO = transform_initialization.output_dimensionality;
noutputdimensions = NDIMSO;
grad_h = zeros(size(x,1)*noutputdimensions,NPARAMS);

%grad_h(:,0*NDIMSO + (1:NDIMSO) ) = g_tx(x, tx,ty,tz,a,b,c);
%grad_h(:,1*NDIMSO + (1:NDIMSO) ) = g_ty(x, tx,ty,tz,a,b,c);
%grad_h(:,2*NDIMSO + (1:NDIMSO) ) = g_tz(x, tx,ty,tz,a,b,c);
%grad_h(:,3*NDIMSO + (1:NDIMSO) ) = g_a(x, tx,ty,tz,a,b,c);
%grad_h(:,4*NDIMSO + (1:NDIMSO) ) = g_b(x, tx,ty,tz,a,b,c);
%grad_h(:,5*NDIMSO + (1:NDIMSO) ) = g_c(x, tx,ty,tz,a,b,c);

grad_h(:,1 ) = g_tx(x, tx,ty,tz,a,b,c);
grad_h(:,2 ) = g_ty(x, tx,ty,tz,a,b,c);
grad_h(:,3 ) = g_tz(x, tx,ty,tz,a,b,c);
grad_h(:,4 ) = g_a(x, tx,ty,tz,a,b,c);
grad_h(:,5) = g_b(x, tx,ty,tz,a,b,c);
grad_h(:,6 ) = g_c(x, tx,ty,tz,a,b,c);

 
% grad_h has as many rows as input data points x, and as many columns as
% noutputdimension x nparams

end


function val = g_tx( x, tx,ty,tz,a,b,c )
    %val = repmat([1 0 0],size(x,1),1);
    val = repmat([1 0 0]',size(x,1),1);
end
%
function val = g_ty( x, tx,ty,tz,a,b,c )
    %val = repmat([0 1 0],size(x,1),1);
    val = repmat([0 1 0]',size(x,1),1);
end
%
function val = g_tz( x, tx,ty,tz,a,b,c )
    %val = repmat([0 0 1],size(x,1),1);
    val = repmat([0 0 1]',size(x,1),1);
end
%
function val = g_a(x, tx,ty,tz,a,b,c)
    val = [ cos(a)*cos(b)*sin(c)*x(:,3)+sin(a)*sin(c)*x(:,2)-cos(a)*sin(b)*sin(c)*x(:,1) ...
        -cos(a)*cos(b)*cos(c)*x(:,3)-sin(a)*cos(c)*x(:,2)+cos(a)*sin(b)*cos(c)*x(:,1)  ...
        -sin(a)*cos(b)*x(:,3)+cos(a)*x(:,2)+sin(a)*sin(b)*x(:,1)]';
    val = val(:);
end
%
function val = g_b(x, tx,ty,tz,a,b,c)
    val = [(cos(b)*cos(c)-sin(a)*sin(b)*sin(c))*x(:,3)+(-sin(a)*cos(b)*sin(c)-sin(b)*cos(c))*x(:,1) ...
        (cos(b)*sin(c)+sin(a)*sin(b)*cos(c))*x(:,3)+(sin(a)*cos(b)*cos(c)-sin(b)*sin(c))*x(:,1) ...
        -cos(a)*sin(b)*x(:,3)-cos(a)*cos(b)*x(:,1)]';
    val = val(:);
end
%
function val = g_c(x, tx,ty,tz,a,b,c)
    val = [(sin(a)*cos(b)*cos(c)-sin(b)*sin(c))*x(:,3)-cos(a)*cos(c)*x(:,2)+(-cos(b)*sin(c)-sin(a)*sin(b)*cos(c))*x(:,1) ...
        (sin(a)*cos(b)*sin(c)+sin(b)*cos(c))*x(:,3)-cos(a)*sin(c)*x(:,2)+(cos(b)*cos(c)-sin(a)*sin(b)*sin(c))*x(:,1) ...
        0*x(:,1)]';
    val = val(:);
end
