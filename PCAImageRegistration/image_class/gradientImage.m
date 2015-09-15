function out  = gradientImage(im,varargin)
% out = gradientImage(im)
%
% Fast gradient
% Works for ImageType
%
% 			 
% Author: Alberto Gomez, Biomedical Engineering, KCL, 2011

difforder = 1;
dbg=false;
for i=1:size(varargin,2)
    if (strcmp(varargin{i},'dbg'))
        dbg=true;
    elseif (strcmp(varargin{i},'order'))
        difforder=varargin{1};
    end
    
end
%----------------------------

out=[];
% Using the derivative of a Gaussian to guarantee smoothness
filter_size = 5;
sigma = 3; % in pixels
for d=1:im.ndimensions
    d_gaussian_kernel = derivativeOfGaussian(filter_size, sigma, d, im.ndimensions);
    dx = convn(im.data,d_gaussian_kernel,'same');
    out = cat(im.ndimensions+1,out,dx);
end

end

function d_gaussian_kernel = derivativeOfGaussian(filter_size, sigma, d, ndimensions)
% example of inputs:
%filter_size = 5;
%sigma = 3; % in pixels
% d=2 (direction along which we do the derivative)

filter_radius = floor(filter_size/2);


% Below will give me the gaussian function along direction d. For example,
% for a 3D image and along direction y:
%d_gaussian_func_2 =@(x,y,z) exp(-(x.^2)/(2*sigma^2)).*-y/(sigma^2).*exp(-(y.^2)/(2*sigma^2)).*exp(-(z.^2)/(2*sigma^2));
str1='';
str2='';

for i=1:ndimensions
    str2=[str2 'x' num2str(i) ','];
    if i==d
        str1 = [str1 '-x' num2str(i) '/(sigma^2).*exp(-(x' num2str(i) '.^2)/(2*sigma^2)).*'];
    else
        str1 = [str1 'exp(-(x' num2str(i) '.^2)/(2*sigma^2)).*'];
    end
end
str3 = repmat('-filter_radius:filter_radius,',1,ndimensions);
eval([ '[' str2(1:end-1) ']=ndgrid(' str3(1:end-1) ');'])
eval(['d_gaussian_kernel=  ' str1(1:end-2) ';'])

end