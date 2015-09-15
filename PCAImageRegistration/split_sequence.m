function [  ] = split_sequence( filename4dsequence, output_dir, varargin )
%split_sequence Split the input 4D sequence into 3D volumes
%   Split the input 4D sequence into 3D volumes. Volumes are written to the
%   output_dir using the prefix provided (default im_).
% 
% Author: Devis Peressutti, Biomedical Engineering, KCL, 2015

prefix = 'im_';
% Parse arguments
for i=1:size(varargin,2)
    if (strcmp(varargin{i},'prefix'))
        prefix = varargin{i+1};
        i = i+1;
    end
end
% Create output folder if not existing
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
% Read in sequence
sequence = read_mhd( filename4dsequence );
n_frames = sequence.size(end);
% Write 3D images in output folder
for n=0:n_frames-1
    write_mhd([output_dir prefix sprintf('%02d',n) '.mhd'], sequence.extractFrame(n+1));
end

end

