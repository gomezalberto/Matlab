%% Script to transform an image sequence given a rigid transformation

close all
clear
clc

% Set number of degrees of freedom
N_DOF   = 6;
N_DIM   = 3;
spacing = [2 2 2];

%% Set up directories path and file names

input_dir           = '/data/dp11/Data/PCAUS/';     % path to parent directory
data_dir            = [input_dir 'data/vol_B/'];    % path to data
dof_dir             = [input_dir 'dof/'];           % path to dof
target_filename     = [data_dir 'aw1.mhd'];         % target sequence filename
source_filename     = [data_dir 'aw2.mhd'];         % source sequence filename
dof_filename        = [dof_dir 'aw_2_to_1.txt'];    % dof filename
warped_filename     = [data_dir 'aw2to1.mhd'];      % warped sequence filename

%% Read in files

target  = read_mhd( target_filename );
source  = read_mhd( source_filename );
dof     = read_dof_file( dof_filename, N_DOF );

%% Transform source sequence to target sequence

% Centre of rotation is centre of image
CoR = (target.size-1).*target.spacing/2+target.origin;
% Warp source sequence given the input rigid parameters
warped = transform_rigid_ND( source, dof, 'centreOfRotation', CoR(1:N_DIM) );

%% Write warped sequence to file

write_mhd( warped_filename, warped, 'ElementType', 'int16' );

%% Display result on first frame
warped = resampleImage(warped, target, 'interpolation', 'linear');

figure,
subplot(2,3,1); imagesc(target.extractFrame(1).data(:,:,round(end/2))); title('Target')
subplot(2,3,2); imagesc(source.extractFrame(1).data(:,:,round(end/2))); title('Source')
subplot(2,3,3); imagesc(warped.extractFrame(1).data(:,:,round(end/2))); title('Warped source')
subplot(2,3,5); imagesc(target.extractFrame(1).data(:,:,round(end/2))-...
    source.extractFrame(1).data(:,:,round(end/2))); title('Difference')
subplot(2,3,6); imagesc(target.extractFrame(1).data(:,:,round(end/2))-...
    warped.extractFrame(1).data(:,:,round(end/2))); title('Difference')
colormap(gray)
