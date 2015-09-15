%% Script to register two echocardiophy sequences

close all
clear
clc

% Set number of degrees of freedom
N_DOF   = 6;
N_DIM   = 3;
spacing = [2 2 2];

%% Set up directories path and file names

input_dir           = '/media/ag09_local/AGOMEZ_3/other_data/Devis/pca_registration/data/demo_data/';     % path to parent directory
data_dir            = [input_dir 'data/vol_B/'];    % path to data
dof_dir             = [input_dir 'dof/'];           % path to dof
target_filename     = [data_dir 'vol_B_aw_1/vol_B_aw_1.mhd'];         % target sequence filename
source_filename     = [data_dir 'vol_B_aw_2/vol_B_aw_2.mhd'];         % source sequence filename
dof_filename        = [dof_dir 'pca_2_to_1.txt'];    % output dof filename

%% Read in files

target  = read_mhd( target_filename );
source  = read_mhd( source_filename );

%% Smooth and resample target and source sequences

target_sres = resampleImage( target, [], 'spacing', [spacing 1], ...
    'interpolation','linear','blur',[3,3,3,1],[1,1,1,1] );
source_sres = resampleImage( source, [], 'spacing', [spacing 1], ...
    'interpolation','linear','blur',[3,3,3,1],[1,1,1,1] );

%% Resample sequences temporally to have the same number of frames

% target and source number of frames
Nt = target.size(N_DIM+1);
Ns = source.size(N_DIM+1);
% new number of frames
[ N, iN ] = max( [Nt, Ns] );
% Centre of rotation is centre of image
CoR = (target.size-1).*target.spacing/2+target.origin;
% Resample using NearestNeighbourhood interpolation
if (iN==2)
    target_stres = resampleImage( target_sres, [], 'spacing', [target_sres.spacing(1:N_DIM)', Nt/N],...
        'interpolation','NN' );
    target_stres.spacing = [target_sres.spacing(1:N_DIM)', 1];
    source_stres = source_sres;
elseif (iN==1)
    source_stres = resampleImage( source_sres, [], 'spacing', [target_sres.spacing(1:N_DIM)', Ns/N],...
        'interpolation','NN' );
    source_stres.spacing = [target_sres.spacing(1:N_DIM)',1];
    target_stres = target_sres;
end

clear target_sres source_sres

%% Compute masks on target and source sequences

% The mask is created thresholding intensities. The mask at the first frame
% is used for all the other frames
% Target sequence
target_mask = ImageType( target_stres );
target_mask.data = zeros(target_mask.size');
target_mask.data(target_stres.data(:,:,:,1)>0) = 1; 
target_mask.data = repmat(target_mask.data(:,:,:,1), [1,1,1,N]);

%% Retrieve transform model

% Set up transform
transform_params = [];
transform_initialization = transform_rigid_initialize( source_stres, transform_params);
transform_object.transform_function = @(im,params) transform_rigid_ND( ...
    source_stres, params, 'init',transform_initialization, 'centreOfRotation', CoR(1:N_DIM) );
transform_object.transform_gradient_function = [];

%% Retrieve similarity measure

f_similarity = @(params) similarityMetric_PCA(  source_stres, target_stres,...
    transform_object, params, 'transformInit', transform_initialization,...
    'targetMask',target_mask);

%% Retrieve optimizer
options.DisplayIter = true;
options.th_Afvalue = 1E-06;
options.th_Aparams = 5E-03;
options.initialStepNorm = 10;
options.lambda_decrease_rate = 0.75;
options.max_iterations = 20;
options.use_gradient = true;
options.Scaling = [1 1 1 2*pi/450 2*pi/450 2*pi/450]';
f_optimizer = @optimizer_gradientDescent_normalized;

%% Run optimization

[ c, output ] = f_optimizer(f_similarity, transform_initialization.params0, options);

%% Write parameters to filename (compliant to IRTK)

write_dof_file( dof_filename, c );

%% Display result on first frame
warped_source = transform_object.transform_function(source,c);
warped_source = resampleImage(warped_source, target, 'interpolation', 'linear');

figure,
subplot(2,3,1); imagesc(target.extractFrame(1).data(:,:,round(end/2))); title('Target')
subplot(2,3,2); imagesc(source.extractFrame(1).data(:,:,round(end/2))); title('Source')
subplot(2,3,3); imagesc(warped_source.extractFrame(1).data(:,:,round(end/2))); title('Warped source')
subplot(2,3,5); imagesc(target.extractFrame(1).data(:,:,round(end/2))-source.extractFrame(1).data(:,:,round(end/2))); title('Difference')
subplot(2,3,6); imagesc(target.extractFrame(1).data(:,:,round(end/2))-warped_source.extractFrame(1).data(:,:,round(end/2))); title('Difference')
colormap(gray)
