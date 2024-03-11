%% include NODDI toolbox, nifti_matlab, SPM12 in directory
addpath('/usr/local/NODDI_toolbox_v1.05/')
addpath('/usr/local/nifti_matlab/')
addpath('/usr/local/spm12/')
addpath('/usr/local/ROCKETSHIP/')
addpath('/usr/local/NODDI_toolbox_v1.05/fitting/CreateROI2.m')
addpath('/root/MATLAB Add-Ons/Collections/imtool3D')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pre-processing + Eddy on NODDI sequences %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set FSL environment
setenv('PATH', [getenv('PATH') ':/usr/local/fsl/bin']);
addpath(genpath('/usr/local/fsl/bin'))

%% go to dataset directory
dataset_directory = '/mnt/c/WSL2_dir/NODDISAH_11/NODDI_processing';
cd(dataset_directory)

%%
NODDI_nii_list = {'DICOM_AX_DTI_NODDI_1_20240109112943_801',...
    'DICOM_AX_DTI_NODDI_2_20240109112943_901',...
    'DICOM_AX_DTI_NODDI_3_20240109112943_1001',...
    'DICOM_AX_DTI_NODDI_4_20240109112943_1101'};

calibration = 'DICOM_AX_DTI_Calibration_20240109112943_701';
t2 = 'DICOM_AX_T2W_CSENSE_20240109112943_601';
anat_seq = 'DICOM_Sag_MP-Rage_20240109112943_201';

%% rename bvec and bval files

for noddi_files = 1:length(NODDI_nii_list)
    noddi_file = NODDI_nii_list{noddi_files};
    orig_bval = [noddi_file '.bval'];
    rename_bval = [noddi_file '_bval.txt'];
    orig_bvec = [noddi_file '.bvec'];
    rename_bvec = [noddi_file '_bvec.txt'];

    rename = ['mv ' orig_bval ' ' rename_bval];
    system(rename)

    rename = ['mv ' orig_bvec ' ' rename_bvec];
    system(rename)

end

%% extract b0 images and merge with calibration
i = 1;
for noddi_files = 1:length(NODDI_nii_list)
    noddi_file = NODDI_nii_list{noddi_files};
    
    input = noddi_file;
    output = [noddi_file '_b0_' num2str(i)];

    fslroi = ['fslroi' ' ' input ' ' output ' ' '0 ' '1'];
    system(fslroi)

    if i == 1
        output_1 = output;
    elseif i == 2
        output_2 = output;
    elseif i == 3
        output_3 = output;
    else 
        output_4 = output;
    end
    i = i+1;
end

%% fsl merge b0 scans
fslmerge = ['fslmerge -t b0 ' output_1 ' ' output_2 ' ' output_3 ' ' output_4 ' ' calibration ];

system(fslmerge)

%% generate acqparams.txt
i = 1;
for noddi_files = 1:length(NODDI_nii_list)
    noddi_file = NODDI_nii_list{noddi_files};

    filename = [noddi_file '.json'];
    metadata = jsondecode(fileread(filename));
    
    if i == 1
        total_readout_time_1 = getfield(metadata,"EstimatedTotalReadoutTime");
    elseif i == 2
        total_readout_time_2 = getfield(metadata,"EstimatedTotalReadoutTime");
    elseif i == 3
        total_readout_time_3 = getfield(metadata,"EstimatedTotalReadoutTime");
    else
        total_readout_time_4 = getfield(metadata,"EstimatedTotalReadoutTime");
    end
i = i + 1;
end 

filename = [calibration '.json'];
metadata = jsondecode(fileread(filename));
total_readout_time_cal = getfield(metadata,"EstimatedTotalReadoutTime");

acqparams = fopen('acqparams.txt', 'wt');
txt = ['0 ' '1 ' '0 ' num2str(total_readout_time_1) '\n' '0 ' '1 ' '0 ' num2str(total_readout_time_2)... 
'\n' '0 ' '1 ' '0 ' num2str(total_readout_time_3) '\n' '0 ' '1 ' '0 ' num2str(total_readout_time_4)...
'\n' '0 ' '-1 ' '0 ' num2str(total_readout_time_cal)];
fprintf(acqparams, txt);
fclose(acqparams);

%% top up
topup = ['topup ' '--imain=b0.nii --datain=acqparams.txt --out=my_output --fout=my_field --iout=my_unwarped_images --verbose --nthr=4'];
system(topup)

%% bet 
%% (adjust f and g, higher f (0-1) value is more stringent, higher g (-1-1) means more stringent at top, more liberal at bottom)

bet = ['bet ' 'my_unwarped_images ' 'nodif_brain_mask ' '-A2 ' t2 ' -R -f 0.5 -v'];
%bet = ['bet ' 'my_unwarped_images ' 'nodif_brain_mask ' '-f 0.7'];
system(bet)
brain_mask = 'nodif_brain_mask';

V = niftiread('my_unwarped_images.nii.gz');
mask = niftiread('nodif_brain_mask.nii.gz');
tool = imtool3D(V);
tool.setMask(mask); 

%% merge all 4 NODDI sequences
fslmerge = ['fslmerge -t data ' NODDI_nii_list{1} ' ' NODDI_nii_list{2} ' ' NODDI_nii_list{3} ' ' NODDI_nii_list{4}];
system(fslmerge)

%% merge bvec and bval files
i = 1;
for noddi_files = 1:length(NODDI_nii_list)
    noddi_file = NODDI_nii_list{noddi_files};

    bval_txt = fileread([noddi_file '_bval.txt']);
    bvec_txt = fileread([noddi_file '_bvec.txt']);

    if i == 1  
        bval_file = fopen('data_bval.txt', 'wt');
        fprintf(bval_file, bval_txt);

        bvec_file = fopen('data_bvec.txt', 'wt');
        fprintf(bvec_file, bvec_txt);
    else 
        fprintf(bval_file, bval_txt);
        fprintf(bvec_file, bvec_txt);
    end 
i = i+1;

end 
fclose(bvec_file);
fclose(bval_file);

%% Manual step %%
% When combining all 4 NODDI sequences before running eddy 
% - combined bval file should have no new lines (all values in one line)
% - bvec file should have 4x of each vector (ie. each line should have 4 repetitions of vectors for x, then y, 
% then z dimension), so only 3 lines of values total in the file
% If any of these are not followed, will run into error

%% generate index.txt
nii = niftiread("data.nii.gz");
header = whos ("nii");
full_dimension = getfield(header, "size");
dimension = full_dimension(4);

indx = strings;
for i = 1:dimension
    if i == 1
        indx = strcat(indx, '1');
    else
        indx = strcat(indx,{' '},'1');
    end
end

index = fopen('index.txt', 'wt');
fprintf(index, indx);
fclose(index);

%% run eddy
eddy = ['eddy_cuda10.2 ' '--imain=data' ' --mask=' brain_mask ' --index=index.txt' ...
        ' --acqp=acqparams.txt' ...
        ' --bvecs=data_bvec.txt' ...
        ' --bvals=data_bval.txt' ...
        ' --topup=my_output --out=data_eddy_unwarped' ' --repol --estimate_move_by_susceptibility --very_verbose'];
system(eddy)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DTI analysis %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% rename eddy_rotated bval/bvec
orig_rotated_bvec = strcat('data_eddy_unwarped', '.eddy_rotated_bvecs');
rename_rotated_bvec = strcat('data_eddy_unwarped', '_eddy_rotated_bvecs.txt');
rename = ['mv ', orig_rotated_bvec, ' ', rename_rotated_bvec];
%%
system(rename)
%%
shortened_base_file_name = erase('data_eddy_unwarped','_eddy_unwarped');
bval_filename = strcat(shortened_base_file_name,'_bval.txt');
bvec_filename = strcat(shortened_base_file_name,'_eddy_unwarped_eddy_rotated_bvecs.txt');

%% DTI Fit
dtifit = ['dtifit --data=data_eddy_unwarped.nii' ' --mask=' brain_mask ' --out=dti'...
    ' --bvecs=' bvec_filename ...
    ' --bvals=' bval_filename ' --save_tensor --sse --verbose'];
system(dtifit)
tensors = 'dti_tensor';

%% split dti_tensor to individual components
input = tensors;
tensor_dir = {'Dxx','Dxy','Dxz','Dyy','Dyz','Dzz'};

for i = 1:length(tensor_dir)
    output = [tensors '_' tensor_dir{i}];
    fslroi = ['fslroi' ' ' input ' ' output ' ' num2str(i-1) ' ' num2str(1)];
    system(fslroi)
end

%% coregister to atlas JHU FA atlas
for i = 1:length(tensor_dir)
    ref_seq = ['/mnt/c/WSL2_dir/Atlases/JHU-ICBM-FA-1mm.nii'];
    source_seq = ['dti_tensor_' tensor_dir{i}];
    if i == 1
        flirt_coreg = ['flirt -in ' source_seq ' -ref ' ref_seq ' -out r' source_seq ' -omat invol2refvol.mat -v'];
        system(flirt_coreg)
    else
        flirt_coreg = ['flirt -in ' source_seq ' -ref ' ref_seq ' -out r' source_seq ' -init invol2refvol.mat -applyxfm -v'];
        system(flirt_coreg)
    end
end
%%
gunzip('rdti_tensor.nii.gz')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NODDI Toolbox analysis %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(dataset_directory)

%% NODDI Toolbox analysis %%

b0_roi = ['fslroi ' 'data_eddy_unwarped ' 'NODDI_data_b0 ' '0 1'];
system(b0_roi)
%%
NODDI_bet = ['bet ' 'NODDI_data_b0 ' 'b0_bet_mask' ' -A2 ' t2 ' -R -f 0.5 -v'];
system(NODDI_bet)
b0_bet_mask = 'b0_bet_mask';

gunzip('NODDI_data_b0.nii.gz')
gunzip('b0_bet_mask.nii.gz')
gunzip('data_eddy_unwarped.nii.gz')

%%
CreateROI2('data_eddy_unwarped.nii','b0_bet_mask','NODDI_roi.mat');

Protocol = FSL2Protocol(bval_filename, bvec_filename);
noddi = MakeModel('WatsonSHStickTortIsoV_B0');
batch_fitting('NODDI_roi.mat', Protocol, noddi, 'FittedParams.mat');

SaveParamsAsNIfTI('FittedParams.mat', 'NODDI_roi.mat', 'b0_bet_mask.nii', 'Case1')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Coregister to brain atlas %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Extract brain from MPRAGE anatomical sequence %%
anat_seq = 'DICOM_Sag_MP-Rage_20240109112943_201';

bet = ['bet ' anat_seq ' anat_seq_brain_mask ' '-A2 ' t2 ' -R -f 0.5 -v'];
system(bet)
brain_mask = 'anat_seq_brain_mask';

V = niftiread([anat_seq '.nii.gz']);
mask = niftiread('anat_seq_brain_mask.nii.gz');
tool = imtool3D(V);
tool.setMask(mask); 

gunzip('anat_seq_brain_mask.nii.gz')

%% extract atlases %%

cd('/usr/local/fsl/data/standard')

gunzip('MNI152_T1_1mm_brain.nii.gz')

cd(dataset_directory)

%% Set origin of MPRAGE sequence to anterior commissure with SPM %%
spm

%% coregister MPRAGE and NODDI seqs to MNI atlas using FLIRT %% ***
source_seq = ['anat_seq_brain_mask.nii'];
ref_seq = ['/usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii'];
flirt_coreg = ['flirt -in ' source_seq ' -ref ' ref_seq ' -out r' source_seq ' -omat invol2refvol.mat -v' ];
system(flirt_coreg)
gunzip('ranat_seq_brain_mask.nii.gz')

source_seq = ['Case1_ficvf.nii'];
ref_seq = ['/usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii'];
flirt_coreg = ['flirt -in ' source_seq ' -ref ' ref_seq ' -out r' source_seq ' -omat invol2refvol.mat -v' ];
system(flirt_coreg)
gunzip('rCase1_ficvf.nii.gz')

source_seq = ['Case1_fiso.nii'];
ref_seq = ['/usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii'];
flirt_coreg = ['flirt -in ' source_seq ' -ref ' ref_seq ' -out r' source_seq ' -omat invol2refvol.mat -v' ];
system(flirt_coreg)
gunzip('rCase1_fiso.nii.gz')

source_seq = ['Case1_odi.nii'];
ref_seq = ['/usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii'];
flirt_coreg = ['flirt -in ' source_seq ' -ref ' ref_seq ' -out r' source_seq ' -omat invol2refvol.mat -v' ];
system(flirt_coreg)
gunzip('rCase1_odi.nii.gz')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DCE coregistration on SPM12 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Go to dataset directory %%

cd('/mnt/c/WSL2_dir/NODDISAH_12/DCE_processing')

%% refrence sequence (T1 10 deg) %%
ref_seq = cellstr('DICOM_T1map_10_deg_20240119125543_2001.nii');

%% source sequences %%
t1_5_deg = cellstr('DICOM_T1map_5_deg_20240119125543_1901.nii');
t1_2_deg = cellstr('DICOM_T1map_2_deg_20240119125543_1801.nii');
dce_seq = cellstr('DICOM_DCE_5sec_50phases_20240119125543_2201.nii');

%% Extract nii.gz files %%
gunzip(strcat(ref_seq, '.gz'))
gunzip(strcat(t1_5_deg, '.gz'))
gunzip(strcat(t1_2_deg, '.gz'))
gunzip(strcat(dce_seq, '.gz'))

%% coregister T1 5 deg sequence %%
nrun = 1; % enter the number of runs here
jobfile = {'/mnt/c/Github/NODDI/batch_coreg_estimate_reslice_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(2, nrun);
for crun = 1:nrun
    inputs{1, crun} = ref_seq; % Coregister: Estimate & Reslice: Reference Image - cfg_files
    inputs{2, crun} = t1_5_deg; % Coregister: Estimate & Reslice: Source Image - cfg_files
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});

% coregister T1 2 deg sequence %
nrun = 1; % enter the number of runs here
jobfile = {'/mnt/c/Github/NODDI/batch_coreg_estimate_reslice_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(2, nrun);
for crun = 1:nrun
    inputs{1, crun} = ref_seq; % Coregister: Estimate & Reslice: Reference Image - cfg_files
    inputs{2, crun} = t1_2_deg; % Coregister: Estimate & Reslice: Source Image - cfg_files
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});

% coregister DCE sequences %

% specify # of volumes %
n = 50;

for i = 1:n
    dce_seq_vol = cellstr(strcat(dce_seq, ',', string(i)));
    nrun = 1; % enter the number of runs here
    jobfile = {'/mnt/c/Github/NODDI/batch_coreg_estimate_reslice_job.m'};
    jobs = repmat(jobfile, 1, nrun);
    inputs = cell(2, nrun);
    for crun = 1:nrun
        inputs{1, crun} = ref_seq; % Coregister: Estimate & Reslice: Reference Image - cfg_files
        inputs{2, crun} = dce_seq_vol; % Coregister: Estimate & Reslice: Source Image - cfg_files
    end
    spm('defaults', 'FMRI');
    spm_jobman('run', jobs, inputs{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ROCKETSHIP DCE analysis %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run_parametric

%%

gunzip(strcat('r', erase(dce_seq, '.nii'), '_ROI.nii.gz'));
gunzip(strcat('r', erase(dce_seq, '.nii'), '_AIF.nii.gz'));

%%

dce



%%

