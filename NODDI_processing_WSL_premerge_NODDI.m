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
dataset_directory = '/mnt/c/WSL2_dir/NODDISAH_12/NODDI_processing';
cd(dataset_directory)

%%
NODDI_nii_list = {'DICOM_AX_DTI_NODDI_1_20240119125543_1401',...
    'DICOM_AX_DTI_NODDI_2_20240119125543_1501',...
    'DICOM_AX_DTI_NODDI_3_20240119125543_1601',...
    'DICOM_AX_DTI_NODDI_4_20240119125543_1701'};

calibration = 'DICOM_AX_DTI_Calibration_20240119125543_1301';
t2 = 'DICOM_AX_T2W_CSENSE_20240119125543_601';

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

orig_rotated_bvec = strcat('data_eddy_unwarped', '.eddy_rotated_bvecs');
rename_rotated_bvec = strcat('data_eddy_unwarped', '_eddy_rotated_bvecs.txt');
rename = ['mv ', orig_rotated_bvec, ' ', rename_rotated_bvec];
%%
system(rename)
%%
shortened_base_file_name = erase('data_eddy_unwarped','_eddy_unwarped');
bval_filename = strcat(shortened_base_file_name,'_bval.txt');
bvec_filename = strcat(shortened_base_file_name,'_eddy_unwarped_eddy_rotated_bvecs.txt');

Protocol = FSL2Protocol(bval_filename, bvec_filename);
noddi = MakeModel('WatsonSHStickTortIsoV_B0');
batch_fitting('NODDI_roi.mat', Protocol, noddi, 'FittedParams.mat');

SaveParamsAsNIfTI('FittedParams.mat', 'NODDI_roi.mat', 'b0_bet_mask.nii', 'Case1')

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

