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
dataset_directory = '/mnt/c/WSL2_dir/VG/DICOM';
cd(dataset_directory)

%%
NODDI_nii_list = {'DICOM_AX_DTI_NODDI_1_20231030112145_1401',...
    'DICOM_AX_DTI_NODDI_2_20231030112145_1501',...
    'DICOM_AX_DTI_NODDI_3_20231030112145_1601',...
    'DICOM_AX_DTI_NODDI_4_20231030112145_1701'};

calibration = 'DICOM_AX_DTI_Calibration_20231030112145_1301';
t2 = 'DICOM_AX_T2W_CSENSE_20231030112145_401';

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

%%

for noddi_files = 1:length(NODDI_nii_list)
    noddi_file = NODDI_nii_list{noddi_files};
    
    %% Pre-processing
    input = noddi_file;
    output = [noddi_file '_b0'];
    fslroi = ['fslroi' ' ' input ' ' output ' ' '0 ' '1'];
    system(fslroi)
    
    fslmerge = ['fslmerge -t b0' ' ' output ' ' calibration ];
    system(fslmerge)

    %% Generate acqparams.txt
    filename = [noddi_file '.json'];
    metadata = jsondecode(fileread(filename));

    %echotime = getfield(metadata, 'EchoTime');
    %grappa = getfield(metadata, 'ParallelReductionFactorInPlane');
    %phase_encode_steps = getfield(metadata,"PhaseEncodingSteps");

    %echo_spacing = (echotime / grappa) / 1000;
    %total_readout_time = echo_spacing * phase_encode_steps;

    total_readout_time = getfield(metadata,"EstimatedTotalReadoutTime");

    acqparams = fopen('acqparams.txt', 'wt');
    txt = ['0 ' '-1 ' '0 ' num2str(total_readout_time) ' ' '\n' '0 ' '1 ' '0 ' num2str(total_readout_time)];
    fprintf(acqparams, txt);
    fclose(acqparams);

    %% Convert bvec bval files to txt
    %bval_files = dir('*.bval');
    %bvec_files = dir('*.bvec');

    %for id = 1:length(bval_files)
        % Get the file name 
    %    [~, f,ext] = fileparts(files(DICOM_AX_DTI_NODDI_2_20230518105839_801.bval).name);

        %[~, f,ext] = fileparts(files(id).name)
    %    rename = strcat(f,'_',ext) ; 
    %    movefile(files(id).name, rename); 
    %end

    %% top up
    topup = ['topup ' '--imain=b0.nii --datain=acqparams.txt --out=my_output --fout=my_field --iout=my_unwarped_images --verbose --nthr=6'];
    system(topup)

    %% visualize top up correction images
    %applytopup = ['applytopup ' '--imain=' output ',' calibration ' --datain=acqparams.txt --inindex=1,2 --topup=my_output --out=my_hifi_images'];
    %system(applytopup)

    %% apply brain extraction tool
    my_unwarped_images_b0=['fslroi ' 'my_unwarped_images ' 'my_unwarped_images_b0 ' '0 ' '1'];
    system(my_unwarped_images_b0)
    bet = ['bet ' 'my_unwarped_images ' 'nodif_brain_mask ' '-A2 ' t2 ' -R -f 0.5 -v'];
    system(bet)
    brain_mask = 'nodif_brain_mask';

    V = niftiread('my_unwarped_images.nii.gz');
    mask = niftiread('nodif_brain_mask.nii.gz');
    tool = imtool3D(V);
    tool.setMask(mask);

    %% DTI Fit
    %dtifit = ['dtifit --data=' input ' --mask=' brain_mask ' --out=dti'...
    %    ' --bvecs=DICOM_AX_DTI_NODDI_1_20230518105839_701_bvec.txt' ...
    %    ' --bvals=DICOM_AX_DTI_NODDI_1_20230518105839_701_bval.txt'];
    %system(dtifit)

    %% generate index.txt file
    nii = niftiread(input);
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

    bvec_file = [noddi_file '_bvec.txt'];
    bval_file = [noddi_file '_bval.txt'];

    %% run eddy
    eddy = ['eddy_cuda10.2 ' '--imain=' input ' --mask=' brain_mask ' --index=index.txt' ...
        ' --acqp=acqparams.txt' ...
        ' --bvecs=' bvec_file ...
        ' --bvals=' bval_file ...
        ' --topup=my_output --out=' noddi_file '_eddy_unwarped' ' --very_verbose --estimate_move_by_susceptibility --repol'];
    system(eddy)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NODDI Toolbox analysis %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% extract nii.gz files %%
cd(dataset_directory)
n = 4;
scan_n = {'901','1001','1101','1201'};
accession_num = '_20231018105619_';

%% NODDI avg 4 sequences %%
for i = 1:n
    scan_num = scan_n(i);
    base_file_name = strcat('DICOM_AX_DTI_NODDI_',string(i),string(accession_num),string(scan_num),'_eddy_unwarped');
    gunzip(strcat(base_file_name, '.nii.gz'))

    if i == 1
        NODDIdata_1 = niftiread(strcat(base_file_name, '.nii'));
    end

    if i == 2
        NODDIdata_2 = niftiread(strcat(base_file_name, '.nii'));
    end

    if i == 3
        NODDIdata_3 = niftiread(strcat(base_file_name, '.nii'));
    end

    if i == 4
        NODDIdata_4 = niftiread(strcat(base_file_name, '.nii'));
       scan_num = scan_n(1);
        i = 1;
        base_file_name = strcat('DICOM_AX_DTI_NODDI_',string(i),string(accession_num) ,string(scan_num),'_eddy_unwarped');
    end
end

xsize = size(NODDIdata_1,1);
ysize = size(NODDIdata_1,2);
zsize = size(NODDIdata_1,3);
ndirs = size(NODDIdata_1,4);

NODDI_data = zeros(xsize, ysize, zsize, ndirs);

for i = 1:xsize
    for j = 1:ysize
        for k = 1:zsize
            for m = 1:ndirs
                NODDI_data(i,j,k,m) = (NODDIdata_1(i,j,k,m) + ...
                    NODDIdata_2(i,j,k,m) + NODDIdata_3(i,j,k,m) + ...
                    NODDIdata_4(i,j,k,m))/4;
            end
        end
    end
end

NODDI_data = single(NODDI_data);

info_NODDI = niftiinfo(strcat(base_file_name, '.nii'));
niftiwrite(NODDI_data,'NODDI_data',info_NODDI)

%% NODDI Toolbox analysis %%

b0_roi = ['fslroi ' 'NODDI_data ' 'NODDI_data_b0 ' '0 1'];
system(b0_roi)

%%
NODDI_bet = ['bet ' 'NODDI_data_b0 ' 'b0_bet_mask ' '-A2 ' t2 ' -R -f 0.7 -v -g -0.01'];
system(NODDI_bet)
b0_bet_mask = 'b0_bet_mask';

V = niftiread('NODDI_data_b0.nii.gz');
    mask = niftiread('b0_bet_mask.nii.gz');
    tool = imtool3D(V);
    tool.setMask(mask);

%%
gunzip('NODDI_data_b0.nii.gz')
gunzip('b0_bet_mask.nii.gz')

orig_rotated_bvec = strcat(base_file_name, '.eddy_rotated_bvecs');
rename_rotated_bvec = strcat(base_file_name, '_eddy_rotated_bvecs.txt');
rename = ['mv ', orig_rotated_bvec, ' ', rename_rotated_bvec];
%%
system(strjoin(rename))
%%
CreateROI2('NODDI_data.nii','b0_bet_mask','NODDI_roi.mat');

shortened_base_file_name = erase(base_file_name,'_eddy_unwarped');
bval_filename = strcat(shortened_base_file_name,'_bval.txt');
bvec_filename = strcat(shortened_base_file_name,'_eddy_unwarped_eddy_rotated_bvecs.txt');

Protocol = FSL2Protocol(bval_filename, bvec_filename);
noddi = MakeModel('WatsonSHStickTortIsoV_B0');
batch_fitting('NODDI_roi.mat', Protocol, noddi, 'FittedParams.mat');

SaveParamsAsNIfTI('FittedParams.mat', 'NODDI_roi.mat', 'b0_bet_mask.nii', 'Case1')

delete('NODDI_data_b0.nii')
delete('b0_bet_mask.nii')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DCE coregistration on SPM12 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Go to dataset directory %%

cd('/mnt/c/WSL2_dir/Patient 6 2023-10-19/NODDISAH6/MRI/DICOM/DCE_processing')

%% refrence sequence (T1 10 deg) %%
ref_seq = cellstr('DICOM_T1map_10_deg_20231018105619_1501.nii');

%% source sequences %%
t1_5_deg = cellstr('DICOM_T1map_5_deg_20231018105619_1401.nii');
t1_2_deg = cellstr('DICOM_T1map_2_deg_20231018105619_1301.nii');
dce_seq = cellstr('DICOM_DCE_5sec_50phases_20231018105619_1601.nii');

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


nrun = 1; % enter the number of runs here %% coregister T1 2 deg sequence
jobfile = {'/mnt/c/Github/NODDI/batch_coreg_estimate_reslice_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(2, nrun);
for crun = 1:nrun
    inputs{1, crun} = ref_seq; % Coregister: Estimate & Reslice: Reference Image - cfg_files
    inputs{2, crun} = t1_2_deg; % Coregister: Estimate & Reslice: Source Image - cfg_files
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});

n = 50; % Coregister DCE sequences, specify # of volumes %%

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

%% extract ROI and AIF %%
gunzip(strcat('r', erase(dce_seq, '.nii'), '_ROI.nii.gz'));
gunzip(strcat('r', erase(dce_seq, '.nii'), '_AIF.nii.gz'));
%%
dce
%%

