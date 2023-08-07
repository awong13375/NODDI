%% include NODDI toolbox, nifti_matlab, SPM12 in directory
addpath('/usr/local/NODDI_toolbox_v1.05/')
addpath('/usr/local/nifti_matlab/')
addpath('/usr/local/spm12/')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pre-processing + Eddy on NODDI sequences %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set FSL environment
setenv('PATH', [getenv('PATH') ':/usr/local/fsl/bin']);
addpath(genpath('/usr/local/fsl/bin'))

%% go to dataset directory
cd('/mnt/c/WSL2_dir/Patient 2 2023-08-01/DICOM/NODDI_processing')

%%
NODDI_nii_list = {'DICOM_AX_DTI_NODDI_1_20230801120712_901','DICOM_AX_DTI_NODDI_2_20230801120712_1001','DICOM_AX_DTI_NODDI_3_20230801120712_1101','DICOM_AX_DTI_NODDI_4_20230801120712_1201'};
for noddi_files = 1:length(NODDI_nii_list)
    noddi_file = NODDI_nii_list{noddi_files};

    %% Pre-processing
    input = noddi_file;
    output = [noddi_file '_b0'];
    fslroi = ['fslroi' ' ' input ' ' output ' ' '0 ' '1'];
    system(fslroi)

    calibration = 'DICOM_AX_DTI_Calibration_20230801120712_801';
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
    topup = ['topup ' '--imain=b0.nii --datain=acqparams.txt --out=my_output --fout=my_field --iout=my_unwarped_images'];
    system(topup)

    %% visualize top up correction images
    %applytopup = ['applytopup ' '--imain=' output ',' calibration ' --datain=acqparams.txt --inindex=1,2 --topup=my_output --out=my_hifi_images'];
    %system(applytopup)

    %% apply brain extraction tool
    my_unwarped_images_b0=['fslroi ' 'my_unwarped_images ' 'my_unwarped_images_b0 ' '0 ' '1'];
    system(my_unwarped_images_b0)
    bet = ['bet ' 'my_unwarped_images_b0 ' 'nodif_brain_mask'];
    system(bet)
    brain_mask = 'nodif_brain_mask';

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
        ' --topup=my_output --out=' noddi_file '_eddy_unwarped' ' --very_verbose'];
    system(eddy)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NODDI Toolbox analysis %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DCE coregistration on SPM12 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Go to dataset directory %%

cd('/mnt/c/WSL2_dir/Patient 2 2023-08-01/DICOM/DCE_processing/')

%% refrence sequence (T1 10 deg) %%
ref_seq = cellstr('/mnt/c/WSL2_dir/Patient 2 2023-08-01/DICOM/DCE_processing/DICOM_T1map_10_deg_20230801120712_1601.nii')

%% source sequences %%
t1_5_deg = cellstr('/mnt/c/WSL2_dir/Patient 2 2023-08-01/DICOM/DCE_processing/DICOM_T1map_5_deg_20230801120712_1501.nii')
t1_2_deg = cellstr('/mnt/c/WSL2_dir/Patient 2 2023-08-01/DICOM/DCE_processing/DICOM_T1map_2_deg_20230801120712_1401.nii')
dce_seq = cellstr('/mnt/c/WSL2_dir/Patient 2 2023-08-01/DICOM/DCE_processing/DICOM_DCE_5sec_50phases_20230801120712_1801.nii')

%% Extract nii.gz files %%
gunzip(strcat(t1_5_deg, '.gz'))
gunzip(strcat(t1_2_deg, '.gz'))
gunzip(strcat(dce_seq, '.gz'))

%% coregister T1 5 deg sequence %%
nrun = 1; % enter the number of runs here
jobfile = {'/mnt/c/WSL2_dir/batch_coreg_estimate_reslice_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(2, nrun);
for crun = 1:nrun
    inputs{1, crun} = ref_seq; % Coregister: Estimate & Reslice: Reference Image - cfg_files
    inputs{2, crun} = t1_5_deg; % Coregister: Estimate & Reslice: Source Image - cfg_files
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});

%% coregister T1 2 deg sequence %%
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

%% coregister DCE sequences %%

%% specify # of volumes %%
n = 50

for i = 1:n
    dce_seq_vol = cellstr(strcat(dce_seq, ',', string(i)))
    nrun = 1; % enter the number of runs here
    jobfile = {'/mnt/c/WSL2_dir/batch_coreg_estimate_reslice_job.m'};
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



