%% include NODDI toolbox in directory
% addpath(genpath('/usr/local/NODDI_tool'))

%% set FSL environment
setenv('PATH', [getenv('PATH') ':/usr/local/fsl/bin']);
addpath(genpath('/usr/local/fsl/bin'))

%% go to dataset directory
cd('/mnt/c/WSL2_dir/NODDI_post_eddy_2023-07-01/SEBE/DICOM/')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%SAH NODDI file list: {'DICOM_AX_DTI_NODDI_1_20230518105839_701','DICOM_AX_DTI_NODDI_2_20230518105839_801','DICOM_AX_DTI_NODDI_3_20230518105839_901','DICOM_AX_DTI_NODDI_4_20230518105839_1001'}

NODDI_nii_list = {'DICOM_AX_DTI_NODDI_1_20230518105839_701','DICOM_AX_DTI_NODDI_2_20230518105839_801',...
    'DICOM_AX_DTI_NODDI_3_20230518105839_901','DICOM_AX_DTI_NODDI_4_20230518105839_1001'};
for noddi_files = 1:length(NODDI_nii_list)
    noddi_file = NODDI_nii_list{noddi_files};

    %% Pre-processing
    input = noddi_file;
    output = [noddi_file '_b0'];
    fslroi = ['fslroi' ' ' input ' ' output ' ' '0 ' '1'];
    system(fslroi)

    calibration = 'DICOM_AX_DTI_Calibration_20230518105839_1101';
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


