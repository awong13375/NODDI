%% include NODDI toolbox in directory
addpath(genpath('/usr/local/NODDI_tool'))

%% set FSL environment
setenv('PATH', [getenv('PATH') ':/Users/alexw/fsl/bin']);
addpath(genpath('/Users/alexw/fsl/bin/fslroi.sh'))


%% go to dataset directory
cd('/Users/alexw/Downloads/NODDI Project/wetransfer_noddi_2023-06-02_0338/SAH_NODDI/DICOM/')

%% Pre-processing
input = 'DICOM_AX_DTI_NODDI_1_20230518105839_701';
output = 'DTI_NODDI_1_b0';
fslroi = ['fslroi' ' ' input ' ' output ' ' '0 ' '1'];
system(fslroi)

calibration = 'DICOM_AX_DTI_Calibration_20230518105839_1101';
fslmerge = ['fslmerge -t b0' ' ' output ' ' calibration ];
system(fslmerge)

%% Generate acqparams.txt
filename = 'DICOM_AX_DTI_NODDI_1_20230518105839_701.json';
metadata = jsondecode(fileread(filename));

echotime = getfield(metadata, 'EchoTime');
grappa = getfield(metadata, 'ParallelReductionFactorInPlane');
phase_encode_steps = getfield(metadata,"PhaseEncodingSteps");

echo_spacing=(echotime / grappa) / 1000;
total_readout_time=echo_spacing * phase_encode_steps;

acqparams = fopen('acqparams.txt', 'wt');
txt = ['0 ' '-1 ' '0 ' num2str(total_readout_time) ' ' '\n' '0 ' '1 ' '0 ' num2str(total_readout_time)];
fprintf(acqparams, txt)
fclose(acqparams)

%% top up
topup = ['topup ' '--imain=b0.nii --datain=acqparams.txt --out=my_output --fout=my_field --iout=my_unwarped_images'];
system(topup)

%% visualize top up correction images
applytopup = ['applytopup ' '--imain=' output ',' calibration ' --datain=acqparams.txt --inindex=1,2 --topup=my_output --out=my_hifi_images'];
system(applytopup)

%% apply brain extraction tool
bet = ['bet ' output ' ' 'nodif_brain_mask'];
system(bet)

%% eddy
eddy = ['eddy ' '--imain=' input ' ' '--mask=nodif_brain_mask --index=index.txt --acqp=acqparams2.txt --bvecs=bvecs.txt --bvals=bvals.txt --topup=my_output --out=eddy_unwarped'];
system(eddy)

