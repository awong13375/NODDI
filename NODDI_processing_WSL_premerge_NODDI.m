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

%% Coreg FA map
fa_map = 'dti_FA';
ref_seq = ['/mnt/c/WSL2_dir/Atlases/JHU-ICBM-FA-1mm.nii'];

flirt_coreg = ['flirt -in ' fa_map ' -ref ' ref_seq ' -out r' fa_map ' -omat invol2refvol.mat -v'];
system(flirt_coreg)

%% Coreg diffusivity maps based on FA registration matrix
for i = 1:length(tensor_dir)
    source_seq = ['dti_tensor_' tensor_dir{i}];

    flirt_coreg = ['flirt -in ' source_seq ' -ref ' ref_seq ' -out r' source_seq ' -init invol2refvol.mat -applyxfm -v'];
    system(flirt_coreg)
end

%% Create 5mm spherical ROI and extract values (SLF=association, SCR=projection)

% Single voxel ROI
roi_L_SLF = '128 1 110 1 99 1 0 1';
roi_L_SCR = '116 1 110 1 99 1 0 1';
roi_R_SLF = '51 1 110 1 99 1 0 1';
roi_R_SCR = '64 1 110 1 99 1 0 1';

Dxx = 'rdti_tensor_Dxx.nii.gz';
Dyy = 'rdti_tensor_Dyy.nii.gz';
Dzz = 'rdti_tensor_Dzz.nii.gz';

% Left Dxx projection
fsl_maths = ['fslmaths ' ref_seq ' -mul 0 -add 1 -roi ' roi_L_SCR ' ROI.nii.gz -odt float'];
system(fsl_maths)

fsl_maths = ['fslmaths ROI.nii.gz -kernel sphere 5 -fmean sphere.nii.gz -odt float'];
system(fsl_maths)

fsl_maths = ['fslmaths sphere.nii.gz -bin sphere_bin.nii.gz'];
system(fsl_maths)

fsl_meants = ['fslmeants -i ' Dxx ' -m sphere_bin.nii.gz -o out.txt'];
system(fsl_meants)
Dxx_assoc = importdata('out.txt');

% Left Dxx association
fsl_maths = ['fslmaths ' ref_seq ' -mul 0 -add 1 -roi ' roi_L_SLF ' ROI.nii.gz -odt float'];
system(fsl_maths)

fsl_maths = ['fslmaths ROI.nii.gz -kernel sphere 5 -fmean sphere.nii.gz -odt float'];
system(fsl_maths)

fsl_maths = ['fslmaths sphere.nii.gz -bin sphere_bin.nii.gz'];
system(fsl_maths)

fsl_meants = ['fslmeants -i ' Dxx ' -m sphere_bin.nii.gz -o out.txt'];
system(fsl_meants)
Dxx_proj = importdata('out.txt');

% Left Dyy projection
fsl_maths = ['fslmaths ' ref_seq ' -mul 0 -add 1 -roi ' roi_L_SCR ' ROI.nii.gz -odt float'];
system(fsl_maths)

fsl_maths = ['fslmaths ROI.nii.gz -kernel sphere 5 -fmean sphere.nii.gz -odt float'];
system(fsl_maths)

fsl_maths = ['fslmaths sphere.nii.gz -bin sphere_bin.nii.gz'];
system(fsl_maths)

fsl_meants = ['fslmeants -i ' Dyy ' -m sphere_bin.nii.gz -o out.txt'];
system(fsl_meants)
Dyy_proj = importdata('out.txt');

% Left Dzz association
fsl_maths = ['fslmaths ' ref_seq ' -mul 0 -add 1 -roi ' roi_L_SLF ' ROI.nii.gz -odt float'];
system(fsl_maths)

fsl_maths = ['fslmaths ROI.nii.gz -kernel sphere 5 -fmean sphere.nii.gz -odt float'];
system(fsl_maths)

fsl_maths = ['fslmaths sphere.nii.gz -bin sphere_bin.nii.gz'];
system(fsl_maths)

fsl_meants = ['fslmeants -i ' Dzz ' -m sphere_bin.nii.gz -o out.txt'];
system(fsl_meants)
Dzz_assoc = importdata('out.txt');

DTI_ALPS_left = ((Dxx_proj+Dxx_assoc)/2)/((Dyy_proj+Dzz_assoc)/2);


% Right Dxx projection
fsl_maths = ['fslmaths ' ref_seq ' -mul 0 -add 1 -roi ' roi_R_SCR ' ROI.nii.gz -odt float'];
system(fsl_maths)

fsl_maths = ['fslmaths ROI.nii.gz -kernel sphere 5 -fmean sphere.nii.gz -odt float'];
system(fsl_maths)

fsl_maths = ['fslmaths sphere.nii.gz -bin sphere_bin.nii.gz'];
system(fsl_maths)

fsl_meants = ['fslmeants -i ' Dxx ' -m sphere_bin.nii.gz -o out.txt'];
system(fsl_meants)
Dxx_assoc = importdata('out.txt');

% Right Dxx association
fsl_maths = ['fslmaths ' ref_seq ' -mul 0 -add 1 -roi ' roi_R_SLF ' ROI.nii.gz -odt float'];
system(fsl_maths)

fsl_maths = ['fslmaths ROI.nii.gz -kernel sphere 5 -fmean sphere.nii.gz -odt float'];
system(fsl_maths)

fsl_maths = ['fslmaths sphere.nii.gz -bin sphere_bin.nii.gz'];
system(fsl_maths)

fsl_meants = ['fslmeants -i ' Dxx ' -m sphere_bin.nii.gz -o out.txt'];
system(fsl_meants)
Dxx_proj = importdata('out.txt');

% Right Dyy projection
fsl_maths = ['fslmaths ' ref_seq ' -mul 0 -add 1 -roi ' roi_R_SCR ' ROI.nii.gz -odt float'];
system(fsl_maths)

fsl_maths = ['fslmaths ROI.nii.gz -kernel sphere 5 -fmean sphere.nii.gz -odt float'];
system(fsl_maths)

fsl_maths = ['fslmaths sphere.nii.gz -bin sphere_bin.nii.gz'];
system(fsl_maths)

fsl_meants = ['fslmeants -i ' Dyy ' -m sphere_bin.nii.gz -o out.txt'];
system(fsl_meants)
Dyy_proj = importdata('out.txt');

% Right Dzz association
fsl_maths = ['fslmaths ' ref_seq ' -mul 0 -add 1 -roi ' roi_R_SLF ' ROI.nii.gz -odt float'];
system(fsl_maths)

fsl_maths = ['fslmaths ROI.nii.gz -kernel sphere 5 -fmean sphere.nii.gz -odt float'];
system(fsl_maths)

fsl_maths = ['fslmaths sphere.nii.gz -bin sphere_bin.nii.gz'];
system(fsl_maths)

fsl_meants = ['fslmeants -i ' Dzz ' -m sphere_bin.nii.gz -o out.txt'];
system(fsl_meants)
Dzz_assoc = importdata('out.txt');

DTI_ALPS_right = ((Dxx_proj+Dxx_assoc)/2)/((Dyy_proj+Dzz_assoc)/2);

Indices = ["DTI_ALPS_left";"DTI_ALPS_right"];
Values = [DTI_ALPS_right; DTI_ALPS_left];

T = table(Indices, Values);
writetable(T, 'DTI_ALPS.csv')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NODDI Toolbox analysis %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(dataset_directory)

%% NODDI Toolbox analysis %%

b0_roi = ['fslroi ' 'data_eddy_unwarped ' 'NODDI_data_b0 ' '0 1'];
system(b0_roi)

%% (adjust f and g, higher f (0-1) value is more stringent, higher g (-1-1) means more stringent at top, more liberal at bottom)
NODDI_bet = ['bet ' 'NODDI_data_b0 ' 'b0_bet_mask' ' -A2 ' t2 ' -R -f 0.5 -v'];
system(NODDI_bet)

V = niftiread(['NODDI_data_b0.nii.gz']);
mask = niftiread('b0_bet_mask.nii.gz');
tool = imtool3D(V);
tool.setMask(mask);

b0_bet_mask = 'b0_bet_mask';
%%
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

%% Go to 3D slicer and crop out neck, save as no_neck

%% Extract brain from MPRAGE anatomical sequence %%
anat_seq = 'no_neck';

bet = ['bet ' anat_seq ' anat_seq_brain_mask ' '-A2 ' t2 ' -R -f 0.5 -g -0.22 -v'];
system(bet)
brain_mask = 'anat_seq_brain_mask';

V = niftiread([anat_seq '.nii.gz']);
mask = niftiread('anat_seq_brain_mask.nii.gz');
tool = imtool3D(V);
tool.setMask(mask); 

%% extract atlases 
% cd('/usr/local/fsl/data/standard')
% gunzip('MNI152_T1_1mm_brain.nii.gz')
% cd(dataset_directory)

%% MNI Analysis %%

% coregister T1 MPRAGE to NODDI
source_seq = ['anat_seq_brain_mask.nii'];
ref_seq = ['Case1_ficvf.nii'];
flirt_coreg = ['flirt -in ' source_seq ' -ref ' ref_seq ' -out r_noddi_' source_seq ' -omat invol2refvol_T1_NODDI.mat -v' ];
system(flirt_coreg)

% coregister NODDI-registered T1 MPRAGE to atlas
source_seq = ['r_noddi_anat_seq_brain_mask.nii'];
ref_seq = ['/mnt/c/WSL2_dir/Atlases/MNI-maxprob-thr0-1mm.nii'];
flirt_coreg = ['flirt -in ' source_seq ' -ref ' ref_seq ' -out r_MNI_' source_seq ' -omat invol2refvol_T1_NODDI_MNI.mat -v' ];
system(flirt_coreg)

% apply transformation matrix to NODDI 
source_seq = ['Case1_ficvf.nii'];
ref_seq = ['/mnt/c/WSL2_dir/Atlases/MNI-maxprob-thr0-1mm.nii'];
flirt_coreg = ['flirt -in ' source_seq ' -ref ' ref_seq ' -out r_MNI_' source_seq ' -init invol2refvol.mat -applyxfm -v' ];
system(flirt_coreg)

source_seq = ['Case1_fiso.nii'];
ref_seq = ['/mnt/c/WSL2_dir/Atlases/MNI-maxprob-thr0-1mm.nii'];
flirt_coreg = ['flirt -in ' source_seq ' -ref ' ref_seq ' -out r_MNI_' source_seq ' -init invol2refvol.mat -applyxfm -v' ];
system(flirt_coreg)

source_seq = ['Case1_odi.nii'];
ref_seq = ['/mnt/c/WSL2_dir/Atlases/MNI-maxprob-thr0-1mm.nii'];
flirt_coreg = ['flirt -in ' source_seq ' -ref ' ref_seq ' -out r_MNI_' source_seq ' -init invol2refvol.mat -applyxfm -v' ];
system(flirt_coreg)

%% Create atlas masks on FSL %%

%% Binarize masks
mask_list = {'Caudate','Cerebellum','Frontal Lobe','Insula','Occipital Lobe','Parietal Lobe',...
    'Putamen','Temporal Lobe','Thalamus'};

for n = 1:length(mask_list)
    mask = mask_list{n};
    instring = ['mni_prob_' mask '.nii.gz'];
    outstring = regexprep(instring, ' ', '_');
    outstring = regexprep(outstring, '(','');
    outstring = regexprep(outstring, ')','');
    outstring = regexprep(outstring, ',','');

    if isfile(instring) == 0
        continue
    elseif instring == outstring
        continue
    else
        movefile(instring, outstring)
    end
end

mask_list = regexprep(mask_list, ' ', '_');
mask_list = regexprep(mask_list, '(','');
mask_list = regexprep(mask_list, ')','');
mask_list = regexprep(mask_list, ',','');

for n = 1:length(mask_list)
    mask = mask_list{n};
    fsl_maths = ['fslmaths mni_prob_' mask '.nii.gz -bin mni_prob_' mask '_bin.nii.gz'];
    system(fsl_maths)
end

%% Extract data from NODDI sequences fiso = FWF, ficvf = NDI
source_seq_list = {'r_MNI_Case1_ficvf','r_MNI_Case1_fiso','r_MNI_Case1_odi'};

for n = 1:length(source_seq_list)
    source_seq = source_seq_list{n};
    if n ==1 
        for m = 1:length(mask_list)
            %generate voxel number, mean, and SD for each brain region
            mask = mask_list{m};
            fslstats = ['fslstats ' source_seq ' -k mni_prob_' mask '_bin.nii.gz -V -M -S >data.txt'];
            system(fslstats)
            %Tabulate data
            if m == 1
                intensities = num2cell(importdata('data.txt'));
                Tndi = cell2table (intensities, 'VariableNames', {'Voxels_ndi','Volumes_ndi','NDI Mean','SD_ndi'});

            else 
                intensities = num2cell(importdata('data.txt'));
                T2 = cell2table (intensities, 'VariableNames', {'Voxels_ndi','Volumes_ndi','NDI Mean','SD_ndi'});
                Tndi = [Tndi;T2];
            end
            
            %Extract raw intensity values per voxel for each brain region
            fsl_meants = ['fslmeants -i ' source_seq ' -m mni_prob_' mask '_bin.nii.gz -o out.txt --showall'];
            system(fsl_meants)
            system('tail -1 out.txt >out1.txt')
            if m == 1
                raw_values_1 = importdata('out1.txt');
            elseif m == 2 
                raw_values_2 = importdata('out1.txt');
            elseif m == 3 
                raw_values_3 = importdata('out1.txt');
            elseif m == 4 
                raw_values_4 = importdata('out1.txt');
            elseif m == 5 
                raw_values_5 = importdata('out1.txt');
            elseif m == 6 
                raw_values_6 = importdata('out1.txt');
            elseif m == 7 
                raw_values_7 = importdata('out1.txt');
            elseif m == 8 
                raw_values_8 = importdata('out1.txt');
            elseif m == 9 
                raw_values_9 = importdata('out1.txt');
            end
        end

        Trv_ndi = num2cell(padcat(raw_values_1,raw_values_2,raw_values_3,raw_values_4, ...
                  raw_values_5,raw_values_6,raw_values_7,raw_values_8,raw_values_9));
        Trv_ndi = cell2table(Trv_ndi',"VariableNames", mask_list);
    
    end

    if n == 2
        for m = 1:length(mask_list)
            mask = mask_list{m};
            fslstats = ['fslstats ' source_seq ' -k mni_prob_' mask '_bin.nii.gz -V -M -S >data.txt'];
            system(fslstats)
            if m == 1
                intensities = num2cell(importdata('data.txt'));
                Tfwf = cell2table (intensities, 'VariableNames', {'Voxels_fwf','Volumes_fwf','FWF Mean','SD_fwf'});
            else 
                intensities = num2cell(importdata('data.txt'));
                T2 = cell2table (intensities, 'VariableNames', {'Voxels_fwf','Volumes_fwf','FWF Mean','SD_fwf'});
                Tfwf = [Tfwf;T2];
            end

            fsl_meants = ['fslmeants -i ' source_seq ' -m mni_prob_' mask '_bin.nii.gz -o out.txt --showall'];
            system(fsl_meants)
            system('tail -1 out.txt >out1.txt')
            if m == 1
                raw_values_1 = importdata('out1.txt');
            elseif m == 2 
                raw_values_2 = importdata('out1.txt');
            elseif m == 3 
                raw_values_3 = importdata('out1.txt');
            elseif m == 4 
                raw_values_4 = importdata('out1.txt');
            elseif m == 5 
                raw_values_5 = importdata('out1.txt');
            elseif m == 6 
                raw_values_6 = importdata('out1.txt');
            elseif m == 7 
                raw_values_7 = importdata('out1.txt');
            elseif m == 8 
                raw_values_8 = importdata('out1.txt');
            elseif m == 9 
                raw_values_9 = importdata('out1.txt');
            end
        end
        Trv_fwf = num2cell(padcat(raw_values_1,raw_values_2,raw_values_3,raw_values_4, ...
                  raw_values_5,raw_values_6,raw_values_7,raw_values_8,raw_values_9));
        Trv_fwf = cell2table(Trv_fwf',"VariableNames",mask_list);
    end

    if n == 3
        for m = 1:length(mask_list)
            mask = mask_list{m};
            fslstats = ['fslstats ' source_seq ' -k mni_prob_' mask '_bin.nii.gz -V -M -S >data.txt'];
            system(fslstats)

            if m == 1
                intensities = num2cell(importdata('data.txt'));
                Todi = cell2table (intensities, 'VariableNames', {'Voxels_odi','Volumes_odi','ODI Mean','SD_odi'});
            else 
                intensities = num2cell(importdata('data.txt'));
                T2 = cell2table (intensities, 'VariableNames', {'Voxels_odi','Volumes_odi','ODI Mean','SD_odi'});
                Todi = [Todi;T2];
            end

            fsl_meants = ['fslmeants -i ' source_seq ' -m mni_prob_' mask '_bin.nii.gz -o out.txt --showall'];
            system(fsl_meants)
            system('tail -1 out.txt >out1.txt')
            if m == 1
                raw_values_1 = importdata('out1.txt');
            elseif m == 2 
                raw_values_2 = importdata('out1.txt');
            elseif m == 3 
                raw_values_3 = importdata('out1.txt');
            elseif m == 4 
                raw_values_4 = importdata('out1.txt');
            elseif m == 5 
                raw_values_5 = importdata('out1.txt');
            elseif m == 6 
                raw_values_6 = importdata('out1.txt');
            elseif m == 7 
                raw_values_7 = importdata('out1.txt');
            elseif m == 8 
                raw_values_8 = importdata('out1.txt');
            elseif m == 9 
                raw_values_9 = importdata('out1.txt');
            end
        end
        Trv_odi = num2cell(padcat(raw_values_1,raw_values_2,raw_values_3,raw_values_4, ...
                  raw_values_5,raw_values_6,raw_values_7,raw_values_8,raw_values_9));
        Trv_odi = cell2table(Trv_odi',"VariableNames",mask_list);

    end
end

mask_list_key = mask_list';

key = cell2table(mask_list_key, 'VariableNames',{'Region'});

Tndi = [key, Tndi];
Todi = [key, Todi];
Tfwf = [key, Tfwf];

Tjoin = join(Tndi, Todi);
T_combined = join(Tjoin, Tfwf);
%%
writetable(T_combined, 'MNI_NODDI_indices.csv')
%%
writetable(Trv_ndi, 'MNI_NDI_raw_intensities.csv')
writetable(Trv_fwf, 'MNI_FWF_raw_intensities.csv')
writetable(Trv_odi, 'MNI_ODI_raw_intensities.csv')

%% Segmentation for Havard Oxford and John Hopkins analysis

% FAST segmentation (0 = CSF, 1 = GM, 2 = WM)
source_seq = ['anat_seq_brain_mask.nii.gz'];
fast_seg = ['fast -B -S 1 -n 3 -t 1 -v ' source_seq];
system(fast_seg)

%% Harvard Oxford Analysis %%

%% Cortical/Grey Matter Analysis

% Binarize and apply transformation matrix from MNI transformations to GM mask
source_seq = ['anat_seq_brain_mask_pve_1'];
fsl_maths = ['fslmaths ' source_seq ' -bin ' source_seq '_bin.nii.gz'];
system(fsl_maths)

source_seq = ['anat_seq_brain_mask_pve_1_bin'];
ref_seq = ['Case1_ficvf'];

flirt_coreg = ['flirt -in ' source_seq ' -ref ' ref_seq ' -out r_noddi_' source_seq ' -init invol2refvol_T1_NODDI.mat -applyxfm -v' ];
system(flirt_coreg)

source_seq = ['r_noddi_anat_seq_brain_mask_pve_1_bin'];
ref_seq = ['/mnt/c/WSL2_dir/Atlases/MNI-maxprob-thr0-1mm.nii'];

flirt_coreg = ['flirt -in ' source_seq ' -ref ' ref_seq ' -out r_MNI_' source_seq ' -init invol2refvol_T1_NODDI_MNI.mat -applyxfm -v' ];
system(flirt_coreg)


%% Rebinarize mask and manually adjust threshold (higher is more stringent)
GM_mask = ['r_MNI_r_noddi_anat_seq_brain_mask_pve_1_bin'];
fsl_maths = ['fslmaths ' GM_mask ' -thr 0.55 -bin ' GM_mask '_bin.nii.gz'];
system(fsl_maths)

%% Extract GM from NODDI using GM mask
GM_mask = ['r_MNI_r_noddi_anat_seq_brain_mask_pve_1_bin_bin'];
ref_seq = ['r_MNI_Case1_ficvf'];
fslmaths = ['fslmaths ' ref_seq ' -mul ' GM_mask ' ' ref_seq '_GM'];
system(fslmaths)

ref_seq = ['r_MNI_Case1_fiso'];
fslmaths = ['fslmaths ' ref_seq ' -mul ' GM_mask ' ' ref_seq '_GM'];
system(fslmaths)

ref_seq = ['r_MNI_Case1_odi'];
fslmaths = ['fslmaths ' ref_seq ' -mul ' GM_mask ' ' ref_seq '_GM'];
system(fslmaths)


%%

%% Create atlas masks on FSL %%

%% Binarize masks

%% WILL BREAK CODE -> Manually edit Heschel's gyrus and remove apostrophe %%

% Remove symbols so code won't break
mask_list = {'Frontal Pole','Parahippocampal Gyrus, posterior division','Parahippocampal Gyrus, anterior division','Frontal Orbital Cortex',...
    'Cuneal Cortex','Precuneous Cortex','Cingulate Gyrus, posterior division','Cingulate Gyrus, anterior division','Paracingulate Gyrus',...
    'Subcallosal Cortex','Juxtapositional Lobule Cortex (formerly Supplementary Motor Cortex)','Frontal Medial Cortex','Intracalcarine Cortex',...
    'Lateral Occipital Cortex, inferior division','Lateral Occipital Cortex, superior division','Angular Gyrus','Supramarginal Gyrus, posterior division',...
    'Supramarginal Gyrus, anterior division','Superior Parietal Lobule','Postcentral Gyrus','Inferior Temporal Gyrus, temporooccipital part',...
    'Inferior Temporal Gyrus, posterior division','Inferior Temporal Gyrus, anterior division','Middle Temporal Gyrus, temporooccipital part',...
    'Middle Temporal Gyrus, posterior division','Middle Temporal Gyrus, anterior division','Superior Temporal Gyrus, posterior division',...
    'Superior Temporal Gyrus, anterior division','Temporal Pole','Precentral Gyrus','Inferior Frontal Gyrus, pars opercularis',...
    'Inferior Frontal Gyrus, pars triangularis','Middle Frontal Gyrus','Superior Frontal Gyrus','Insular Cortex','Occipital Pole',...
    'Supracalcarine Cortex','Planum Temporale','Heschls Gyrus (includes H1 and H2)','Planum Polare','Parietal Operculum Cortex',...
    'Central Opercular Cortex','Frontal Operculum Cortex','Occipital Fusiform Gyrus','Temporal Occipital Fusiform Cortex',...
    'Temporal Fusiform Cortex, posterior division','Temporal Fusiform Cortex, anterior division','Lingual Gyrus'};

for n = 1:length(mask_list)
    mask = mask_list{n};
    instring = ['harvardoxford-cortical_prob_' mask '.nii.gz'];
    outstring = regexprep(instring, ' ', '_');
    outstring = regexprep(outstring, '(','');
    outstring = regexprep(outstring, ')','');
    outstring = regexprep(outstring, ',','');

    if isfile(instring) == 0
        continue
    else
        movefile(instring, outstring)
    end
end

mask_list = regexprep(mask_list, ' ', '_');
mask_list = regexprep(mask_list, '(','');
mask_list = regexprep(mask_list, ')','');
mask_list = regexprep(mask_list, ',','');

% binarize masks
for n = 1:length(mask_list)
    mask = mask_list{n};
    fsl_maths = ['fslmaths harvardoxford-cortical_prob_' mask '.nii.gz -bin harvardoxford-cortical_prob_' mask '_bin.nii.gz'];
    system(fsl_maths)
end

%% Extract data from NODDI sequences 
source_seq_list = {'r_MNI_Case1_ficvf_GM','r_MNI_Case1_fiso_GM','r_MNI_Case1_odi_GM'};

for n = 1:length(source_seq_list)
    source_seq = source_seq_list{n};
    if n ==1 
        for m = 1:length(mask_list)
            %generate voxel number, mean, and SD for each brain region
            mask = mask_list{m};
            fslstats = ['fslstats ' source_seq ' -k harvardoxford-cortical_prob_' mask '_bin.nii.gz -V -M -S >data.txt'];
            system(fslstats)
            %Tabulate data
            if m == 1
                intensities = num2cell(importdata('data.txt'));
                Tndi = cell2table (intensities, 'VariableNames', {'Voxels_ndi','Volumes_ndi','NDI Mean','SD_ndi'});

            else 
                intensities = num2cell(importdata('data.txt'));
                T2 = cell2table (intensities, 'VariableNames', {'Voxels_ndi','Volumes_ndi','NDI Mean','SD_ndi'});
                Tndi = [Tndi;T2];
            end
            
            %Extract raw intensity values per voxel for each brain region
            fsl_meants = ['fslmeants -i ' source_seq ' -m harvardoxford-cortical_prob_' mask '_bin.nii.gz -o out.txt --showall'];
            system(fsl_meants)
            system('tail -1 out.txt >out1.txt')
            if m == 1
                raw_values_1 = importdata('out1.txt');
            elseif m == 2 
                raw_values_2 = importdata('out1.txt');
            elseif m == 3 
                raw_values_3 = importdata('out1.txt');
            elseif m == 4 
                raw_values_4 = importdata('out1.txt');
            elseif m == 5 
                raw_values_5 = importdata('out1.txt');
            elseif m == 6 
                raw_values_6 = importdata('out1.txt');
            elseif m == 7 
                raw_values_7 = importdata('out1.txt');
            elseif m == 8 
                raw_values_8 = importdata('out1.txt');
            elseif m == 9 
                raw_values_9 = importdata('out1.txt');
            elseif m == 10 
                raw_values_10 = importdata('out1.txt');
            elseif m == 11 
                raw_values_11 = importdata('out1.txt');
            elseif m == 12 
                raw_values_12 = importdata('out1.txt');
            elseif m == 13 
                raw_values_13 = importdata('out1.txt');
            elseif m == 14 
                raw_values_14 = importdata('out1.txt');
            elseif m == 15 
                raw_values_15 = importdata('out1.txt');
            elseif m == 16 
                raw_values_16 = importdata('out1.txt');
            elseif m == 17 
                raw_values_17 = importdata('out1.txt');
            elseif m == 18 
                raw_values_18 = importdata('out1.txt');
            elseif m == 19 
                raw_values_19 = importdata('out1.txt');
            elseif m == 20 
                raw_values_20 = importdata('out1.txt');
            elseif m == 21 
                raw_values_21 = importdata('out1.txt');
            elseif m == 22 
                raw_values_22 = importdata('out1.txt');
            elseif m == 23
                raw_values_23 = importdata('out1.txt');
            elseif m == 24 
                raw_values_24 = importdata('out1.txt');
            elseif m == 25 
                raw_values_25 = importdata('out1.txt');
            elseif m == 26 
                raw_values_26 = importdata('out1.txt');
            elseif m == 27 
                raw_values_27 = importdata('out1.txt');
            elseif m == 28 
                raw_values_28 = importdata('out1.txt');
            elseif m == 29 
                raw_values_29 = importdata('out1.txt');
            elseif m == 30 
                raw_values_30 = importdata('out1.txt');
            elseif m == 31 
                raw_values_31 = importdata('out1.txt');
            elseif m == 32 
                raw_values_32 = importdata('out1.txt');
            elseif m == 33 
                raw_values_33 = importdata('out1.txt');
            elseif m == 34
                raw_values_34 = importdata('out1.txt');
            elseif m == 35  
                raw_values_35 = importdata('out1.txt');
            elseif m == 36 
                raw_values_36 = importdata('out1.txt');
            elseif m == 37 
                raw_values_37 = importdata('out1.txt');
            elseif m == 38 
                raw_values_38 = importdata('out1.txt');
            elseif m == 39 
                raw_values_39 = importdata('out1.txt');
            elseif m == 40 
                raw_values_40 = importdata('out1.txt');
            elseif m == 41 
                raw_values_41 = importdata('out1.txt');
            elseif m == 42 
                raw_values_42 = importdata('out1.txt');
            elseif m == 43 
                raw_values_43 = importdata('out1.txt');
            elseif m == 44 
                raw_values_44 = importdata('out1.txt');
            elseif m == 45 
                raw_values_45 = importdata('out1.txt');
            elseif m == 46 
                raw_values_46 = importdata('out1.txt');
            elseif m == 47 
                raw_values_47 = importdata('out1.txt');
            elseif m == 48 
                raw_values_48 = importdata('out1.txt');
            end
        end

        Trv_ndi = num2cell(padcat(raw_values_1,raw_values_2,raw_values_3,raw_values_4, ...
                  raw_values_5,raw_values_6,raw_values_7,raw_values_8,raw_values_9,raw_values_10,raw_values_11,raw_values_12,raw_values_13,...
                  raw_values_14,raw_values_15,raw_values_16,raw_values_17,raw_values_18,raw_values_19,raw_values_20,raw_values_21,raw_values_22,...
                  raw_values_23,raw_values_24,raw_values_25,raw_values_26,raw_values_27,raw_values_28,raw_values_29,raw_values_30,raw_values_31,...
                  raw_values_32,raw_values_33,raw_values_34,raw_values_35,raw_values_36,raw_values_37,raw_values_38,raw_values_39,raw_values_40,...
                  raw_values_41,raw_values_42,raw_values_43,raw_values_44,raw_values_45,raw_values_46,raw_values_47,raw_values_48));
        %Truncate names to max length
        mask_list_trunc = regexp(mask_list, '^.{1,50}','match','once');

        Trv_ndi = cell2table(Trv_ndi',"VariableNames", mask_list_trunc);
    
    end

    if n == 2
        for m = 1:length(mask_list)
            mask = mask_list{m};
            fslstats = ['fslstats ' source_seq ' -k harvardoxford-cortical_prob_' mask '_bin.nii.gz -V -M -S >data.txt'];
            system(fslstats)
            if m == 1
                intensities = num2cell(importdata('data.txt'));
                Tfwf = cell2table (intensities, 'VariableNames', {'Voxels_fwf','Volumes_fwf','FWF Mean','SD_fwf'});
            else 
                intensities = num2cell(importdata('data.txt'));
                T2 = cell2table (intensities, 'VariableNames', {'Voxels_fwf','Volumes_fwf','FWF Mean','SD_fwf'});
                Tfwf = [Tfwf;T2];
            end

            fsl_meants = ['fslmeants -i ' source_seq ' -m harvardoxford-cortical_prob_' mask '_bin.nii.gz -o out.txt --showall'];
            system(fsl_meants)
            system('tail -1 out.txt >out1.txt')
            if m == 1
                raw_values_1 = importdata('out1.txt');
            elseif m == 2 
                raw_values_2 = importdata('out1.txt');
            elseif m == 3 
                raw_values_3 = importdata('out1.txt');
            elseif m == 4 
                raw_values_4 = importdata('out1.txt');
            elseif m == 5 
                raw_values_5 = importdata('out1.txt');
            elseif m == 6 
                raw_values_6 = importdata('out1.txt');
            elseif m == 7 
                raw_values_7 = importdata('out1.txt');
            elseif m == 8 
                raw_values_8 = importdata('out1.txt');
            elseif m == 9 
                raw_values_9 = importdata('out1.txt');
            elseif m == 10 
                raw_values_10 = importdata('out1.txt');
            elseif m == 11 
                raw_values_11 = importdata('out1.txt');
            elseif m == 12 
                raw_values_12 = importdata('out1.txt');
            elseif m == 13 
                raw_values_13 = importdata('out1.txt');
            elseif m == 14 
                raw_values_14 = importdata('out1.txt');
            elseif m == 15 
                raw_values_15 = importdata('out1.txt');
            elseif m == 16 
                raw_values_16 = importdata('out1.txt');
            elseif m == 17 
                raw_values_17 = importdata('out1.txt');
            elseif m == 18 
                raw_values_18 = importdata('out1.txt');
            elseif m == 19 
                raw_values_19 = importdata('out1.txt');
            elseif m == 20 
                raw_values_20 = importdata('out1.txt');
            elseif m == 21 
                raw_values_21 = importdata('out1.txt');
            elseif m == 22 
                raw_values_22 = importdata('out1.txt');
            elseif m == 23
                raw_values_23 = importdata('out1.txt');
            elseif m == 24 
                raw_values_24 = importdata('out1.txt');
            elseif m == 25 
                raw_values_25 = importdata('out1.txt');
            elseif m == 26 
                raw_values_26 = importdata('out1.txt');
            elseif m == 27 
                raw_values_27 = importdata('out1.txt');
            elseif m == 28 
                raw_values_28 = importdata('out1.txt');
            elseif m == 29 
                raw_values_29 = importdata('out1.txt');
            elseif m == 30 
                raw_values_30 = importdata('out1.txt');
            elseif m == 31 
                raw_values_31 = importdata('out1.txt');
            elseif m == 32 
                raw_values_32 = importdata('out1.txt');
            elseif m == 33 
                raw_values_33 = importdata('out1.txt');
            elseif m == 34
                raw_values_34 = importdata('out1.txt');
            elseif m == 35  
                raw_values_35 = importdata('out1.txt');
            elseif m == 36 
                raw_values_36 = importdata('out1.txt');
            elseif m == 37 
                raw_values_37 = importdata('out1.txt');
            elseif m == 38 
                raw_values_38 = importdata('out1.txt');
            elseif m == 39 
                raw_values_39 = importdata('out1.txt');
            elseif m == 40 
                raw_values_40 = importdata('out1.txt');
            elseif m == 41 
                raw_values_41 = importdata('out1.txt');
            elseif m == 42 
                raw_values_42 = importdata('out1.txt');
            elseif m == 43 
                raw_values_43 = importdata('out1.txt');
            elseif m == 44 
                raw_values_44 = importdata('out1.txt');
            elseif m == 45 
                raw_values_45 = importdata('out1.txt');
            elseif m == 46 
                raw_values_46 = importdata('out1.txt');
            elseif m == 47 
                raw_values_47 = importdata('out1.txt');
            elseif m == 48 
                raw_values_48 = importdata('out1.txt');
            end
        end

        Trv_fwf = num2cell(padcat(raw_values_1,raw_values_2,raw_values_3,raw_values_4, ...
                  raw_values_5,raw_values_6,raw_values_7,raw_values_8,raw_values_9,raw_values_10,raw_values_11,raw_values_12,raw_values_13,...
                  raw_values_14,raw_values_15,raw_values_16,raw_values_17,raw_values_18,raw_values_19,raw_values_20,raw_values_21,raw_values_22,...
                  raw_values_23,raw_values_24,raw_values_25,raw_values_26,raw_values_27,raw_values_28,raw_values_29,raw_values_30,raw_values_31,...
                  raw_values_32,raw_values_33,raw_values_34,raw_values_35,raw_values_36,raw_values_37,raw_values_38,raw_values_39,raw_values_40,...
                  raw_values_41,raw_values_42,raw_values_43,raw_values_44,raw_values_45,raw_values_46,raw_values_47,raw_values_48));
        
        %Truncate names to max length
        mask_list_trunc = regexp(mask_list, '^.{1,50}','match','once');
        
        Trv_fwf = cell2table(Trv_fwf',"VariableNames", mask_list_trunc);
    end

    if n == 3
        for m = 1:length(mask_list)
            mask = mask_list{m};
            fslstats = ['fslstats ' source_seq ' -k harvardoxford-cortical_prob_' mask '_bin.nii.gz -V -M -S >data.txt'];
            system(fslstats)

            if m == 1
                intensities = num2cell(importdata('data.txt'));
                Todi = cell2table (intensities, 'VariableNames', {'Voxels_odi','Volumes_odi','ODI Mean','SD_odi'});
            else 
                intensities = num2cell(importdata('data.txt'));
                T2 = cell2table (intensities, 'VariableNames', {'Voxels_odi','Volumes_odi','ODI Mean','SD_odi'});
                Todi = [Todi;T2];
            end

            fsl_meants = ['fslmeants -i ' source_seq ' -m harvardoxford-cortical_prob_' mask '_bin.nii.gz -o out.txt --showall'];
            system(fsl_meants)
            system('tail -1 out.txt >out1.txt')
            if m == 1
                raw_values_1 = importdata('out1.txt');
            elseif m == 2 
                raw_values_2 = importdata('out1.txt');
            elseif m == 3 
                raw_values_3 = importdata('out1.txt');
            elseif m == 4 
                raw_values_4 = importdata('out1.txt');
            elseif m == 5 
                raw_values_5 = importdata('out1.txt');
            elseif m == 6 
                raw_values_6 = importdata('out1.txt');
            elseif m == 7 
                raw_values_7 = importdata('out1.txt');
            elseif m == 8 
                raw_values_8 = importdata('out1.txt');
            elseif m == 9 
                raw_values_9 = importdata('out1.txt');
            elseif m == 10 
                raw_values_10 = importdata('out1.txt');
            elseif m == 11 
                raw_values_11 = importdata('out1.txt');
            elseif m == 12 
                raw_values_12 = importdata('out1.txt');
            elseif m == 13 
                raw_values_13 = importdata('out1.txt');
            elseif m == 14 
                raw_values_14 = importdata('out1.txt');
            elseif m == 15 
                raw_values_15 = importdata('out1.txt');
            elseif m == 16 
                raw_values_16 = importdata('out1.txt');
            elseif m == 17 
                raw_values_17 = importdata('out1.txt');
            elseif m == 18 
                raw_values_18 = importdata('out1.txt');
            elseif m == 19 
                raw_values_19 = importdata('out1.txt');
            elseif m == 20 
                raw_values_20 = importdata('out1.txt');
            elseif m == 21 
                raw_values_21 = importdata('out1.txt');
            elseif m == 22 
                raw_values_22 = importdata('out1.txt');
            elseif m == 23
                raw_values_23 = importdata('out1.txt');
            elseif m == 24 
                raw_values_24 = importdata('out1.txt');
            elseif m == 25 
                raw_values_25 = importdata('out1.txt');
            elseif m == 26 
                raw_values_26 = importdata('out1.txt');
            elseif m == 27 
                raw_values_27 = importdata('out1.txt');
            elseif m == 28 
                raw_values_28 = importdata('out1.txt');
            elseif m == 29 
                raw_values_29 = importdata('out1.txt');
            elseif m == 30 
                raw_values_30 = importdata('out1.txt');
            elseif m == 31 
                raw_values_31 = importdata('out1.txt');
            elseif m == 32 
                raw_values_32 = importdata('out1.txt');
            elseif m == 33 
                raw_values_33 = importdata('out1.txt');
            elseif m == 34
                raw_values_34 = importdata('out1.txt');
            elseif m == 35  
                raw_values_35 = importdata('out1.txt');
            elseif m == 36 
                raw_values_36 = importdata('out1.txt');
            elseif m == 37 
                raw_values_37 = importdata('out1.txt');
            elseif m == 38 
                raw_values_38 = importdata('out1.txt');
            elseif m == 39 
                raw_values_39 = importdata('out1.txt');
            elseif m == 40 
                raw_values_40 = importdata('out1.txt');
            elseif m == 41 
                raw_values_41 = importdata('out1.txt');
            elseif m == 42 
                raw_values_42 = importdata('out1.txt');
            elseif m == 43 
                raw_values_43 = importdata('out1.txt');
            elseif m == 44 
                raw_values_44 = importdata('out1.txt');
            elseif m == 45 
                raw_values_45 = importdata('out1.txt');
            elseif m == 46 
                raw_values_46 = importdata('out1.txt');
            elseif m == 47 
                raw_values_47 = importdata('out1.txt');
            elseif m == 48 
                raw_values_48 = importdata('out1.txt');
            end
        end

        Trv_odi = num2cell(padcat(raw_values_1,raw_values_2,raw_values_3,raw_values_4, ...
                  raw_values_5,raw_values_6,raw_values_7,raw_values_8,raw_values_9,raw_values_10,raw_values_11,raw_values_12,raw_values_13,...
                  raw_values_14,raw_values_15,raw_values_16,raw_values_17,raw_values_18,raw_values_19,raw_values_20,raw_values_21,raw_values_22,...
                  raw_values_23,raw_values_24,raw_values_25,raw_values_26,raw_values_27,raw_values_28,raw_values_29,raw_values_30,raw_values_31,...
                  raw_values_32,raw_values_33,raw_values_34,raw_values_35,raw_values_36,raw_values_37,raw_values_38,raw_values_39,raw_values_40,...
                  raw_values_41,raw_values_42,raw_values_43,raw_values_44,raw_values_45,raw_values_46,raw_values_47,raw_values_48));
        
        %Truncate names to max length
        mask_list_trunc = regexp(mask_list, '^.{1,50}','match','once');
       
        Trv_odi = cell2table(Trv_odi',"VariableNames",mask_list_trunc);
    end
end

mask_list_key = mask_list';

key = cell2table(mask_list_key, 'VariableNames',{'Region'});

Tndi = [key, Tndi];
Todi = [key, Todi];
Tfwf = [key, Tfwf];

Tjoin = join(Tndi, Todi);
T_combined = join(Tjoin, Tfwf);
%%
writetable(T_combined, 'HO_cort_NODDI_indices.csv')
%%
writetable(Trv_ndi, 'HO_cort_NDI_raw_intensities.csv')
writetable(Trv_fwf, 'HO_cort_FWF_raw_intensities.csv')
writetable(Trv_odi, 'HO_cort_ODI_raw_intensities.csv')


%% Subcortical/White Matter Analysis

% Binarize and apply transformation matrix from MNI transformations to GM mask
source_seq = ['anat_seq_brain_mask_pve_2'];
fsl_maths = ['fslmaths ' source_seq ' -bin ' source_seq '_bin.nii.gz'];
system(fsl_maths)

source_seq = ['anat_seq_brain_mask_pve_2_bin'];
ref_seq = ['Case1_ficvf'];

flirt_coreg = ['flirt -in ' source_seq ' -ref ' ref_seq ' -out r_noddi_' source_seq ' -init invol2refvol_T1_NODDI.mat -applyxfm -v' ];
system(flirt_coreg)

source_seq = ['r_noddi_anat_seq_brain_mask_pve_2_bin'];
ref_seq = ['/mnt/c/WSL2_dir/Atlases/MNI-maxprob-thr0-1mm.nii'];

flirt_coreg = ['flirt -in ' source_seq ' -ref ' ref_seq ' -out r_MNI_' source_seq ' -init invol2refvol_T1_NODDI_MNI.mat -applyxfm -v' ];
system(flirt_coreg)


% Rebinarize mask and manually adjust threshold (higher is more stringent)
WM_mask = ['r_MNI_r_noddi_anat_seq_brain_mask_pve_2_bin'];
fsl_maths = ['fslmaths ' WM_mask ' -thr 0.55 -bin ' WM_mask '_bin.nii.gz'];
system(fsl_maths)

% Extract GM from NODDI using GM mask
WM_mask = ['r_MNI_r_noddi_anat_seq_brain_mask_pve_2_bin_bin'];
ref_seq = ['r_MNI_Case1_ficvf'];
fslmaths = ['fslmaths ' ref_seq ' -mul ' WM_mask ' ' ref_seq '_WM'];
system(fslmaths)

ref_seq = ['r_MNI_Case1_fiso'];
fslmaths = ['fslmaths ' ref_seq ' -mul ' WM_mask ' ' ref_seq '_WM'];
system(fslmaths)

ref_seq = ['r_MNI_Case1_odi'];
fslmaths = ['fslmaths ' ref_seq ' -mul ' WM_mask ' ' ref_seq '_WM'];
system(fslmaths)
%%

%% Create atlas masks on FSL %%

%% Binarize masks

%% WILL BREAK CODE -> Manually edit Heschel's gyrus and remove apostrophe %%

% Remove symbols so code won't break
mask_list = {'Left Cerebral White Matter','Left Cerebral Cortex','Left Lateral Ventricle','Left Thalamus','Left Caudate','Left Putamen',...
    'Left Pallidum','Brain-Stem','Left Hippocampus','Left Amygdala','Left Accumbens','Right Cerebral White Matter','Right Cerebral Cortex',...
    'Right Lateral Ventricle','Right Thalamus','Right Caudate','Right Putamen','Right Pallidum','Right Hippocampus','Right Amygdala',...
    'Right Accumbens'};

for n = 1:length(mask_list)
    mask = mask_list{n};
    instring = ['harvardoxford-subcortical_prob_' mask '.nii.gz'];
    outstring = regexprep(instring, ' ', '_');
    outstring = regexprep(outstring, '(','');
    outstring = regexprep(outstring, ')','');
    outstring = regexprep(outstring, ',','');

    if isfile(instring) == 0
        continue
    elseif instring == outstring
        continue
    else
        movefile(instring, outstring)
    end
end

mask_list = regexprep(mask_list, ' ', '_');
mask_list = regexprep(mask_list, '(','');
mask_list = regexprep(mask_list, ')','');
mask_list = regexprep(mask_list, ',','');

% binarize masks
for n = 1:length(mask_list)
    mask = mask_list{n};
    fsl_maths = ['fslmaths harvardoxford-subcortical_prob_' mask '.nii.gz -bin harvardoxford-subcortical_prob_' mask '_bin.nii.gz'];
    system(fsl_maths)
end

%% Extract data from NODDI sequences 
source_seq_list = {'r_MNI_Case1_ficvf_WM','r_MNI_Case1_fiso_WM','r_MNI_Case1_odi_WM'};

for n = 1:length(source_seq_list)
    source_seq = source_seq_list{n};
    if n ==1 
        for m = 1:length(mask_list)
            %generate voxel number, mean, and SD for each brain region
            mask = mask_list{m};
            fslstats = ['fslstats ' source_seq ' -k harvardoxford-subcortical_prob_' mask '_bin.nii.gz -V -M -S >data.txt'];
            system(fslstats)
            %Tabulate data
            if m == 1
                intensities = num2cell(importdata('data.txt'));
                Tndi = cell2table (intensities, 'VariableNames', {'Voxels_ndi','Volumes_ndi','NDI Mean','SD_ndi'});

            else 
                intensities = num2cell(importdata('data.txt'));
                T2 = cell2table (intensities, 'VariableNames', {'Voxels_ndi','Volumes_ndi','NDI Mean','SD_ndi'});
                Tndi = [Tndi;T2];
            end
            
            %Extract raw intensity values per voxel for each brain region
            fsl_meants = ['fslmeants -i ' source_seq ' -m harvardoxford-subcortical_prob_' mask '_bin.nii.gz -o out.txt --showall'];
            system(fsl_meants)
            system('tail -1 out.txt >out1.txt')
            if m == 1
                raw_values_1 = importdata('out1.txt');
            elseif m == 2 
                raw_values_2 = importdata('out1.txt');
            elseif m == 3 
                raw_values_3 = importdata('out1.txt');
            elseif m == 4 
                raw_values_4 = importdata('out1.txt');
            elseif m == 5 
                raw_values_5 = importdata('out1.txt');
            elseif m == 6 
                raw_values_6 = importdata('out1.txt');
            elseif m == 7 
                raw_values_7 = importdata('out1.txt');
            elseif m == 8 
                raw_values_8 = importdata('out1.txt');
            elseif m == 9 
                raw_values_9 = importdata('out1.txt');
            elseif m == 10 
                raw_values_10 = importdata('out1.txt');
            elseif m == 11 
                raw_values_11 = importdata('out1.txt');
            elseif m == 12 
                raw_values_12 = importdata('out1.txt');
            elseif m == 13 
                raw_values_13 = importdata('out1.txt');
            elseif m == 14 
                raw_values_14 = importdata('out1.txt');
            elseif m == 15 
                raw_values_15 = importdata('out1.txt');
            elseif m == 16 
                raw_values_16 = importdata('out1.txt');
            elseif m == 17 
                raw_values_17 = importdata('out1.txt');
            elseif m == 18 
                raw_values_18 = importdata('out1.txt');
            elseif m == 19 
                raw_values_19 = importdata('out1.txt');
            elseif m == 20 
                raw_values_20 = importdata('out1.txt');
            elseif m == 21 
                raw_values_21 = importdata('out1.txt');
            end
        end

        Trv_ndi = num2cell(padcat(raw_values_1,raw_values_2,raw_values_3,raw_values_4, ...
                  raw_values_5,raw_values_6,raw_values_7,raw_values_8,raw_values_9,raw_values_10,raw_values_11,raw_values_12,raw_values_13,...
                  raw_values_14,raw_values_15,raw_values_16,raw_values_17,raw_values_18,raw_values_19,raw_values_20,raw_values_21));
        %Truncate names to max length
        mask_list_trunc = regexp(mask_list, '^.{1,50}','match','once');

        Trv_ndi = cell2table(Trv_ndi',"VariableNames", mask_list_trunc);
    
    end

    if n == 2
        for m = 1:length(mask_list)
            mask = mask_list{m};
            fslstats = ['fslstats ' source_seq ' -k harvardoxford-subcortical_prob_' mask '_bin.nii.gz -V -M -S >data.txt'];
            system(fslstats)
            if m == 1
                intensities = num2cell(importdata('data.txt'));
                Tfwf = cell2table (intensities, 'VariableNames', {'Voxels_fwf','Volumes_fwf','FWF Mean','SD_fwf'});
            else 
                intensities = num2cell(importdata('data.txt'));
                T2 = cell2table (intensities, 'VariableNames', {'Voxels_fwf','Volumes_fwf','FWF Mean','SD_fwf'});
                Tfwf = [Tfwf;T2];
            end

            fsl_meants = ['fslmeants -i ' source_seq ' -m harvardoxford-subcortical_prob_' mask '_bin.nii.gz -o out.txt --showall'];
            system(fsl_meants)
            system('tail -1 out.txt >out1.txt')
            if m == 1
                raw_values_1 = importdata('out1.txt');
            elseif m == 2 
                raw_values_2 = importdata('out1.txt');
            elseif m == 3 
                raw_values_3 = importdata('out1.txt');
            elseif m == 4 
                raw_values_4 = importdata('out1.txt');
            elseif m == 5 
                raw_values_5 = importdata('out1.txt');
            elseif m == 6 
                raw_values_6 = importdata('out1.txt');
            elseif m == 7 
                raw_values_7 = importdata('out1.txt');
            elseif m == 8 
                raw_values_8 = importdata('out1.txt');
            elseif m == 9 
                raw_values_9 = importdata('out1.txt');
            elseif m == 10 
                raw_values_10 = importdata('out1.txt');
            elseif m == 11 
                raw_values_11 = importdata('out1.txt');
            elseif m == 12 
                raw_values_12 = importdata('out1.txt');
            elseif m == 13 
                raw_values_13 = importdata('out1.txt');
            elseif m == 14 
                raw_values_14 = importdata('out1.txt');
            elseif m == 15 
                raw_values_15 = importdata('out1.txt');
            elseif m == 16 
                raw_values_16 = importdata('out1.txt');
            elseif m == 17 
                raw_values_17 = importdata('out1.txt');
            elseif m == 18 
                raw_values_18 = importdata('out1.txt');
            elseif m == 19 
                raw_values_19 = importdata('out1.txt');
            elseif m == 20 
                raw_values_20 = importdata('out1.txt');
            elseif m == 21 
                raw_values_21 = importdata('out1.txt');
            end
        end

        Trv_fwf = num2cell(padcat(raw_values_1,raw_values_2,raw_values_3,raw_values_4, ...
                  raw_values_5,raw_values_6,raw_values_7,raw_values_8,raw_values_9,raw_values_10,raw_values_11,raw_values_12,raw_values_13,...
                  raw_values_14,raw_values_15,raw_values_16,raw_values_17,raw_values_18,raw_values_19,raw_values_20,raw_values_21));
        
        %Truncate names to max length
        mask_list_trunc = regexp(mask_list, '^.{1,50}','match','once');
        
        Trv_fwf = cell2table(Trv_fwf',"VariableNames", mask_list_trunc);
    end

    if n == 3
        for m = 1:length(mask_list)
            mask = mask_list{m};
            fslstats = ['fslstats ' source_seq ' -k harvardoxford-subcortical_prob_' mask '_bin.nii.gz -V -M -S >data.txt'];
            system(fslstats)

            if m == 1
                intensities = num2cell(importdata('data.txt'));
                Todi = cell2table (intensities, 'VariableNames', {'Voxels_odi','Volumes_odi','ODI Mean','SD_odi'});
            else 
                intensities = num2cell(importdata('data.txt'));
                T2 = cell2table (intensities, 'VariableNames', {'Voxels_odi','Volumes_odi','ODI Mean','SD_odi'});
                Todi = [Todi;T2];
            end

            fsl_meants = ['fslmeants -i ' source_seq ' -m harvardoxford-subcortical_prob_' mask '_bin.nii.gz -o out.txt --showall'];
            system(fsl_meants)
            system('tail -1 out.txt >out1.txt')
            if m == 1
                raw_values_1 = importdata('out1.txt');
            elseif m == 2 
                raw_values_2 = importdata('out1.txt');
            elseif m == 3 
                raw_values_3 = importdata('out1.txt');
            elseif m == 4 
                raw_values_4 = importdata('out1.txt');
            elseif m == 5 
                raw_values_5 = importdata('out1.txt');
            elseif m == 6 
                raw_values_6 = importdata('out1.txt');
            elseif m == 7 
                raw_values_7 = importdata('out1.txt');
            elseif m == 8 
                raw_values_8 = importdata('out1.txt');
            elseif m == 9 
                raw_values_9 = importdata('out1.txt');
            elseif m == 10 
                raw_values_10 = importdata('out1.txt');
            elseif m == 11 
                raw_values_11 = importdata('out1.txt');
            elseif m == 12 
                raw_values_12 = importdata('out1.txt');
            elseif m == 13 
                raw_values_13 = importdata('out1.txt');
            elseif m == 14 
                raw_values_14 = importdata('out1.txt');
            elseif m == 15 
                raw_values_15 = importdata('out1.txt');
            elseif m == 16 
                raw_values_16 = importdata('out1.txt');
            elseif m == 17 
                raw_values_17 = importdata('out1.txt');
            elseif m == 18 
                raw_values_18 = importdata('out1.txt');
            elseif m == 19 
                raw_values_19 = importdata('out1.txt');
            elseif m == 20 
                raw_values_20 = importdata('out1.txt');
            elseif m == 21 
                raw_values_21 = importdata('out1.txt');
            end
        end

        Trv_odi = num2cell(padcat(raw_values_1,raw_values_2,raw_values_3,raw_values_4, ...
                  raw_values_5,raw_values_6,raw_values_7,raw_values_8,raw_values_9,raw_values_10,raw_values_11,raw_values_12,raw_values_13,...
                  raw_values_14,raw_values_15,raw_values_16,raw_values_17,raw_values_18,raw_values_19,raw_values_20,raw_values_21));
        
        %Truncate names to max length
        mask_list_trunc = regexp(mask_list, '^.{1,50}','match','once');
       
        Trv_odi = cell2table(Trv_odi',"VariableNames",mask_list_trunc);
    end
end

mask_list_key = mask_list';

key = cell2table(mask_list_key, 'VariableNames',{'Region'});

Tndi = [key, Tndi];
Todi = [key, Todi];
Tfwf = [key, Tfwf];

Tjoin = join(Tndi, Todi);
T_combined = join(Tjoin, Tfwf);
%%
writetable(T_combined, 'HO_subcort_NODDI_indices.csv')
%%
writetable(Trv_ndi, 'HO_subcort_NDI_raw_intensities.csv')
writetable(Trv_fwf, 'HO_subcort_FWF_raw_intensities.csv')
writetable(Trv_odi, 'HO_subcort_ODI_raw_intensities.csv')


%% John Hopkins Analysis (DTI based atlases) %%

%% Labels Analysis

%%
% coregister DTI tensor to JHU atlas
source_seq = ['dti_FA'];
ref_seq = ['/mnt/c/WSL2_dir/Atlases/JHU-ICBM-FA-1mm'];
flirt_coreg = ['flirt -in ' source_seq ' -ref ' ref_seq ' -out r_JHU_' source_seq ' -omat invol2refvol_DTI_JHU.mat -v' ];
system(flirt_coreg)

% coregister T1 MPRAGE to DTI tensor
source_seq = ['anat_seq_brain_mask'];
ref_seq = ['dti_FA'];
flirt_coreg = ['flirt -in ' source_seq ' -ref ' ref_seq ' -out r_dti_' source_seq ' -omat invol2refvol_T1_DTI.mat -v' ];
system(flirt_coreg)

% apply transformation matrix to WM mask
source_seq = ['anat_seq_brain_mask_pve_2'];
ref_seq = ['dti_FA'];
flirt_coreg = ['flirt -in ' source_seq ' -ref ' ref_seq ' -out r_dti_' source_seq ' -init invol2refvol_T1_DTI.mat -applyxfm -v' ];
system(flirt_coreg)

source_seq = ['r_dti_anat_seq_brain_mask_pve_2'];
ref_seq = ['/mnt/c/WSL2_dir/Atlases/JHU-ICBM-FA-1mm'];
flirt_coreg = ['flirt -in ' source_seq ' -ref ' ref_seq ' -out r_JHU_' source_seq ' -init invol2refvol_DTI_JHU.mat -applyxfm -v' ];
system(flirt_coreg)

% apply transformation matrix to NODDI 
ref_seq = ['/mnt/c/WSL2_dir/Atlases/JHU-ICBM-FA-1mm'];

source_seq = ['Case1_ficvf'];
flirt_coreg = ['flirt -in ' source_seq ' -ref ' ref_seq ' -out r_JHU_' source_seq ' -init invol2refvol_DTI_JHU.mat -applyxfm -v' ];
system(flirt_coreg)

source_seq = ['Case1_fiso'];
flirt_coreg = ['flirt -in ' source_seq ' -ref ' ref_seq ' -out r_JHU_' source_seq ' -init invol2refvol_DTI_JHU.mat -applyxfm -v' ];
system(flirt_coreg)

source_seq = ['Case1_odi'];
flirt_coreg = ['flirt -in ' source_seq ' -ref ' ref_seq ' -out r_JHU_' source_seq ' -init invol2refvol_DTI_JHU.mat -applyxfm -v' ];
system(flirt_coreg)

% Extract WM from NODDI using WM mask
WM_mask = ['r_JHU_r_dti_anat_seq_brain_mask_pve_2'];
fslmaths = ['fslmaths ' WM_mask ' -thr 0.4 -bin ' WM_mask '_bin.nii.gz'];
system(fslmaths)

ref_seq = ['r_JHU_Case1_ficvf'];
fslmaths = ['fslmaths ' ref_seq ' -mul ' WM_mask '_bin ' ref_seq '_WM'];
system(fslmaths)

ref_seq = ['r_JHU_Case1_fiso'];
fslmaths = ['fslmaths ' ref_seq ' -mul ' WM_mask '_bin ' ref_seq '_WM'];
system(fslmaths)

ref_seq = ['r_JHU_Case1_odi'];
fslmaths = ['fslmaths ' ref_seq ' -mul ' WM_mask '_bin ' ref_seq '_WM'];
system(fslmaths)
%%

%% Create atlas masks on FSL %%

%% Binarize masks

%% WILL BREAK CODE -> Manually cut short file name for fornix cres, superior frontal occipital fasciculus, sagittal stratum  %%

% Remove symbols so code won't break
mask_list = {'Anterior corona radiata L','Anterior corona radiata R','Superior corona radiata R','Superior corona radiata L',...
    'Retrolenticular part of internal capsule L','Retrolenticular part of internal capsule R',...
    'Posterior limb of internal capsule L','Posterior limb of internal capsule R','Anterior limb of internal capsule L',...
    'Anterior limb of internal capsule R','Cerebral peduncle L','Cerebral peduncle R','Superior cerebellar peduncle L','Superior cerebellar peduncle R',...
    'Inferior cerebellar peduncle L','Inferior cerebellar peduncle R','Medial lemniscus L','Medial lemniscus R','Corticospinal tract L',...
    'Corticospinal tract R','Fornix (column and body of fornix)','Splenium of corpus callosum','Body of corpus callosum','Genu of corpus callosum',...
    'Pontine crossing tract (a part of MCP)','Middle cerebellar peduncle','Tapetum L','Tapetum R','Uncinate fasciculus L','Uncinate fasciculus R'...
    'Inferior fronto-occipital fasciculus L','Inferior fronto-occipital fasciculus R','Superior fronto-occipital fasciculus L',...
    'Superior fronto-occipital fasciculus R','Superior longitudinal fasciculus L','Superior longitudinal fasciculus R',...
    'Fornix (cres) L','Fornix (cres) R',...
    'Cingulum (hippocampus) L','Cingulum (hippocampus) R','Cingulum (cingulate gyrus) L','Cingulum (cingulate gyrus) R','External capsule L',...
    'External capsule R','Sagittal stratum R',...
    'Sagittal stratum L',...
    'Posterior thalamic radiation (include optic radiation) L','Posterior thalamic radiation (include optic radiation) R'...
    'Posterior corona radiata L','Posterior corona radiata R'};

for n = 1:length(mask_list)
    mask = mask_list{n};
    instring = ['jhu-labels_label_' mask '.nii.gz'];
    outstring = regexprep(instring, ' ', '_');
    outstring = regexprep(outstring, '(','');
    outstring = regexprep(outstring, ')','');
    outstring = regexprep(outstring, ',','');
    if isfile(instring) == 0
        continue
    elseif instring == outstring
        continue
    else
        movefile(instring, outstring)
    end

end

mask_list = regexprep(mask_list, ' ', '_');
mask_list = regexprep(mask_list, '(','');
mask_list = regexprep(mask_list, ')','');
mask_list = regexprep(mask_list, ',','');

% binarize masks
for n = 1:length(mask_list)
    mask = mask_list{n};
    fsl_maths = ['fslmaths jhu-labels_label_' mask '.nii.gz -bin jhu-labels_label_' mask '_bin.nii.gz'];
    system(fsl_maths)
end

%% Extract data from NODDI sequences 
source_seq_list = {'r_JHU_Case1_ficvf_WM','r_JHU_Case1_fiso_WM','r_JHU_Case1_odi_WM'};

for n = 1:length(source_seq_list)
    source_seq = source_seq_list{n};
    if n ==1 
        for m = 1:length(mask_list)
            %generate voxel number, mean, and SD for each brain region
            mask = mask_list{m};
            fslstats = ['fslstats ' source_seq ' -k jhu-labels_label_' mask '_bin.nii.gz -V -M -S >data.txt'];
            system(fslstats)
            %Tabulate data
            if m == 1
                intensities = num2cell(importdata('data.txt'));
                Tndi = cell2table (intensities, 'VariableNames', {'Voxels_ndi','Volumes_ndi','NDI Mean','SD_ndi'});

            else 
                intensities = num2cell(importdata('data.txt'));
                T2 = cell2table (intensities, 'VariableNames', {'Voxels_ndi','Volumes_ndi','NDI Mean','SD_ndi'});
                Tndi = [Tndi;T2];
            end
            
            %Extract raw intensity values per voxel for each brain region
            fsl_meants = ['fslmeants -i ' source_seq ' -m jhu-labels_label_' mask '_bin.nii.gz -o out.txt --showall'];
            system(fsl_meants)
            system('tail -1 out.txt >out1.txt')
            if m == 1
                raw_values_1 = importdata('out1.txt');
            elseif m == 2 
                raw_values_2 = importdata('out1.txt');
            elseif m == 3 
                raw_values_3 = importdata('out1.txt');
            elseif m == 4 
                raw_values_4 = importdata('out1.txt');
            elseif m == 5 
                raw_values_5 = importdata('out1.txt');
            elseif m == 6 
                raw_values_6 = importdata('out1.txt');
            elseif m == 7 
                raw_values_7 = importdata('out1.txt');
            elseif m == 8 
                raw_values_8 = importdata('out1.txt');
            elseif m == 9 
                raw_values_9 = importdata('out1.txt');
            elseif m == 10 
                raw_values_10 = importdata('out1.txt');
            elseif m == 11 
                raw_values_11 = importdata('out1.txt');
            elseif m == 12 
                raw_values_12 = importdata('out1.txt');
            elseif m == 13 
                raw_values_13 = importdata('out1.txt');
            elseif m == 14 
                raw_values_14 = importdata('out1.txt');
            elseif m == 15 
                raw_values_15 = importdata('out1.txt');
            elseif m == 16 
                raw_values_16 = importdata('out1.txt');
            elseif m == 17 
                raw_values_17 = importdata('out1.txt');
            elseif m == 18 
                raw_values_18 = importdata('out1.txt');
            elseif m == 19 
                raw_values_19 = importdata('out1.txt');
            elseif m == 20 
                raw_values_20 = importdata('out1.txt');
            elseif m == 21 
                raw_values_21 = importdata('out1.txt');
            elseif m == 22 
                raw_values_22 = importdata('out1.txt');
            elseif m == 23
                raw_values_23 = importdata('out1.txt');
            elseif m == 24 
                raw_values_24 = importdata('out1.txt');
            elseif m == 25 
                raw_values_25 = importdata('out1.txt');
            elseif m == 26 
                raw_values_26 = importdata('out1.txt');
            elseif m == 27 
                raw_values_27 = importdata('out1.txt');
            elseif m == 28 
                raw_values_28 = importdata('out1.txt');
            elseif m == 29 
                raw_values_29 = importdata('out1.txt');
            elseif m == 30 
                raw_values_30 = importdata('out1.txt');
            elseif m == 31 
                raw_values_31 = importdata('out1.txt');
            elseif m == 32 
                raw_values_32 = importdata('out1.txt');
            elseif m == 33 
                raw_values_33 = importdata('out1.txt');
            elseif m == 34
                raw_values_34 = importdata('out1.txt');
            elseif m == 35  
                raw_values_35 = importdata('out1.txt');
            elseif m == 36 
                raw_values_36 = importdata('out1.txt');
            elseif m == 37 
                raw_values_37 = importdata('out1.txt');
            elseif m == 38 
                raw_values_38 = importdata('out1.txt');
            elseif m == 39 
                raw_values_39 = importdata('out1.txt');
            elseif m == 40 
                raw_values_40 = importdata('out1.txt');
            elseif m == 41 
                raw_values_41 = importdata('out1.txt');
            elseif m == 42 
                raw_values_42 = importdata('out1.txt');
            elseif m == 43 
                raw_values_43 = importdata('out1.txt');
            elseif m == 44 
                raw_values_44 = importdata('out1.txt');
            elseif m == 45 
                raw_values_45 = importdata('out1.txt');
            elseif m == 46 
                raw_values_46 = importdata('out1.txt');
            elseif m == 47 
                raw_values_47 = importdata('out1.txt');
            elseif m == 48 
                raw_values_48 = importdata('out1.txt');
            elseif m == 49 
                raw_values_49 = importdata('out1.txt');
            elseif m == 50 
                raw_values_50 = importdata('out1.txt');

            end
        end

        Trv_ndi = num2cell(padcat(raw_values_1,raw_values_2,raw_values_3,raw_values_4, ...
                  raw_values_5,raw_values_6,raw_values_7,raw_values_8,raw_values_9,raw_values_10,raw_values_11,raw_values_12,raw_values_13,...
                  raw_values_14,raw_values_15,raw_values_16,raw_values_17,raw_values_18,raw_values_19,raw_values_20,raw_values_21,raw_values_22,...
                  raw_values_23,raw_values_24,raw_values_25,raw_values_26,raw_values_27,raw_values_28,raw_values_29,raw_values_30,raw_values_31,...
                  raw_values_32,raw_values_33,raw_values_34,raw_values_35,raw_values_36,raw_values_37,raw_values_38,raw_values_39,raw_values_40,...
                  raw_values_41,raw_values_42,raw_values_43,raw_values_44,raw_values_45,raw_values_46,raw_values_47,raw_values_48,raw_values_49,...
                  raw_values_50));
        %Truncate names to max length
        mask_list_trunc = regexp(mask_list, '^.{1,63}','match','once');

        Trv_ndi = cell2table(Trv_ndi',"VariableNames", mask_list_trunc);
    
    end

    if n == 2
        for m = 1:length(mask_list)
            mask = mask_list{m};
            fslstats = ['fslstats ' source_seq ' -k jhu-labels_label_' mask '_bin.nii.gz -V -M -S >data.txt'];
            system(fslstats)
            if m == 1
                intensities = num2cell(importdata('data.txt'));
                Tfwf = cell2table (intensities, 'VariableNames', {'Voxels_fwf','Volumes_fwf','FWF Mean','SD_fwf'});
            else 
                intensities = num2cell(importdata('data.txt'));
                T2 = cell2table (intensities, 'VariableNames', {'Voxels_fwf','Volumes_fwf','FWF Mean','SD_fwf'});
                Tfwf = [Tfwf;T2];
            end

            fsl_meants = ['fslmeants -i ' source_seq ' -m jhu-labels_label_' mask '_bin.nii.gz -o out.txt --showall'];
            system(fsl_meants)
            system('tail -1 out.txt >out1.txt')
            if m == 1
                raw_values_1 = importdata('out1.txt');
            elseif m == 2 
                raw_values_2 = importdata('out1.txt');
            elseif m == 3 
                raw_values_3 = importdata('out1.txt');
            elseif m == 4 
                raw_values_4 = importdata('out1.txt');
            elseif m == 5 
                raw_values_5 = importdata('out1.txt');
            elseif m == 6 
                raw_values_6 = importdata('out1.txt');
            elseif m == 7 
                raw_values_7 = importdata('out1.txt');
            elseif m == 8 
                raw_values_8 = importdata('out1.txt');
            elseif m == 9 
                raw_values_9 = importdata('out1.txt');
            elseif m == 10 
                raw_values_10 = importdata('out1.txt');
            elseif m == 11 
                raw_values_11 = importdata('out1.txt');
            elseif m == 12 
                raw_values_12 = importdata('out1.txt');
            elseif m == 13 
                raw_values_13 = importdata('out1.txt');
            elseif m == 14 
                raw_values_14 = importdata('out1.txt');
            elseif m == 15 
                raw_values_15 = importdata('out1.txt');
            elseif m == 16 
                raw_values_16 = importdata('out1.txt');
            elseif m == 17 
                raw_values_17 = importdata('out1.txt');
            elseif m == 18 
                raw_values_18 = importdata('out1.txt');
            elseif m == 19 
                raw_values_19 = importdata('out1.txt');
            elseif m == 20 
                raw_values_20 = importdata('out1.txt');
            elseif m == 21 
                raw_values_21 = importdata('out1.txt');
            elseif m == 22 
                raw_values_22 = importdata('out1.txt');
            elseif m == 23
                raw_values_23 = importdata('out1.txt');
            elseif m == 24 
                raw_values_24 = importdata('out1.txt');
            elseif m == 25 
                raw_values_25 = importdata('out1.txt');
            elseif m == 26 
                raw_values_26 = importdata('out1.txt');
            elseif m == 27 
                raw_values_27 = importdata('out1.txt');
            elseif m == 28 
                raw_values_28 = importdata('out1.txt');
            elseif m == 29 
                raw_values_29 = importdata('out1.txt');
            elseif m == 30 
                raw_values_30 = importdata('out1.txt');
            elseif m == 31 
                raw_values_31 = importdata('out1.txt');
            elseif m == 32 
                raw_values_32 = importdata('out1.txt');
            elseif m == 33 
                raw_values_33 = importdata('out1.txt');
            elseif m == 34
                raw_values_34 = importdata('out1.txt');
            elseif m == 35  
                raw_values_35 = importdata('out1.txt');
            elseif m == 36 
                raw_values_36 = importdata('out1.txt');
            elseif m == 37 
                raw_values_37 = importdata('out1.txt');
            elseif m == 38 
                raw_values_38 = importdata('out1.txt');
            elseif m == 39 
                raw_values_39 = importdata('out1.txt');
            elseif m == 40 
                raw_values_40 = importdata('out1.txt');
            elseif m == 41 
                raw_values_41 = importdata('out1.txt');
            elseif m == 42 
                raw_values_42 = importdata('out1.txt');
            elseif m == 43 
                raw_values_43 = importdata('out1.txt');
            elseif m == 44 
                raw_values_44 = importdata('out1.txt');
            elseif m == 45 
                raw_values_45 = importdata('out1.txt');
            elseif m == 46 
                raw_values_46 = importdata('out1.txt');
            elseif m == 47 
                raw_values_47 = importdata('out1.txt');
            elseif m == 48 
                raw_values_48 = importdata('out1.txt');
            elseif m == 49 
                raw_values_49 = importdata('out1.txt');
            elseif m == 50 
                raw_values_50 = importdata('out1.txt');
            end
        end

        Trv_fwf = num2cell(padcat(raw_values_1,raw_values_2,raw_values_3,raw_values_4, ...
                  raw_values_5,raw_values_6,raw_values_7,raw_values_8,raw_values_9,raw_values_10,raw_values_11,raw_values_12,raw_values_13,...
                  raw_values_14,raw_values_15,raw_values_16,raw_values_17,raw_values_18,raw_values_19,raw_values_20,raw_values_21,raw_values_22,...
                  raw_values_23,raw_values_24,raw_values_25,raw_values_26,raw_values_27,raw_values_28,raw_values_29,raw_values_30,raw_values_31,...
                  raw_values_32,raw_values_33,raw_values_34,raw_values_35,raw_values_36,raw_values_37,raw_values_38,raw_values_39,raw_values_40,...
                  raw_values_41,raw_values_42,raw_values_43,raw_values_44,raw_values_45,raw_values_46,raw_values_47,raw_values_48,raw_values_49,...
                  raw_values_50));
        
        %Truncate names to max length
        mask_list_trunc = regexp(mask_list, '^.{1,63}','match','once');
        
        Trv_fwf = cell2table(Trv_fwf',"VariableNames", mask_list_trunc);
    end

    if n == 3
        for m = 1:length(mask_list)
            mask = mask_list{m};
            fslstats = ['fslstats ' source_seq ' -k jhu-labels_label_' mask '_bin.nii.gz -V -M -S >data.txt'];
            system(fslstats)

            if m == 1
                intensities = num2cell(importdata('data.txt'));
                Todi = cell2table (intensities, 'VariableNames', {'Voxels_odi','Volumes_odi','ODI Mean','SD_odi'});
            else 
                intensities = num2cell(importdata('data.txt'));
                T2 = cell2table (intensities, 'VariableNames', {'Voxels_odi','Volumes_odi','ODI Mean','SD_odi'});
                Todi = [Todi;T2];
            end

            fsl_meants = ['fslmeants -i ' source_seq ' -m jhu-labels_label_' mask '_bin.nii.gz -o out.txt --showall'];
            system(fsl_meants)
            system('tail -1 out.txt >out1.txt')
            if m == 1
                raw_values_1 = importdata('out1.txt');
            elseif m == 2 
                raw_values_2 = importdata('out1.txt');
            elseif m == 3 
                raw_values_3 = importdata('out1.txt');
            elseif m == 4 
                raw_values_4 = importdata('out1.txt');
            elseif m == 5 
                raw_values_5 = importdata('out1.txt');
            elseif m == 6 
                raw_values_6 = importdata('out1.txt');
            elseif m == 7 
                raw_values_7 = importdata('out1.txt');
            elseif m == 8 
                raw_values_8 = importdata('out1.txt');
            elseif m == 9 
                raw_values_9 = importdata('out1.txt');
            elseif m == 10 
                raw_values_10 = importdata('out1.txt');
            elseif m == 11 
                raw_values_11 = importdata('out1.txt');
            elseif m == 12 
                raw_values_12 = importdata('out1.txt');
            elseif m == 13 
                raw_values_13 = importdata('out1.txt');
            elseif m == 14 
                raw_values_14 = importdata('out1.txt');
            elseif m == 15 
                raw_values_15 = importdata('out1.txt');
            elseif m == 16 
                raw_values_16 = importdata('out1.txt');
            elseif m == 17 
                raw_values_17 = importdata('out1.txt');
            elseif m == 18 
                raw_values_18 = importdata('out1.txt');
            elseif m == 19 
                raw_values_19 = importdata('out1.txt');
            elseif m == 20 
                raw_values_20 = importdata('out1.txt');
            elseif m == 21 
                raw_values_21 = importdata('out1.txt');
            elseif m == 22 
                raw_values_22 = importdata('out1.txt');
            elseif m == 23
                raw_values_23 = importdata('out1.txt');
            elseif m == 24 
                raw_values_24 = importdata('out1.txt');
            elseif m == 25 
                raw_values_25 = importdata('out1.txt');
            elseif m == 26 
                raw_values_26 = importdata('out1.txt');
            elseif m == 27 
                raw_values_27 = importdata('out1.txt');
            elseif m == 28 
                raw_values_28 = importdata('out1.txt');
            elseif m == 29 
                raw_values_29 = importdata('out1.txt');
            elseif m == 30 
                raw_values_30 = importdata('out1.txt');
            elseif m == 31 
                raw_values_31 = importdata('out1.txt');
            elseif m == 32 
                raw_values_32 = importdata('out1.txt');
            elseif m == 33 
                raw_values_33 = importdata('out1.txt');
            elseif m == 34
                raw_values_34 = importdata('out1.txt');
            elseif m == 35  
                raw_values_35 = importdata('out1.txt');
            elseif m == 36 
                raw_values_36 = importdata('out1.txt');
            elseif m == 37 
                raw_values_37 = importdata('out1.txt');
            elseif m == 38 
                raw_values_38 = importdata('out1.txt');
            elseif m == 39 
                raw_values_39 = importdata('out1.txt');
            elseif m == 40 
                raw_values_40 = importdata('out1.txt');
            elseif m == 41 
                raw_values_41 = importdata('out1.txt');
            elseif m == 42 
                raw_values_42 = importdata('out1.txt');
            elseif m == 43 
                raw_values_43 = importdata('out1.txt');
            elseif m == 44 
                raw_values_44 = importdata('out1.txt');
            elseif m == 45 
                raw_values_45 = importdata('out1.txt');
            elseif m == 46 
                raw_values_46 = importdata('out1.txt');
            elseif m == 47 
                raw_values_47 = importdata('out1.txt');
            elseif m == 48 
                raw_values_48 = importdata('out1.txt');
            elseif m == 49 
                raw_values_49 = importdata('out1.txt');
            elseif m == 50 
                raw_values_50 = importdata('out1.txt');
            end
        end

        Trv_odi = num2cell(padcat(raw_values_1,raw_values_2,raw_values_3,raw_values_4, ...
                  raw_values_5,raw_values_6,raw_values_7,raw_values_8,raw_values_9,raw_values_10,raw_values_11,raw_values_12,raw_values_13,...
                  raw_values_14,raw_values_15,raw_values_16,raw_values_17,raw_values_18,raw_values_19,raw_values_20,raw_values_21,raw_values_22,...
                  raw_values_23,raw_values_24,raw_values_25,raw_values_26,raw_values_27,raw_values_28,raw_values_29,raw_values_30,raw_values_31,...
                  raw_values_32,raw_values_33,raw_values_34,raw_values_35,raw_values_36,raw_values_37,raw_values_38,raw_values_39,raw_values_40,...
                  raw_values_41,raw_values_42,raw_values_43,raw_values_44,raw_values_45,raw_values_46,raw_values_47,raw_values_48,raw_values_49,...
                  raw_values_50));
        
        %Truncate names to max length
        mask_list_trunc = regexp(mask_list, '^.{1,63}','match','once');
       
        Trv_odi = cell2table(Trv_odi',"VariableNames",mask_list_trunc);
    end
end

mask_list_key = mask_list';

key = cell2table(mask_list_key, 'VariableNames',{'Region'});

Tndi = [key, Tndi];
Todi = [key, Todi];
Tfwf = [key, Tfwf];

Tjoin = join(Tndi, Todi);
T_combined = join(Tjoin, Tfwf);
%%
writetable(T_combined, 'JHU_label_NODDI_indices.csv')
%%
writetable(Trv_ndi, 'JHU_label_NDI_raw_intensities.csv')
writetable(Trv_fwf, 'JHU_label_FWF_raw_intensities.csv')
writetable(Trv_odi, 'JHU_label_ODI_raw_intensities.csv')


%% WM Tractography Analysis

%% Use coregistered NODDI files from previous step

%% Create atlas masks on FSL %%

%% Binarize masks

%% WILL BREAK CODE -> Manually cut short file name for fornix cres, superior frontal occipital fasciculus, sagittal stratum  %%

% Remove symbols so code won't break
mask_list = {'Superior longitudinal fasciculus (temporal part) R','Superior longitudinal fasciculus (temporal part) L','Uncinate fasciculus R',...
    'Uncinate fasciculus L','Superior longitudinal fasciculus R','Superior longitudinal fasciculus L',...
    'Inferior longitudinal fasciculus R','Inferior longitudinal fasciculus L','Inferior fronto-occipital fasciculus R','Inferior fronto-occipital fasciculus L',...
    'Forceps minor','Forceps major','Cingulum (hippocampus) R','Cingulum (hippocampus) L','Cingulum (cingulate gyrus) R','Cingulum (cingulate gyrus) L',...
    'Corticospinal tract R','Corticospinal tract L','Anterior thalamic radiation R','Anterior thalamic radiation L'};

for n = 1:length(mask_list)
    mask = mask_list{n};
    instring = ['jhu-tracts_prob_' mask '.nii.gz'];
    outstring = regexprep(instring, ' ', '_');
    outstring = regexprep(outstring, '(','');
    outstring = regexprep(outstring, ')','');
    outstring = regexprep(outstring, ',','');
    if isfile(instring) == 0
        continue
    else
        movefile(instring, outstring)
    end

end

mask_list = regexprep(mask_list, ' ', '_');
mask_list = regexprep(mask_list, '(','');
mask_list = regexprep(mask_list, ')','');
mask_list = regexprep(mask_list, ',','');

% binarize masks
for n = 1:length(mask_list)
    mask = mask_list{n};
    fsl_maths = ['fslmaths jhu-tracts_prob_' mask '.nii.gz -bin jhu-tracts_prob_' mask '_bin.nii.gz'];
    system(fsl_maths)
end

%% Extract data from NODDI sequences 
source_seq_list = {'r_JHU_Case1_ficvf_WM','r_JHU_Case1_fiso_WM','r_JHU_Case1_odi_WM'};

for n = 1:length(source_seq_list)
    source_seq = source_seq_list{n};
    if n ==1 
        for m = 1:length(mask_list)
            %generate voxel number, mean, and SD for each brain region
            mask = mask_list{m};
            fslstats = ['fslstats ' source_seq ' -k jhu-tracts_prob_' mask '_bin.nii.gz -V -M -S >data.txt'];
            system(fslstats)
            %Tabulate data
            if m == 1
                intensities = num2cell(importdata('data.txt'));
                Tndi = cell2table (intensities, 'VariableNames', {'Voxels_ndi','Volumes_ndi','NDI Mean','SD_ndi'});

            else 
                intensities = num2cell(importdata('data.txt'));
                T2 = cell2table (intensities, 'VariableNames', {'Voxels_ndi','Volumes_ndi','NDI Mean','SD_ndi'});
                Tndi = [Tndi;T2];
            end
            
            %Extract raw intensity values per voxel for each brain region
            fsl_meants = ['fslmeants -i ' source_seq ' -m jhu-tracts_prob_' mask '_bin.nii.gz -o out.txt --showall'];
            system(fsl_meants)
            system('tail -1 out.txt >out1.txt')
            if m == 1
                raw_values_1 = importdata('out1.txt');
            elseif m == 2 
                raw_values_2 = importdata('out1.txt');
            elseif m == 3 
                raw_values_3 = importdata('out1.txt');
            elseif m == 4 
                raw_values_4 = importdata('out1.txt');
            elseif m == 5 
                raw_values_5 = importdata('out1.txt');
            elseif m == 6 
                raw_values_6 = importdata('out1.txt');
            elseif m == 7 
                raw_values_7 = importdata('out1.txt');
            elseif m == 8 
                raw_values_8 = importdata('out1.txt');
            elseif m == 9 
                raw_values_9 = importdata('out1.txt');
            elseif m == 10 
                raw_values_10 = importdata('out1.txt');
            elseif m == 11 
                raw_values_11 = importdata('out1.txt');
            elseif m == 12 
                raw_values_12 = importdata('out1.txt');
            elseif m == 13 
                raw_values_13 = importdata('out1.txt');
            elseif m == 14 
                raw_values_14 = importdata('out1.txt');
            elseif m == 15 
                raw_values_15 = importdata('out1.txt');
            elseif m == 16 
                raw_values_16 = importdata('out1.txt');
            elseif m == 17 
                raw_values_17 = importdata('out1.txt');
            elseif m == 18 
                raw_values_18 = importdata('out1.txt');
            elseif m == 19 
                raw_values_19 = importdata('out1.txt');
            elseif m == 20 
                raw_values_20 = importdata('out1.txt');
            end
        end

        Trv_ndi = num2cell(padcat(raw_values_1,raw_values_2,raw_values_3,raw_values_4, ...
                  raw_values_5,raw_values_6,raw_values_7,raw_values_8,raw_values_9,raw_values_10,raw_values_11,raw_values_12,raw_values_13,...
                  raw_values_14,raw_values_15,raw_values_16,raw_values_17,raw_values_18,raw_values_19,raw_values_20));
        %Truncate names to max length
        mask_list_trunc = regexp(mask_list, '^.{1,63}','match','once');

        Trv_ndi = cell2table(Trv_ndi',"VariableNames", mask_list_trunc);
    
    end

    if n == 2
        for m = 1:length(mask_list)
            mask = mask_list{m};
            fslstats = ['fslstats ' source_seq ' -k jhu-tracts_prob_' mask '_bin.nii.gz -V -M -S >data.txt'];
            system(fslstats)
            if m == 1
                intensities = num2cell(importdata('data.txt'));
                Tfwf = cell2table (intensities, 'VariableNames', {'Voxels_fwf','Volumes_fwf','FWF Mean','SD_fwf'});
            else 
                intensities = num2cell(importdata('data.txt'));
                T2 = cell2table (intensities, 'VariableNames', {'Voxels_fwf','Volumes_fwf','FWF Mean','SD_fwf'});
                Tfwf = [Tfwf;T2];
            end

            fsl_meants = ['fslmeants -i ' source_seq ' -m jhu-tracts_prob_' mask '_bin.nii.gz -o out.txt --showall'];
            system(fsl_meants)
            system('tail -1 out.txt >out1.txt')
            if m == 1
                raw_values_1 = importdata('out1.txt');
            elseif m == 2 
                raw_values_2 = importdata('out1.txt');
            elseif m == 3 
                raw_values_3 = importdata('out1.txt');
            elseif m == 4 
                raw_values_4 = importdata('out1.txt');
            elseif m == 5 
                raw_values_5 = importdata('out1.txt');
            elseif m == 6 
                raw_values_6 = importdata('out1.txt');
            elseif m == 7 
                raw_values_7 = importdata('out1.txt');
            elseif m == 8 
                raw_values_8 = importdata('out1.txt');
            elseif m == 9 
                raw_values_9 = importdata('out1.txt');
            elseif m == 10 
                raw_values_10 = importdata('out1.txt');
            elseif m == 11 
                raw_values_11 = importdata('out1.txt');
            elseif m == 12 
                raw_values_12 = importdata('out1.txt');
            elseif m == 13 
                raw_values_13 = importdata('out1.txt');
            elseif m == 14 
                raw_values_14 = importdata('out1.txt');
            elseif m == 15 
                raw_values_15 = importdata('out1.txt');
            elseif m == 16 
                raw_values_16 = importdata('out1.txt');
            elseif m == 17 
                raw_values_17 = importdata('out1.txt');
            elseif m == 18 
                raw_values_18 = importdata('out1.txt');
            elseif m == 19 
                raw_values_19 = importdata('out1.txt');
            elseif m == 20 
                raw_values_20 = importdata('out1.txt');
            end
        end

        Trv_fwf = num2cell(padcat(raw_values_1,raw_values_2,raw_values_3,raw_values_4, ...
                  raw_values_5,raw_values_6,raw_values_7,raw_values_8,raw_values_9,raw_values_10,raw_values_11,raw_values_12,raw_values_13,...
                  raw_values_14,raw_values_15,raw_values_16,raw_values_17,raw_values_18,raw_values_19,raw_values_20));
        
        %Truncate names to max length
        mask_list_trunc = regexp(mask_list, '^.{1,63}','match','once');
        
        Trv_fwf = cell2table(Trv_fwf',"VariableNames", mask_list_trunc);
    end

    if n == 3
        for m = 1:length(mask_list)
            mask = mask_list{m};
            fslstats = ['fslstats ' source_seq ' -k jhu-tracts_prob_' mask '_bin.nii.gz -V -M -S >data.txt'];
            system(fslstats)

            if m == 1
                intensities = num2cell(importdata('data.txt'));
                Todi = cell2table (intensities, 'VariableNames', {'Voxels_odi','Volumes_odi','ODI Mean','SD_odi'});
            else 
                intensities = num2cell(importdata('data.txt'));
                T2 = cell2table (intensities, 'VariableNames', {'Voxels_odi','Volumes_odi','ODI Mean','SD_odi'});
                Todi = [Todi;T2];
            end

            fsl_meants = ['fslmeants -i ' source_seq ' -m jhu-tracts_prob_' mask '_bin.nii.gz -o out.txt --showall'];
            system(fsl_meants)
            system('tail -1 out.txt >out1.txt')
            if m == 1
                raw_values_1 = importdata('out1.txt');
            elseif m == 2 
                raw_values_2 = importdata('out1.txt');
            elseif m == 3 
                raw_values_3 = importdata('out1.txt');
            elseif m == 4 
                raw_values_4 = importdata('out1.txt');
            elseif m == 5 
                raw_values_5 = importdata('out1.txt');
            elseif m == 6 
                raw_values_6 = importdata('out1.txt');
            elseif m == 7 
                raw_values_7 = importdata('out1.txt');
            elseif m == 8 
                raw_values_8 = importdata('out1.txt');
            elseif m == 9 
                raw_values_9 = importdata('out1.txt');
            elseif m == 10 
                raw_values_10 = importdata('out1.txt');
            elseif m == 11 
                raw_values_11 = importdata('out1.txt');
            elseif m == 12 
                raw_values_12 = importdata('out1.txt');
            elseif m == 13 
                raw_values_13 = importdata('out1.txt');
            elseif m == 14 
                raw_values_14 = importdata('out1.txt');
            elseif m == 15 
                raw_values_15 = importdata('out1.txt');
            elseif m == 16 
                raw_values_16 = importdata('out1.txt');
            elseif m == 17 
                raw_values_17 = importdata('out1.txt');
            elseif m == 18 
                raw_values_18 = importdata('out1.txt');
            elseif m == 19 
                raw_values_19 = importdata('out1.txt');
            elseif m == 20 
                raw_values_20 = importdata('out1.txt');
            end
        end

        Trv_odi = num2cell(padcat(raw_values_1,raw_values_2,raw_values_3,raw_values_4, ...
                  raw_values_5,raw_values_6,raw_values_7,raw_values_8,raw_values_9,raw_values_10,raw_values_11,raw_values_12,raw_values_13,...
                  raw_values_14,raw_values_15,raw_values_16,raw_values_17,raw_values_18,raw_values_19,raw_values_20));
        
        %Truncate names to max length
        mask_list_trunc = regexp(mask_list, '^.{1,63}','match','once');
       
        Trv_odi = cell2table(Trv_odi',"VariableNames",mask_list_trunc);
    end
end

mask_list_key = mask_list';

key = cell2table(mask_list_key, 'VariableNames',{'Region'});

Tndi = [key, Tndi];
Todi = [key, Todi];
Tfwf = [key, Tfwf];

Tjoin = join(Tndi, Todi);
T_combined = join(Tjoin, Tfwf);
%%
writetable(T_combined, 'JHU_tract_NODDI_indices.csv')
%%
writetable(Trv_ndi, 'JHU_tract_NDI_raw_intensities.csv')
writetable(Trv_fwf, 'JHU_tract_FWF_raw_intensities.csv')
writetable(Trv_odi, 'JHU_tract_ODI_raw_intensities.csv')

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




%% template codes %%

%% epireg technique for coregistration
% Register WM to Atlas
source_seq = ['anat_seq_brain_mask_pve_2'];
ref_seq = ['/mnt/c/WSL2_dir/Atlases/HarvardOxford-sub-maxprob-thr0-1mm'];
flirt_coreg = ['flirt -in ' source_seq ' -ref ' ref_seq ' -out r_cort_' source_seq ' -omat invol2refvol.mat -v' ];
system(flirt_coreg)

% Register NODDI to T1
source_seq_list = {'Case1_odi','Case1_ficvf','Case1_fiso'};
for n = 1:length(source_seq_list)

    source_seq = source_seq_list{n};

    epi_reg = ['epi_reg --epi=' source_seq ' --t1=' anat_seq ' --t1brain=anat_seq_brain_mask --wmseg=anat_seq_brain_mask_pve_2 --out=t1r_' source_seq ' -v'];
    system(epi_reg)
end

% Apply atlas matrix to NODDI 
source_seq_list = {'t1r_Case1_odi','t1r_Case1_ficvf','t1r_Case1_fiso'};
for n = 1:length(source_seq_list)

    source_seq = source_seq_list{n};
    ref_seq = ['/mnt/c/WSL2_dir/Atlases/HarvardOxford-sub-maxprob-thr0-1mm'];
    flirt_coreg = ['flirt -in ' source_seq ' -ref ' ref_seq ' -out r_cort_' source_seq ' -init invol2refvol.mat -applyxfm -v' ];
    system(flirt_coreg)

end

% Extract WM from NODDI using mask
WM_mask = ['r_cort_anat_seq_brain_mask_pve_2'];
fsl_maths = ['fslmaths ' WM_mask ' -bin ' WM_mask '_bin.nii.gz'];
system(fsl_maths)

ref_seq = ['r_cort_t1r_Case1_ficvf'];
fslmaths = ['fslmaths ' ref_seq ' -mul ' WM_mask '_bin ' ref_seq '_WM'];
system(fslmaths)

ref_seq = ['r_cort_t1r_Case1_fiso'];
fslmaths = ['fslmaths ' ref_seq ' -mul ' WM_mask '_bin ' ref_seq '_WM'];
system(fslmaths)

ref_seq = ['r_cort_t1r_Case1_odi'];
fslmaths = ['fslmaths ' ref_seq ' -mul ' WM_mask '_bin ' ref_seq '_WM'];
system(fslmaths)


