% Making maps of ethanol NODDI differences with one control scan at once
% Jesse M Klostranec MD PhD FRCPC
% January 17, 2024

clc
clear

disp(' Ethanol-NODDI DTI Program ')
disp(' ')
disp(' January 17, 2024 ')
disp(' By Jesse M Klostranec MD PhD FRCPC ')
disp(' ')
disp(' Press any key to continue ...')
pause

disp('---------')
disp(' Ethanol NODDI Analysis Using the Harvard Oxford Atlas ')
disp('---------')
disp(' ')

% Load the NODDI data in NIFTI format for averaging the three control runs
% pre-EtOH including only the ROI defined by the post-EtOH nodif mask

mask = niftiread('rnodif_brain_mask.nii');

% Ask whether the first, second, or third control acquisitions will be used

control_num = input('Which control scan would you like to compare with? (1, 2, or 3): ');

NDI_name_pre = strcat('rPre-EtOH_',num2str(control_num),'_ficvf.nii');
ODI_name_pre = strcat('rPre-EtOH_',num2str(control_num),'_odi.nii');
FWF_name_pre = strcat('rPre-EtOH_',num2str(control_num),'_fiso.nii');

% Import data for the control run

    % NDI_pre
place_holder = niftiread(NDI_name_pre);
image_info = niftiinfo(NDI_name_pre);
Sizer = size(place_holder);
Pre_NDI = double(zeros(Sizer));
for i = 1: Sizer(1,1)
    for j = 1:Sizer(1,2)
        for k=1:Sizer(1,3)
            if mask(i,j,k) == 0
                Pre_NDI(i,j,k) = 0;
            else
                Pre_NDI(i,j,k) = place_holder(i,j,k);
            end
        end
    end
end
Max_ph = max(Pre_NDI, [], "all");
Min_ph = min(Pre_NDI, [], "all");
for i = 1: Sizer(1,1)
    for j = 1:Sizer(1,2)
        for k=1:Sizer(1,3)
            if mask(i,j,k) == 0
                Pre_NDI(i,j,k) = 0;
            else
                Pre_NDI(i,j,k) = (Pre_NDI(i,j,k) - Min_ph)/(Max_ph - Min_ph);
            end
        end
    end
end
Sizer = [];
place_holder = [];

    % ODI_pre
place_holder = niftiread(ODI_name_pre);
Sizer = size(place_holder);
Pre_ODI = double(zeros(Sizer));
for i = 1: Sizer(1,1)
    for j = 1:Sizer(1,2)
        for k=1:Sizer(1,3)
            if mask(i,j,k) == 0
                Pre_ODI(i,j,k) = 0;
            else
                Pre_ODI(i,j,k) = place_holder(i,j,k);
            end
        end
    end
end
Max_ph = max(Pre_ODI, [], "all");
Min_ph = min(Pre_ODI, [], "all");
for i = 1: Sizer(1,1)
    for j = 1:Sizer(1,2)
        for k=1:Sizer(1,3)
            if mask(i,j,k) == 0
                Pre_ODI(i,j,k) = 0;
            else
                Pre_ODI(i,j,k) = (Pre_ODI(i,j,k) - Min_ph)/(Max_ph - Min_ph);
            end
        end
    end
end
Sizer = [];
place_holder = [];

    % FWF_pre
place_holder = niftiread(FWF_name_pre);
Sizer = size(place_holder);
Pre_FWF = double(zeros(Sizer));
for i = 1: Sizer(1,1)
    for j = 1:Sizer(1,2)
        for k=1:Sizer(1,3)
            if mask(i,j,k) == 0
                Pre_FWF(i,j,k) = 0;
            else
                Pre_FWF(i,j,k) = place_holder(i,j,k);
            end
        end
    end
end
Max_ph = max(Pre_FWF, [], "all");
Min_ph = min(Pre_FWF, [], "all");
for i = 1: Sizer(1,1)
    for j = 1:Sizer(1,2)
        for k=1:Sizer(1,3)
            if mask(i,j,k) == 0
                Pre_FWF(i,j,k) = 0;
            else
                Pre_FWF(i,j,k) = (Pre_FWF(i,j,k) - Min_ph)/(Max_ph - Min_ph);
            end
        end
    end
end
Sizer = [];
place_holder = [];

    % Pre-Ethanol Diffusion tensor import
place_holder = niftiread(strcat('rdti_tensor_',num2str(control_num),'.nii');
Sizer = size(place_holder);
Pre_Dxx = double(zeros(Sizer(1,1),Sizer(1,2),Sizer(1,3)));
Pre_Dyy = double(zeros(Sizer(1,1),Sizer(1,2),Sizer(1,3)));
Pre_Dzz = double(zeros(Sizer(1,1),Sizer(1,2),Sizer(1,3)));
for i = 1: Sizer(1,1)
    for j = 1:Sizer(1,2)
        for k=1:Sizer(1,3)
            if mask(i,j,k) == 0
                Pre_Dxx(i,j,k) = 0;
                Pre_Dyy(i,j,k) = 0;
                Pre_Dzz(i,j,k) = 0;
            else
                Pre_Dxx(i,j,k) = place_holder(i,j,k,1);
                Pre_Dyy(i,j,k) = place_holder(i,j,k,4);
                Pre_Dzz(i,j,k) = place_holder(i,j,k,6);
            end
        end
    end
end
Sizer = [];
place_holder = [];   

    % Post-Ethanol Diffusion tensor import
place_holder = niftiread(strcat('rdti_tensor.nii');
Sizer = size(place_holder);
Post_Dxx = double(zeros(Sizer(1,1),Sizer(1,2),Sizer(1,3)));
Post_Dyy = double(zeros(Sizer(1,1),Sizer(1,2),Sizer(1,3)));
Post_Dzz = double(zeros(Sizer(1,1),Sizer(1,2),Sizer(1,3)));
for i = 1: Sizer(1,1)
    for j = 1:Sizer(1,2)
        for k=1:Sizer(1,3)
            if mask(i,j,k) == 0
                Post_Dxx(i,j,k) = 0;
                Post_Dyy(i,j,k) = 0;
                Post_Dzz(i,j,k) = 0;
            else
                Post_Dxx(i,j,k) = place_holder(i,j,k,1);
                Post_Dyy(i,j,k) = place_holder(i,j,k,4);
                Post_Dzz(i,j,k) = place_holder(i,j,k,6);
            end
        end
    end
end
Sizer = [];
place_holder = [];  

% Remove all voxels for which the FWF is greater that 35%

Sizer = size(Pre_NDI);
for i = 1: Sizer(1,1)
    for j = 1:Sizer(1,2)
        for k=1:Sizer(1,3)
            if Pre_FWF(i,j,k) >= 0.35
                Pre_NDI(i,j,k) = 0;
                Pre_ODI(i,j,k) = 0;
                Pre_Dxx(i,j,k) = 0;
                Pre_Dyy(i,j,k) = 0;
                Pre_Dzz(i,j,k) = 0;
                Post_Dxx(i,j,k) = 0;
                Post_Dyy(i,j,k) = 0;
                Post_Dzz(i,j,k) = 0;
            else
                Pre_NDI(i,j,k) = Pre_NDI(i,j,k);
                Pre_ODI(i,j,k) = Pre_ODI(i,j,k);
                Pre_Dxx(i,j,k) = Pre_Dxx(i,j,k);
                Pre_Dyy(i,j,k) = Pre_Dyy(i,j,k);
                Pre_Dzz(i,j,k) = Pre_Dzz(i,j,k);
                Post_Dxx(i,j,k) = Post_Dxx(i,j,k);
                Post_Dyy(i,j,k) = Post_Dyy(i,j,k);
                Post_Dzz(i,j,k) = Post_Dzz(i,j,k);
            end
        end
    end
end
Sizer = [];

    % Apply isotropic Gaussian smoothing kernel of 2 standard deviations
 
 sPre_NDI = imgaussfilt(Pre_NDI,1);
 sPre_ODI = imgaussfilt(Pre_ODI,1);
 sPre_FWF = imgaussfilt(Pre_FWF,1);

% Load the Post-EtOH data

    % NDI_post
place_holder = niftiread('rPost-EtOH_ficvf.nii');
info_NDI = niftiinfo('rPost-EtOH_ficvf.nii');
Sizer = size(place_holder);
Post_NDI = double(zeros(Sizer));
for i = 1: Sizer(1,1)
    for j = 1:Sizer(1,2)
        for k=1:Sizer(1,3)
            if mask(i,j,k) == 0
                Post_NDI(i,j,k) = 0;
            else
                Post_NDI(i,j,k) = place_holder(i,j,k);
            end
        end
    end
end
Max_ph = max(Post_NDI, [], "all");
Min_ph = min(Post_NDI, [], "all");
for i = 1: Sizer(1,1)
    for j = 1:Sizer(1,2)
        for k=1:Sizer(1,3)
            if mask(i,j,k) == 0
                Post_NDI(i,j,k) = 0;
            else
                Post_NDI(i,j,k) = (Post_NDI(i,j,k) - Min_ph)/(Max_ph - Min_ph);
            end
        end
    end
end
Sizer = [];
place_holder = [];

    % ODI_post
place_holder = niftiread('rPost-EtOH_odi.nii');
info_ODI = niftiinfo('rPost-EtOH_odi.nii');
Sizer = size(place_holder);
Post_ODI = double(zeros(Sizer));
for i = 1: Sizer(1,1)
    for j = 1:Sizer(1,2)
        for k=1:Sizer(1,3)
            if mask(i,j,k) == 0
                Post_ODI(i,j,k) = 0;
            else
                Post_ODI(i,j,k) = place_holder(i,j,k);
            end
        end
    end
end
Max_ph = max(Post_ODI, [], "all");
Min_ph = min(Post_ODI, [], "all");
for i = 1: Sizer(1,1)
    for j = 1:Sizer(1,2)
        for k=1:Sizer(1,3)
            if mask(i,j,k) == 0
                Post_ODI(i,j,k) = 0;
            else
                Post_ODI(i,j,k) = (Post_ODI(i,j,k) - Min_ph)/(Max_ph - Min_ph);
            end
        end
    end
end
Sizer = [];
place_holder = [];

    % FWF_post
place_holder = niftiread('rPost-EtOH_fiso.nii');
info_FWF = niftiinfo('rPost-EtOH_fiso.nii');
Sizer = size(place_holder);
Post_FWF = double(zeros(Sizer));
for i = 1: Sizer(1,1)
    for j = 1:Sizer(1,2)
        for k=1:Sizer(1,3)
            if mask(i,j,k) == 0
                Post_FWF(i,j,k) = 0;
            else
                Post_FWF(i,j,k) = place_holder(i,j,k);
            end
        end
    end
end
Max_ph = max(Post_FWF, [], "all");
Min_ph = min(Post_FWF, [], "all");
for i = 1: Sizer(1,1)
    for j = 1:Sizer(1,2)
        for k=1:Sizer(1,3)
            if mask(i,j,k) == 0
                Post_FWF(i,j,k) = 0;
            else
                Post_FWF(i,j,k) = (Post_FWF(i,j,k) - Min_ph)/(Max_ph - Min_ph);
            end
        end
    end
end
Sizer = [];
place_holder = [];

    % Use the Gaussian smoothing kernel

sPost_NDI = imgaussfilt(Post_NDI,1);
sPost_ODI = imgaussfilt(Post_ODI,1);
sPost_FWF = imgaussfilt(Post_FWF,1);

% Load the Harvard-Oxford atlases

HOc = double(niftiread('HarvardOxford-cort-maxprob-thr0-1mm.nii'));
sizeHO = size(HOc);
HOsc = double(niftiread('HarvardOxford-sub-maxprob-thr0-1mm.nii'));

% Load the MNI atlas

MNIt = double(niftiread('MNI-maxprob-thr0-1mm.nii'));

% Load the segmented grey matter image

GM = double(niftiread('T1_masked_gm.nii'));

% Match the orientation and sizes of the two images:

% These are the following coordinates for the different scans:
%   Max: voxel [0 0 0] = (-98 -134 -72)
%              [196 232 188] = (98.5 98.5 116.5)
%
%   MNI: voxel [0 0 0] = (90 -126 -72)
%              [181 217 181] = (-91.5 91.5 109.5)

new_HOc = zeros(size(sPost_NDI));
sizer = size(new_HOc);
new_HOsc = zeros(size(sPost_NDI));
new_MNI = zeros(size(sPost_NDI));

for i = 9:(9 - 1 + sizeHO(1,1))
    for j = 9:(9 - 1 + sizeHO(1,2))
        for k = 1:sizeHO(1,3)
            new_HOc(i,j,k) = HOc((i-8),(j-8),k);
            new_HOsc(i,j,k) = HOsc((i-8),(j-8),k);
            new_MNI(i,j,k) = MNIt((i-8),(j-8),k);
        end
    end
end
sizer = [];

% For the MNI and Harvard Oxford atlases, the following values correspond to the following
% anatomic regions:
%
%   1 - caudates
%   2 - cerebellar hemispheres
%   3 - frontal lobes
%   4 - insula
%   5 - occipital lobes
%   6 - parietal lobes
%   7 - lentiform nuclei
%   8 - temporal lobes
%   9 - thalami
%

% For each atlas, there will be a different number of ROIs, so the matrices
% of collected values will be different, ie. for MNI - 9 ROIs, for Harvard
% Oxford cortical - 48 ROIs, for Harvard Oxford subcortical - 21 ROIs

% Let's start with the MNI atlas:

    % Vector of values selected from a single MNI mask
NDI_pre = [];
ODI_pre = [];
FWF_pre = [];
NDI_post = [];
ODI_post = [];
FWF_post = [];

    % Matrices to store the values selected from all MNI masks
NDI_MNI_pre = [];
ODI_MNI_pre = [];
FWF_MNI_pre = [];
NDI_MNI_post = [];
ODI_MNI_post = [];
FWF_MNI_post = [];

    % Vectors for mean values
Pre_mean_NDI = zeros(1,9);
Pre_mean_ODI = zeros(1,9);
Pre_mean_FWF = zeros(1,9);
Post_mean_NDI = zeros(1,9);
Post_mean_ODI = zeros(1,9);
Post_mean_FWF = zeros(1,9);

    % Vectors for std values
Pre_std_NDI = zeros(1,9);
Pre_std_ODI = zeros(1,9);
Pre_std_FWF = zeros(1,9);
Post_std_NDI = zeros(1,9);
Post_std_ODI = zeros(1,9);
Post_std_FWF = zeros(1,9);

% Counts the number of voxels in each ROI
MNI_NDI_numbers = [];

sizer = size(sPre_NDI);

for z = 1:9
    for i = 1:sizer(1,1)
        for j = 1:sizer(1,2)
            for k = 1:sizer(1,3)
                if new_MNI(i,j,k) == z
                    if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                    elseif GM(i,j,k) == 0

                    else
                        NDI_pre = [NDI_pre sPre_NDI(i,j,k)];
                        NDI_post = [NDI_post sPost_NDI(i,j,k)];
                    end
                end
            end
        end
    end
    
    a = length(NDI_MNI_pre);
    b = length(NDI_pre);

    MNI_NDI_numbers = [MNI_NDI_numbers; b];
    
    Pre_mean_NDI(1,z) = mean(NDI_pre);
    Post_mean_NDI(1,z) = mean(NDI_post);

    Pre_std_NDI(1,z) = std(NDI_pre);
    Post_std_NDI(1,z) = std(NDI_post);

    if a > b

        NDI_MNI_pre = [NDI_MNI_pre; zeros(1,a)];
        NDI_MNI_post = [NDI_MNI_post; zeros(1,a)];

        for i = 1:b

            NDI_MNI_pre(z,i) = NDI_pre(1,i);
            NDI_MNI_post(z,i) = NDI_post(1,i);

        end

    else

        if z == 1

            NDI_MNI_pre(1,:) = NDI_pre(1,:);
            NDI_MNI_post(1,:) = NDI_post(1,:);

        else

            NDI_MNI_ph1 = zeros(z,b);
            NDI_MNI_ph2 = zeros(z,b);

            for i = 1:(z-1)
                for j = 1:a
                    NDI_MNI_ph1(i,j) = NDI_MNI_pre(i,j);
                    NDI_MNI_ph2(i,j) = NDI_MNI_post(i,j);
                end
            end

            for i = 1:a

                NDI_MNI_ph1(z,i) = NDI_pre(1,i);
                NDI_MNI_ph2(z,i) = NDI_post(1,i);

            end

            NDI_MNI_pre = NDI_MNI_ph1;
            NDI_MNI_post = NDI_MNI_ph2;

        end
end

    NDI_pre = [];
    NDI_post = [];

end


% Counts the number of voxels in each ROI
MNI_ODI_numbers = [];

sizer = size(sPre_NDI);

for z = 1:9
    for i = 1:sizer(1,1)
        for j = 1:sizer(1,2)
            for k = 1:sizer(1,3)
                if new_MNI(i,j,k) == z
                    if sPre_ODI(i,j,k) == 0 || sPost_ODI(i,j,k) == 0 || isnan(sPre_ODI(i,j,k)) == 1 || isnan(sPost_ODI(i,j,k)) == 1

                    elseif GM(i,j,k) == 0

                    else
                        ODI_pre = [ODI_pre sPre_ODI(i,j,k)];
                        ODI_post = [ODI_post sPost_ODI(i,j,k)];
                    end
                end
            end
        end
    end
    
    a = length(ODI_MNI_pre);
    b = length(ODI_pre);

    MNI_ODI_numbers = [MNI_ODI_numbers; b];
    
    Pre_mean_ODI(1,z) = mean(ODI_pre);
    Post_mean_ODI(1,z) = mean(ODI_post);

    Pre_std_ODI(1,z) = std(ODI_pre);
    Post_std_ODI(1,z) = std(ODI_post);

    if a > b

        ODI_MNI_pre = [ODI_MNI_pre; zeros(1,a)];
        ODI_MNI_post = [ODI_MNI_post; zeros(1,a)];

        for i = 1:b

            ODI_MNI_pre(z,i) = ODI_pre(1,i);
            ODI_MNI_post(z,i) = ODI_post(1,i);

        end

    else

        if z == 1

            ODI_MNI_pre(1,:) = ODI_pre(1,:);
            ODI_MNI_post(1,:) = ODI_post(1,:);

        else

            ODI_MNI_ph1 = zeros(z,b);
            ODI_MNI_ph2 = zeros(z,b);

            for i = 1:(z-1)
                for j = 1:a
                    ODI_MNI_ph1(i,j) = ODI_MNI_pre(i,j);
                    ODI_MNI_ph2(i,j) = ODI_MNI_post(i,j);
                end
            end

            for i = 1:a

                ODI_MNI_ph1(z,i) = ODI_pre(1,i);
                ODI_MNI_ph2(z,i) = ODI_post(1,i);

            end

            ODI_MNI_pre = ODI_MNI_ph1;
            ODI_MNI_post = ODI_MNI_ph2;

        end
end

    ODI_pre = [];
    ODI_post = [];

end

% Counts the number of voxels in each ROI
MNI_FWF_numbers = [];

sizer = size(sPre_NDI);

for z = 1:9
    for i = 1:sizer(1,1)
        for j = 1:sizer(1,2)
            for k = 1:sizer(1,3)
                if new_MNI(i,j,k) == z
                    if sPre_FWF(i,j,k) == 0 || sPost_FWF(i,j,k) == 0 || isnan(sPre_FWF(i,j,k)) == 1 || isnan(sPost_FWF(i,j,k)) == 1

                    elseif GM(i,j,k) == 0

                    else
                        FWF_pre = [FWF_pre sPre_FWF(i,j,k)];
                        FWF_post = [FWF_post sPost_FWF(i,j,k)];
                    end
                end
            end
        end
    end
    
    a = length(FWF_MNI_pre);
    b = length(FWF_pre);

    MNI_FWF_numbers = [MNI_FWF_numbers; b];
    
    Pre_mean_ECF(1,z) = (1-mean(FWF_pre))*(1-Pre_mean_NDI(1,z));
    Post_mean_ECF(1,z) = (1-mean(FWF_post))*(1-Post_mean_NDI(1,z));

    Pre_std_FWF(1,z) = std(FWF_pre);
    Post_std_FWF(1,z) = std(FWF_post);

    if a > b

        FWF_MNI_pre = [FWF_MNI_pre; zeros(1,a)];
        FWF_MNI_post = [FWF_MNI_post; zeros(1,a)];

        for i = 1:b

            FWF_MNI_pre(z,i) = FWF_pre(1,i);
            FWF_MNI_post(z,i) = FWF_post(1,i);

        end

    else

        if z == 1

            FWF_MNI_pre(1,:) = FWF_pre(1,:);
            FWF_MNI_post(1,:) = FWF_post(1,:);

        else

            FWF_MNI_ph1 = zeros(z,b);
            FWF_MNI_ph2 = zeros(z,b);

            for i = 1:(z-1)
                for j = 1:a
                    FWF_MNI_ph1(i,j) = FWF_MNI_pre(i,j);
                    FWF_MNI_ph2(i,j) = FWF_MNI_post(i,j);
                end
            end

            for i = 1:a

                FWF_MNI_ph1(z,i) = FWF_pre(1,i);
                FWF_MNI_ph2(z,i) = FWF_post(1,i);

            end

            FWF_MNI_pre = FWF_MNI_ph1;
            FWF_MNI_post = FWF_MNI_ph2;

        end
end

    FWF_pre = [];
    FWF_post = [];

end

% Calculate the statistics

fileID = fopen(strcat('ROI_MNI_results_',num2str(control_num),'.txt'),'w');
fprintf(fileID,'%16s %12s %12s %12s %12s %12s %12s %12s\n','ROI','Index','Pre_Mean','Pre_Std','Post_Mean','Post_Std','Difference','p-value');

MNI_NDI_image = zeros(sizer(1,1),sizer(1,2),sizer(1,3));
MNI_ODI_image = zeros(sizer(1,1),sizer(1,2),sizer(1,3));
MNI_ECF_image = zeros(sizer(1,1),sizer(1,2),sizer(1,3));

for z = 1:9

    Delta_NDI = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
    Delta_ODI = (Pre_mean_ODI(1,z) - Post_mean_ODI(1,z));
    Delta_ECF = (Pre_mean_ECF(1,z) - Post_mean_ECF(1,z));

    values_pre = [];
    values_post = [];
    Ovalues_pre = [];
    Ovalues_post = [];

    for i = 1:MNI_NDI_numbers(z,1)
        values_pre = [values_pre NDI_MNI_pre(z,i)];
        values_post = [values_post NDI_MNI_post(z,i)];
        Ovalues_pre = [Ovalues_pre ODI_MNI_pre(z,i)];
        Ovalues_post = [Ovalues_post ODI_MNI_post(z,i)];
    end
    
    [hN1,pN1] = ttest2(values_pre,values_post);
    [hO1,pO1] = ttest2(Ovalues_pre,Ovalues_post);

    if z == 1
        ROI = strcat('Caudates');
        
        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_MNI(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            MNI_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                            MNI_ECF_image(i,j,k) = (Pre_mean_ECF(1,z) - Post_mean_ECF(1,z));
                        end
                    end
                end
            end
        end  

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_MNI(i,j,k) == z
                        if sPre_ODI(i,j,k) == 0 || sPost_ODI(i,j,k) == 0 || isnan(sPre_ODI(i,j,k)) == 1 || isnan(sPost_ODI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            MNI_ODI_image(i,j,k) = (Pre_mean_ODI(1,z) - Post_mean_ODI(1,z));
                        end
                    end
                end
            end
        end  

    elseif z == 2
        ROI = strcat('Cerebellum');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_MNI(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            MNI_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                            MNI_ECF_image(i,j,k) = (Pre_mean_ECF(1,z) - Post_mean_ECF(1,z));
                        end
                    end
                end
            end
        end  

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_MNI(i,j,k) == z
                        if sPre_ODI(i,j,k) == 0 || sPost_ODI(i,j,k) == 0 || isnan(sPre_ODI(i,j,k)) == 1 || isnan(sPost_ODI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            MNI_ODI_image(i,j,k) = (Pre_mean_ODI(1,z) - Post_mean_ODI(1,z));
                        end
                    end
                end
            end
        end  

    elseif z == 3
        ROI = strcat('Frontal Lobes');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_MNI(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            MNI_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                            MNI_ECF_image(i,j,k) = (Pre_mean_ECF(1,z) - Post_mean_ECF(1,z));
                        end
                    end
                end
            end
        end  

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_MNI(i,j,k) == z
                        if sPre_ODI(i,j,k) == 0 || sPost_ODI(i,j,k) == 0 || isnan(sPre_ODI(i,j,k)) == 1 || isnan(sPost_ODI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            MNI_ODI_image(i,j,k) = (Pre_mean_ODI(1,z) - Post_mean_ODI(1,z));
                        end
                    end
                end
            end
        end  

    elseif z == 4
        ROI = strcat('Insula');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_MNI(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            MNI_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                            MNI_ECF_image(i,j,k) = (Pre_mean_ECF(1,z) - Post_mean_ECF(1,z));
                        end
                    end
                end
            end
        end  

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_MNI(i,j,k) == z
                        if sPre_ODI(i,j,k) == 0 || sPost_ODI(i,j,k) == 0 || isnan(sPre_ODI(i,j,k)) == 1 || isnan(sPost_ODI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            MNI_ODI_image(i,j,k) = (Pre_mean_ODI(1,z) - Post_mean_ODI(1,z));
                        end
                    end
                end
            end
        end  

    elseif z == 5
        ROI = strcat('Occipital Lobes');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_MNI(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            MNI_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                            MNI_ECF_image(i,j,k) = (Pre_mean_ECF(1,z) - Post_mean_ECF(1,z));
                        end
                    end
                end
            end
        end 

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_MNI(i,j,k) == z
                        if sPre_ODI(i,j,k) == 0 || sPost_ODI(i,j,k) == 0 || isnan(sPre_ODI(i,j,k)) == 1 || isnan(sPost_ODI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            MNI_ODI_image(i,j,k) = (Pre_mean_ODI(1,z) - Post_mean_ODI(1,z));
                        end
                    end
                end
            end
        end  

    elseif z == 6
        ROI = strcat('Parietal Lobes');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_MNI(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            MNI_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                            MNI_ECF_image(i,j,k) = (Pre_mean_ECF(1,z) - Post_mean_ECF(1,z));
                        end
                    end
                end
            end
        end 

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_MNI(i,j,k) == z
                        if sPre_ODI(i,j,k) == 0 || sPost_ODI(i,j,k) == 0 || isnan(sPre_ODI(i,j,k)) == 1 || isnan(sPost_ODI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            MNI_ODI_image(i,j,k) = (Pre_mean_ODI(1,z) - Post_mean_ODI(1,z));
                        end
                    end
                end
            end
        end  

    elseif z == 7 
        ROI = strcat('Lentiform Nuclei');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_MNI(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            MNI_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                            MNI_ECF_image(i,j,k) = (Pre_mean_ECF(1,z) - Post_mean_ECF(1,z));
                        end
                    end
                end
            end
        end 

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_MNI(i,j,k) == z
                        if sPre_ODI(i,j,k) == 0 || sPost_ODI(i,j,k) == 0 || isnan(sPre_ODI(i,j,k)) == 1 || isnan(sPost_ODI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            MNI_ODI_image(i,j,k) = (Pre_mean_ODI(1,z) - Post_mean_ODI(1,z));
                        end
                    end
                end
            end
        end  

    elseif z == 8
        ROI = strcat('Temporal Lobes');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_MNI(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            MNI_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                            MNI_ECF_image(i,j,k) = (Pre_mean_ECF(1,z) - Post_mean_ECF(1,z));
                        end
                    end
                end
            end
        end 

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_MNI(i,j,k) == z
                        if sPre_ODI(i,j,k) == 0 || sPost_ODI(i,j,k) == 0 || isnan(sPre_ODI(i,j,k)) == 1 || isnan(sPost_ODI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            MNI_ODI_image(i,j,k) = (Pre_mean_ODI(1,z) - Post_mean_ODI(1,z));
                        end
                    end
                end
            end
        end  

    else
        ROI = strcat('Thalami');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_MNI(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            MNI_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                            MNI_ECF_image(i,j,k) = (Pre_mean_ECF(1,z) - Post_mean_ECF(1,z));
                        end
                    end
                end
            end
        end 

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_MNI(i,j,k) == z
                        if sPre_ODI(i,j,k) == 0 || sPost_ODI(i,j,k) == 0 || isnan(sPre_ODI(i,j,k)) == 1 || isnan(sPost_ODI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            MNI_ODI_image(i,j,k) = (Pre_mean_ODI(1,z) - Post_mean_ODI(1,z));
                        end
                    end
                end
            end
        end  


    end

fprintf(fileID,'%16s %12s %5s %1.4f %5s %1.4f %5s %1.4f %5s %1.4f %4s %1.4f %5s %1.4f\n',ROI,'NDI',' ',Pre_mean_NDI(1,z),' ',Pre_std_NDI(1,z),' ',Post_mean_NDI(1,z),' ',Post_std_NDI(1,z),' ',Delta_NDI,' ',pN1);
fprintf(fileID,'%16s %12s %5s %1.4f %5s %1.4f %5s %1.4f %5s %1.4f %4s %1.4f %5s %1.4f\n',' ','ODI',' ',Pre_mean_ODI(1,z),' ',Pre_std_ODI(1,z),' ',Post_mean_ODI(1,z),' ',Post_std_ODI(1,z),' ',Delta_ODI,' ',pO1);
fprintf(fileID,'%16s %12s %5s %1.4f %5s %1.4f %5s %1.4f %5s %1.4f %4s %1.4f %5s %1.4f\n',' ','ECF',' ',Pre_mean_ECF(1,z),' ',' ',' ',Post_mean_ECF(1,z),' ',' ',' ',Delta_ECF,' ',' ');

end

fclose(fileID);

niftiwrite(MNI_NDI_image,'MNI_NDI_image',image_info)
niftiwrite(MNI_ODI_image,'MNI_ODI_image',image_info)
niftiwrite(MNI_ECF_image,'MNI_ECF_image',image_info)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Now, let's do the Harvard Oxford Cortical Brain Atlas

    % Vector of values selected from a single HOc mask
NDI_pre = [];
ODI_pre = [];
FWF_pre = [];
NDI_post = [];
ODI_post = [];
FWF_post = [];

    % Matrices to store the values selected from all HOc masks
NDI_HOc_pre = [];
ODI_HOc_pre = [];
FWF_HOc_pre = [];
NDI_HOc_post = [];
ODI_HOc_post = [];
FWF_HOc_post = [];

    % Vectors for mean values
Pre_mean_NDI = zeros(1,48);
Pre_mean_ODI = zeros(1,48);
Pre_mean_FWF = zeros(1,48);
Post_mean_NDI = zeros(1,48);
Post_mean_ODI = zeros(1,48);
Post_mean_FWF = zeros(1,48);

    % Vectors for std values
Pre_std_NDI = zeros(1,48);
Pre_std_ODI = zeros(1,48);
Pre_std_FWF = zeros(1,48);
Post_std_NDI = zeros(1,48);
Post_std_ODI = zeros(1,48);
Post_std_FWF = zeros(1,48);

% Counts the number of voxels in each ROI
HOc_NDI_numbers = [];

sizer = size(sPre_NDI);

for z = 1:48
    for i = 1:sizer(1,1)
        for j = 1:sizer(1,2)
            for k = 1:sizer(1,3)
                if new_HOc(i,j,k) == z
                    if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                    elseif GM(i,j,k) == 0

                    else
                        NDI_pre = [NDI_pre sPre_NDI(i,j,k)];
                        NDI_post = [NDI_post sPost_NDI(i,j,k)];
                    end
                end
            end
        end
    end
    
    a = length(NDI_HOc_pre);
    b = length(NDI_pre);

    HOc_NDI_numbers = [HOc_NDI_numbers; b];
    
    Pre_mean_NDI(1,z) = mean(NDI_pre);
    Post_mean_NDI(1,z) = mean(NDI_post);

    Pre_std_NDI(1,z) = std(NDI_pre);
    Post_std_NDI(1,z) = std(NDI_post);

    if a > b

        NDI_HOc_pre = [NDI_HOc_pre; zeros(1,a)];
        NDI_HOc_post = [NDI_HOc_post; zeros(1,a)];

        for i = 1:b

            NDI_HOc_pre(z,i) = NDI_pre(1,i);
            NDI_HOc_post(z,i) = NDI_post(1,i);

        end

    else

        if z == 1

            NDI_HOc_pre(1,:) = NDI_pre(1,:);
            NDI_HOc_post(1,:) = NDI_post(1,:);

        else

        NDI_HOc_ph1 = zeros(z,b);
        NDI_HOc_ph2 = zeros(z,b);

            for i = 1:(z-1)
                for j = 1:a
                    NDI_HOc_ph1(i,j) = NDI_HOc_pre(i,j);
                    NDI_HOc_ph2(i,j) = NDI_HOc_post(i,j);
                end
            end

            for i = 1:a

                NDI_HOc_ph1(z,i) = NDI_pre(1,i);
                NDI_HOc_ph2(z,i) = NDI_post(1,i);

            end

            NDI_HOc_pre = NDI_HOc_ph1;
            NDI_HOc_post = NDI_HOc_ph2;

        end
    end

    NDI_pre = [];
    NDI_post = [];

end

% Calculate the statistics

fileID = fopen(strcat('ROI_HOc_results_',num2str(control_num),'.txt'),'w');
fprintf(fileID,'%47s %12s %12s %12s %12s %12s %12s %12s\n','ROI','Index','Pre_Mean','Pre_Std','Post_Mean','Post_Std','Difference','p-value');

HOc_NDI_image = zeros(sizer(1,1),sizer(1,2),sizer(1,3));

for z = 1:48

    Delta_NDI = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));

    values_pre = [];
    values_post = [];

    for i = 1:HOc_NDI_numbers(z,1)
        values_pre = [values_pre NDI_HOc_pre(z,i)];
        values_post = [values_post NDI_HOc_post(z,i)];
    end
    
    [hN1,pN1] = ttest2(values_pre,values_post);

    if z == 1
        ROI = strcat('Frontal Pole');
        
        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end        

    elseif z == 2
        ROI = strcat('Insular Cortex');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end  

    elseif z == 3
        ROI = strcat('Superior Frontal Gyrus');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end  

    elseif z == 4
        ROI = strcat('Middle Frontal gyrus');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end  

    elseif z == 5
        ROI = strcat('Inferior Frontal Gyrus, Pars Triangularis');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 6
        ROI = strcat('Inferior Frontal Gyrus, Pars Opercularis');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 7 
        ROI = strcat('Precentral Gyrus');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 8
        ROI = strcat('Temporal Pole');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 9
        ROI = strcat('Superior Temporal Gyrus, Anterior Division');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 10
        ROI = strcat('Superior Temporal Gyrus, Posterior Division');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 11
        ROI = strcat('Middle Temporal Gyrus, Anterior Division');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 12
        ROI = strcat('Middle Temporal Gyrus, Posterior Division');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end         

    elseif z == 13
        ROI = strcat('Middle Temporal Gyrus, Temporooccipital Part');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 14
        ROI = strcat('Inferior Temporal Gyrus, Anterior Division');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 15
        ROI = strcat('Inferior Temporal Gyrus, Posterior Division');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end         

    elseif z == 16
        ROI = strcat('Inferior Temporal Gyrus, Temporooccipital Part');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 17
        ROI = strcat('Postcentral Gyrus');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 18
        ROI = strcat('Superior Parietal Lobule');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 19
        ROI = strcat('Supramarginal Gyrus, Anterior Division');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 20
        ROI = strcat('Supramarginal Gyrus, Posterior Division');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 21
        ROI = strcat('Angular Gyrus');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 22
        ROI = strcat('Lateral Occipital Cortex, Superior Division');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end       

    elseif z == 23
        ROI = strcat('Lateral Occipital Cortex, Inferior Division');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 24
        ROI = strcat('Intracalcarine Cortex');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 25
        ROI = strcat('Frontal Medial Cortex');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 26
        ROI = strcat('Juxtapositional Lobule Cortex (ie. SMA)');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 27
        ROI = strcat('Subcallosal Cortex');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 28
        ROI = strcat('Paracingulate Gyrus');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 29
        ROI = strcat('Cingulate Gyrus, Anterior');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 30
        ROI = strcat('Cingulate Gyrus, Posterior Division');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 31
        ROI = strcat('Precuneous Cortex');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 32
        ROI = strcat('Cuneal Cortex');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 33
        ROI = strcat('Frontal Orbital Cortex');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 34
        ROI = strcat('Parahippocampal Gyrus, Anterior Division');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 35
        ROI = strcat('Parahippocampal Gyrus, Posterior Division');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 36
        ROI = strcat('Lingual Gyrus');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 37
        ROI = strcat('Temporal Fusiform Cortex, Anterior Division');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 38
        ROI = strcat('Temporal Fusiform Cortex');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 39
        ROI = strcat('Temporal Occipital Fusiform Cortex');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 40
        ROI = strcat('Occipital Fusiform Gyrus');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 41
        ROI = strcat('Frontal Operculum Cortex');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 42
        ROI = strcat('Central Opercular Cortex');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 43
        ROI = strcat('Parietal Operculum Cortex');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 44
        ROI = strcat('Planum Polare');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 45
        ROI = strcat('Heschl s gyrus (includes H1 and H2)');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 46
        ROI = strcat('Planum Temporale');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 47
        ROI = strcat('Supracalcarine Cortex');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 48
        ROI = strcat('Occipital Pole');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    end

fprintf(fileID,'%47s %12s %5s %1.4f %5s %1.4f %5s %1.4f %5s %1.4f %4s %1.4f %5s %1.4f\n',ROI,'NDI',' ',Pre_mean_NDI(1,z),' ',Pre_std_NDI(1,z),' ',Post_mean_NDI(1,z),' ',Post_std_NDI(1,z),' ',Delta_NDI,' ',pN1);
%fprintf(fileID,'%16s %12s %5s %1.4f %5s %1.4f %5s %1.4f %5s %1.4f %4s %1.4f %5s %1.4f\n',' ','ODI',' ',PreO_mean1,' ',PreO_std1,' ',PostO_mean1,' ',PostO_std1,' ',Delta_ODI_1,' ',pO1);

end

fclose(fileID);

niftiwrite(HOc_NDI_image,'HOc_NDI_image',image_info)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Now, let's do the Harvard Oxford Subcortical Brain Atlas

    % Vector of values selected from a single HOsc mask
NDI_pre = [];
ODI_pre = [];
FWF_pre = [];
NDI_post = [];
ODI_post = [];
FWF_post = [];

    % Matrices to store the values selected from all HOsc masks
NDI_HOsc_pre = [];
ODI_HOsc_pre = [];
FWF_HOsc_pre = [];
NDI_HOsc_post = [];
ODI_HOsc_post = [];
FWF_HOsc_post = [];

    % Vectors for mean values
Pre_mean_NDI = zeros(1,21);
Pre_mean_ODI = zeros(1,21);
Pre_mean_FWF = zeros(1,21);
Post_mean_NDI = zeros(1,21);
Post_mean_ODI = zeros(1,21);
Post_mean_FWF = zeros(1,21);

    % Vectors for std values
Pre_std_NDI = zeros(1,21);
Pre_std_ODI = zeros(1,21);
Pre_std_FWF = zeros(1,21);
Post_std_NDI = zeros(1,21);
Post_std_ODI = zeros(1,21);
Post_std_FWF = zeros(1,21);

% Counts the number of voxels in each ROI
HOsc_NDI_numbers = [];

sizer = size(sPre_NDI);

for z = 1:21
    for i = 1:sizer(1,1)
        for j = 1:sizer(1,2)
            for k = 1:sizer(1,3)
                if new_HOsc(i,j,k) == z
                    if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                    elseif GM(i,j,k) == 0

                    else
                        NDI_pre = [NDI_pre sPre_NDI(i,j,k)];
                        NDI_post = [NDI_post sPost_NDI(i,j,k)];
                    end
                end
            end
        end
    end
    
    a = length(NDI_HOsc_pre);
    b = length(NDI_pre);

    HOsc_NDI_numbers = [HOsc_NDI_numbers; b];
    
    Pre_mean_NDI(1,z) = mean(NDI_pre);
    Post_mean_NDI(1,z) = mean(NDI_post);

    Pre_std_NDI(1,z) = std(NDI_pre);
    Post_std_NDI(1,z) = std(NDI_post);

    if a > b

        NDI_HOsc_pre = [NDI_HOsc_pre; zeros(1,a)];
        NDI_HOsc_post = [NDI_HOsc_post; zeros(1,a)];

        for i = 1:b

            NDI_HOsc_pre(z,i) = NDI_pre(1,i);
            NDI_HOsc_post(z,i) = NDI_post(1,i);

        end

    else

        if z == 1

            NDI_HOsc_pre(1,:) = NDI_pre(1,:);
            NDI_HOsc_post(1,:) = NDI_post(1,:);

        else

        NDI_HOsc_ph1 = zeros(z,b);
        NDI_HOsc_ph2 = zeros(z,b);

            for i = 1:(z-1)
                for j = 1:a
                    NDI_HOsc_ph1(i,j) = NDI_HOsc_pre(i,j);
                    NDI_HOsc_ph2(i,j) = NDI_HOsc_post(i,j);
                end
            end

            for i = 1:a

                NDI_HOsc_ph1(z,i) = NDI_pre(1,i);
                NDI_HOsc_ph2(z,i) = NDI_post(1,i);

            end

            NDI_HOsc_pre = NDI_HOsc_ph1;
            NDI_HOsc_post = NDI_HOsc_ph2;

        end
    end

    NDI_pre = [];
    NDI_post = [];

end

% Calculate the statistics

fileID = fopen(strcat('ROI_HOsc_results_',num2str(control_num),'.txt'),'w');
fprintf(fileID,'%47s %12s %12s %12s %12s %12s %12s %12s\n','ROI','Index','Pre_Mean','Pre_Std','Post_Mean','Post_Std','Difference','p-value');

HOsc_NDI_image = zeros(sizer(1,1),sizer(1,2),sizer(1,3));

for z = 1:21

    Delta_NDI = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));

    values_pre = [];
    values_post = [];

    for i = 1:HOsc_NDI_numbers(z,1)
        values_pre = [values_pre NDI_HOsc_pre(z,i)];
        values_post = [values_post NDI_HOsc_post(z,i)];
    end
    
    [hN1,pN1] = ttest2(values_pre,values_post);

    if z == 1
        ROI = strcat('Left Cerebral White Matter');
        
        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOsc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOsc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end        

    elseif z == 2
        ROI = strcat('Left Cerebral Cortex');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOsc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOsc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end  

    elseif z == 3
        ROI = strcat('Left Lateral Ventrilce');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOsc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOsc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end  

    elseif z == 4
        ROI = strcat('Left Thalamus');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOsc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOsc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end  

    elseif z == 5
        ROI = strcat('Left Caudate');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOsc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOsc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 6
        ROI = strcat('Left Putamen');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOsc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOsc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 7 
        ROI = strcat('Left Pallidum');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOsc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOsc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 8
        ROI = strcat('Brainstem');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOsc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOsc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 9
        ROI = strcat('Left Hippocampus');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOsc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOsc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 10
        ROI = strcat('Left Amygdala');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOsc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOsc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 11
        ROI = strcat('Left Accumbens');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOsc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOsc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 12
        ROI = strcat('Right Cerebral White Matter');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOsc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOsc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end         

    elseif z == 13
        ROI = strcat('Right Cerebral Cortex');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOsc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOsc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 14
        ROI = strcat('Right Lateral Ventricle');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOsc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOsc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 15
        ROI = strcat('Right Thalamus');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOsc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOsc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end         

    elseif z == 16
        ROI = strcat('Right Caudate');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOsc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOsc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 17
        ROI = strcat('Right Putamen');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOsc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOsc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 18
        ROI = strcat('Right Pallidum');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOsc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOsc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 19
        ROI = strcat('Right Hippocampus');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOsc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOsc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 20
        ROI = strcat('Right Amygdala');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOsc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOsc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 21
        ROI = strcat('Right Accumbens');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_HOsc(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif GM(i,j,k) == 0

                        else
                            HOsc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 
    end

fprintf(fileID,'%47s %12s %5s %1.4f %5s %1.4f %5s %1.4f %5s %1.4f %4s %1.4f %5s %1.4f\n',ROI,'NDI',' ',Pre_mean_NDI(1,z),' ',Pre_std_NDI(1,z),' ',Post_mean_NDI(1,z),' ',Post_std_NDI(1,z),' ',Delta_NDI,' ',pN1);
%fprintf(fileID,'%16s %12s %5s %1.4f %5s %1.4f %5s %1.4f %5s %1.4f %4s %1.4f %5s %1.4f\n',' ','ODI',' ',PreO_mean1,' ',PreO_std1,' ',PostO_mean1,' ',PostO_std1,' ',Delta_ODI_1,' ',pO1);

end

fclose(fileID);

niftiwrite(HOsc_NDI_image,'HOsc_NDI_image',image_info)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


