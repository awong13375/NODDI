% Making maps of ethanol NODDI differences with one control scan at once
% Jesse M Klostranec MD PhD FRCPC
% February 24, 2024

clc
clear

disp(' Ethanol-NODDI DTI Program for White Matter ')
disp(' ')
disp(' February 24, 2024 ')
disp(' By Jesse M Klostranec MD PhD FRCPC ')
disp(' ')
disp(' Press any key to continue ...')
pause

disp('---------')
disp(' Ethanol NODDI Analysis Using the JHU White Matter Atlas ')
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
place_holder = niftiread(strcat('rdti_tensor_',num2str(control_num),'.nii'));
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
place_holder = niftiread(strcat('rdti_tensor.nii'));
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

% Load the JHU white matter atlases

WMlabels = double(niftiread('JHU-ICBM-labels-1mm.nii'));
sizeWMl = size(WMlabels);
WMtracts = double(niftiread('JHU-ICBM-tracts-maxprob-thr0-1mm.nii'));

% Load the MNI atlas

MNIt = double(niftiread('MNI-maxprob-thr0-1mm.nii'));

% Load the segmented white matter image

WM = double(niftiread('T1_masked_wm.nii'));

% Match the orientation and sizes of the two images:

% These are the following coordinates for the different scans:
%   Max: voxel [0 0 0] = (-98 -134 -72)
%              [196 232 188] = (98.5 98.5 116.5)
%
%   MNI: voxel [0 0 0] = (90 -126 -72)
%              [181 217 181] = (-91.5 91.5 109.5)

new_WMl = zeros(size(sPost_NDI));
sizer = size(new_WMl);
new_WMt = zeros(size(sPost_NDI));
new_MNI = zeros(size(sPost_NDI));

for i = 9:(9 - 1 + sizeWMl(1,1))
    for j = 9:(9 - 1 + sizeWMl(1,2))
        for k = 1:sizeWMl(1,3)
            new_WMl(i,j,k) = WMlabels((i-8),(j-8),k);
            new_WMt(i,j,k) = WMtracts((i-8),(j-8),k);
            new_MNI(i,j,k) = MNIt((i-8),(j-8),k);
        end
    end
end
sizer = [];

% For the MNI and JHU white matter atlases, the following values correspond to the following
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
% of collected values will be different, ie. for MNI - 9 ROIs, for JHU
% white matter labels 50 ROIs, for JHU white matter tracts - 20 ROIs

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

                    elseif WM(i,j,k) == 0

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

                    elseif WM(i,j,k) == 0

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

                    elseif WM(i,j,k) == 0

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

                        elseif WM(i,j,k) == 0

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

                        elseif WM(i,j,k) == 0

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

                        elseif WM(i,j,k) == 0

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

                        elseif WM(i,j,k) == 0

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

                        elseif WM(i,j,k) == 0

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

                        elseif WM(i,j,k) == 0

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

                        elseif WM(i,j,k) == 0

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

                        elseif WM(i,j,k) == 0

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

                        elseif WM(i,j,k) == 0

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

                        elseif WM(i,j,k) == 0

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

                        elseif WM(i,j,k) == 0

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

                        elseif WM(i,j,k) == 0

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

                        elseif WM(i,j,k) == 0

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

                        elseif WM(i,j,k) == 0

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

                        elseif WM(i,j,k) == 0

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

                        elseif WM(i,j,k) == 0

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

                        elseif WM(i,j,k) == 0

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

                        elseif WM(i,j,k) == 0

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


% Now, let's do the JHU white matter lables Brain Atlas

    % Vector of values selected from a single WMlabels mask
NDI_pre = [];
ODI_pre = [];
FWF_pre = [];
NDI_post = [];
ODI_post = [];
FWF_post = [];

    % Matrices to store the values selected from all WMlabels masks
NDI_WMl_pre = [];
ODI_WMl_pre = [];
FWF_WMl_pre = [];
NDI_WMl_post = [];
ODI_WMl_post = [];
FWF_WMl_post = [];

    % Vectors for mean values
Pre_mean_NDI = zeros(1,50);
Pre_mean_ODI = zeros(1,50);
Pre_mean_FWF = zeros(1,50);
Post_mean_NDI = zeros(1,50);
Post_mean_ODI = zeros(1,50);
Post_mean_FWF = zeros(1,50);

    % Vectors for std values
Pre_std_NDI = zeros(1,50);
Pre_std_ODI = zeros(1,50);
Pre_std_FWF = zeros(1,50);
Post_std_NDI = zeros(1,50);
Post_std_ODI = zeros(1,50);
Post_std_FWF = zeros(1,50);

% Counts the number of voxels in each ROI
WMl_NDI_numbers = [];

sizer = size(sPre_NDI);

for z = 1:50
    for i = 1:sizer(1,1)
        for j = 1:sizer(1,2)
            for k = 1:sizer(1,3)
                if new_WMl(i,j,k) == z
                    if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                    elseif WM(i,j,k) == 0

                    else
                        NDI_pre = [NDI_pre sPre_NDI(i,j,k)];
                        NDI_post = [NDI_post sPost_NDI(i,j,k)];
                    end
                end
            end
        end
    end
    
    a = length(NDI_WMl_pre);
    b = length(NDI_pre);

    WMl_NDI_numbers = [WMl_NDI_numbers; b];
    
    Pre_mean_NDI(1,z) = mean(NDI_pre);
    Post_mean_NDI(1,z) = mean(NDI_post);

    Pre_std_NDI(1,z) = std(NDI_pre);
    Post_std_NDI(1,z) = std(NDI_post);

    if a > b

        NDI_WMl_pre = [NDI_WMl_pre; zeros(1,a)];
        NDI_WMl_post = [NDI_WMl_post; zeros(1,a)];

        for i = 1:b

            NDI_WMl_pre(z,i) = NDI_pre(1,i);
            NDI_WMl_post(z,i) = NDI_post(1,i);

        end

    else

        if z == 1

            NDI_WMl_pre(1,:) = NDI_pre(1,:);
            NDI_WMl_post(1,:) = NDI_post(1,:);

        else

        NDI_WMl_ph1 = zeros(z,b);
        NDI_WMl_ph2 = zeros(z,b);

            for i = 1:(z-1)
                for j = 1:a
                    NDI_WMl_ph1(i,j) = NDI_WMl_pre(i,j);
                    NDI_WMl_ph2(i,j) = NDI_WMl_post(i,j);
                end
            end

            for i = 1:a

                NDI_WMl_ph1(z,i) = NDI_pre(1,i);
                NDI_WMl_ph2(z,i) = NDI_post(1,i);

            end

            NDI_WMl_pre = NDI_WMl_ph1;
            NDI_WMl_post = NDI_WMl_ph2;

        end
    end

    NDI_pre = [];
    NDI_post = [];

end

% Calculate the statistics

fileID = fopen(strcat('ROI_WMl_results_',num2str(control_num),'.txt'),'w');
fprintf(fileID,'%47s %12s %12s %12s %12s %12s %12s %12s\n','ROI','Index','Pre_Mean','Pre_Std','Post_Mean','Post_Std','Difference','p-value');

WMl_NDI_image = zeros(sizer(1,1),sizer(1,2),sizer(1,3));

for z = 1:50

    Delta_NDI = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));

    values_pre = [];
    values_post = [];

    for i = 1:WMl_NDI_numbers(z,1)
        values_pre = [values_pre NDI_WMl_pre(z,i)];
        values_post = [values_post NDI_WMl_post(z,i)];
    end
    
    [hN1,pN1] = ttest2(values_pre,values_post);
    [hO1,pO1] = ttest2(Ovalues_pre,Ovalues_post);

    if z == 1
        ROI = strcat('Middle Cerebellar Peduncle');
        
        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end        

    elseif z == 2
        ROI = strcat('Pontine Crossing Tract');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end  

    elseif z == 3
        ROI = strcat('Genu of Corpus Callosum');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end  

    elseif z == 4
        ROI = strcat('Body of Corpus Callosum');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end  

    elseif z == 5
        ROI = strcat('Splenium');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 6
        ROI = strcat('Fornix (column and body of fornix)');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 7 
        ROI = strcat('Right Corticospinal Tract');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 8
        ROI = strcat('Left Corticospinal Tract');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 9
        ROI = strcat('Right Medial Lemniscus');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 10
        ROI = strcat('Left Medial Lemniscus');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 11
        ROI = strcat('Right Inferior Cerebellar Peduncle');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 12
        ROI = strcat('Left Inferior Cerebellar Peduncle');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end         

    elseif z == 13
        ROI = strcat('Right Superior Cerebellar Peduncle');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 14
        ROI = strcat('Left Superior Cerebellar Peduncle');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 15
        ROI = strcat('Right Cerebral Peduncle');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end         

    elseif z == 16
        ROI = strcat('Left Cerebral Peduncle');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 17
        ROI = strcat('Right Anterior Limb Internal Capsule');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 18
        ROI = strcat('Left Anterior Limb Internal Capsule');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 19
        ROI = strcat('Right Posterior Limb Internal Capsule');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 20
        ROI = strcat('Left Posterior Limb Internal Capsule');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 21
        ROI = strcat('Right Retrolenticular Internal Capsule');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 22
        ROI = strcat('Left Retrolenticular Internal Capsule');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end       

    elseif z == 23
        ROI = strcat('Right Anterior Corona Radiata');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 24
        ROI = strcat('Left Anterior Corona Radiata');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 25
        ROI = strcat('Right Superior Corona Radiata');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 26
        ROI = strcat('Left Superior Corona Radiata');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 27
        ROI = strcat('Right Posterior Corona Radiata');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 28
        ROI = strcat('Left Posterior Corona Radiata');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 29
        ROI = strcat('Right Posterior Thalamic Radiation');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 30
        ROI = strcat('Left Posterior Thalamic Radiation');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 31
        ROI = strcat('Right Sagittal Stratum');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 32
        ROI = strcat('Left Sagittal Stratum');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 33
        ROI = strcat('Right External Capsule');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 34
        ROI = strcat('Left External Capsule');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 35
        ROI = strcat('Right Cingulum');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 36
        ROI = strcat('Left Cingulum');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 37
        ROI = strcat('Right Cingulum (hippocampus)');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 38
        ROI = strcat('Left Cingulum (hippocampus)');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 39
        ROI = strcat('Right Fornix (cres)/Stria Terminalis');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 40
        ROI = strcat('Left Fornix (cres)/Stria Terminalis');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 41
        ROI = strcat('Right Superior Longitudinal Fasciculus');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 42
        ROI = strcat('Left Superior Longitudinal Fasciculus');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 43
        ROI = strcat('Right Superior Fronto-occipital Fasciculus');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 44
        ROI = strcat('Left Superior Fronto-occipital Fasciculus');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 45
        ROI = strcat('Right Inferior Fronto-occipital Fasciculus');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 46
        ROI = strcat('Left Inferior Fronto-occipital Fasciculus');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 47
        ROI = strcat('Right Uncinate Fasciculus');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 48
        ROI = strcat('Left Uncinate Fasciculus');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 49
        ROI = strcat('Right Tapetum');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 50
        ROI = strcat('Left Tapetum');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMl(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMl_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
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

niftiwrite(WMl_NDI_image,'WMl_NDI_image',image_info)


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
NDI_WMt_pre = [];
ODI_WMt_pre = [];
FWF_WMt_pre = [];
NDI_WMt_post = [];
ODI_WMt_post = [];
FWF_WMt_post = [];

    % Vectors for mean values
Pre_mean_NDI = zeros(1,20);
Pre_mean_ODI = zeros(1,20);
Pre_mean_FWF = zeros(1,20);
Post_mean_NDI = zeros(1,20);
Post_mean_ODI = zeros(1,20);
Post_mean_FWF = zeros(1,20);

    % Vectors for std values
Pre_std_NDI = zeros(1,20);
Pre_std_ODI = zeros(1,20);
Pre_std_FWF = zeros(1,20);
Post_std_NDI = zeros(1,20);
Post_std_ODI = zeros(1,20);
Post_std_FWF = zeros(1,20);

% Counts the number of voxels in each ROI
WMt_NDI_numbers = [];

sizer = size(sPre_NDI);

for z = 1:21
    for i = 1:sizer(1,1)
        for j = 1:sizer(1,2)
            for k = 1:sizer(1,3)
                if new_WMt(i,j,k) == z
                    if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                    elseif WM(i,j,k) == 0

                    else
                        NDI_pre = [NDI_pre sPre_NDI(i,j,k)];
                        NDI_post = [NDI_post sPost_NDI(i,j,k)];
                    end
                end
            end
        end
    end
    
    a = length(NDI_WMt_pre);
    b = length(NDI_pre);

    WMt_NDI_numbers = [WMt_NDI_numbers; b];
    
    Pre_mean_NDI(1,z) = mean(NDI_pre);
    Post_mean_NDI(1,z) = mean(NDI_post);

    Pre_std_NDI(1,z) = std(NDI_pre);
    Post_std_NDI(1,z) = std(NDI_post);

    if a > b

        NDI_WMt_pre = [NDI_WMt_pre; zeros(1,a)];
        NDI_WMt_post = [NDI_WMt_post; zeros(1,a)];

        for i = 1:b

            NDI_WMt_pre(z,i) = NDI_pre(1,i);
            NDI_WMt_post(z,i) = NDI_post(1,i);

        end

    else

        if z == 1

            NDI_WMt_pre(1,:) = NDI_pre(1,:);
            NDI_WMt_post(1,:) = NDI_post(1,:);

        else

        NDI_WMt_ph1 = zeros(z,b);
        NDI_WMt_ph2 = zeros(z,b);

            for i = 1:(z-1)
                for j = 1:a
                    NDI_WMt_ph1(i,j) = NDI_WMt_pre(i,j);
                    NDI_WMt_ph2(i,j) = NDI_WMt_post(i,j);
                end
            end

            for i = 1:a

                NDI_WMt_ph1(z,i) = NDI_pre(1,i);
                NDI_WMt_ph2(z,i) = NDI_post(1,i);

            end

            NDI_WMt_pre = NDI_WMt_ph1;
            NDI_WMt_post = NDI_WMt_ph2;

        end
    end

    NDI_pre = [];
    NDI_post = [];

end

% Calculate the statistics

fileID = fopen(strcat('ROI_WMt_results_',num2str(control_num),'.txt'),'w');
fprintf(fileID,'%47s %12s %12s %12s %12s %12s %12s %12s\n','ROI','Index','Pre_Mean','Pre_Std','Post_Mean','Post_Std','Difference','p-value');

WMt_NDI_image = zeros(sizer(1,1),sizer(1,2),sizer(1,3));

for z = 1:20

    Delta_NDI = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));

    values_pre = [];
    values_post = [];

    for i = 1:WMt_NDI_numbers(z,1)
        values_pre = [values_pre NDI_WMt_pre(z,i)];
        values_post = [values_post NDI_WMt_post(z,i)];
    end
    
    [hN1,pN1] = ttest2(values_pre,values_post);

    if z == 1
        ROI = strcat('Left Anterior Thalamic Radiation');
        
        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMt(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMt_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end        

    elseif z == 2
        ROI = strcat('Right Anterior Thalamic Radiation');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMt(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMt_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end  

    elseif z == 3
        ROI = strcat('Left Corticospinal Tract');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMt(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMt_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end  

    elseif z == 4
        ROI = strcat('Right Corticospinal Tract');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMt(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            HOsc_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end  

    elseif z == 5
        ROI = strcat('Left Cingulum (Cingulate Gyrus)');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMt(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMt_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 6
        ROI = strcat('Right Cingulum (Cingulate Gyrus)');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMt(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMt_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 7 
        ROI = strcat('Left Cingulum (Hippocampus)');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMt(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMt_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 8
        ROI = strcat('Right Cingulum (Hippocampus)');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMt(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMt_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 9
        ROI = strcat('Forceps Major');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMt(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMt_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 10
        ROI = strcat('Forceps Minor');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMt(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMt_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 11
        ROI = strcat('Left Inferior Fronto-occipital Fasciculus');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMt(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMt_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 12
        ROI = strcat('Right Inferior Fronto-occipital Fasciculus');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMt(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMt_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end         

    elseif z == 13
        ROI = strcat('Left Inferior Longitudinal Fasciculus');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMt(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMt_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 14
        ROI = strcat('Right Inferior Longitudinal Fasciculus');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMt(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMt_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 15
        ROI = strcat('Left Superior Longitudinal Fasciculus');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMt(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMt_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end         

    elseif z == 16
        ROI = strcat('Right Superior Longitudinal Fasciculus');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMt(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMt_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 17
        ROI = strcat('Left Uncinate Fasciculus');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMt(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMt_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 18
        ROI = strcat('Right Uncinate Fasciculus');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMt(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMt_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 19
        ROI = strcat('Left SLF (Temporal Part)');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMt(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMt_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
                        end
                    end
                end
            end
        end 

    elseif z == 20
        ROI = strcat('Right SLF (Temporal Part)');

        for i = 1:sizer(1,1)
            for j = 1:sizer(1,2)
                for k = 1:sizer(1,3)
                    if new_WMt(i,j,k) == z
                        if sPre_NDI(i,j,k) == 0 || sPost_NDI(i,j,k) == 0 || isnan(sPre_NDI(i,j,k)) == 1 || isnan(sPost_NDI(i,j,k)) == 1

                        elseif WM(i,j,k) == 0

                        else
                            WMt_NDI_image(i,j,k) = (Pre_mean_NDI(1,z) - Post_mean_NDI(1,z));
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

niftiwrite(WMt_NDI_image,'WMt_NDI_image',image_info)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

