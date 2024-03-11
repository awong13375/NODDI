%-----------------------------------------------------------------------
% Job saved on 29-Jan-2024 22:56:04 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.spatial.coreg.estimate.ref = {'/usr/local/fsl/data/standard/MNI152_T1_1mm.nii,1'};
matlabbatch{1}.spm.spatial.coreg.estimate.source = {'/mnt/c/WSL2_dir/NODDISAH_11/NODDI_processing/anat_seq_brain_mask.nii,1'};
%%
matlabbatch{1}.spm.spatial.coreg.estimate.other = {
                                                   '/mnt/c/WSL2_dir/NODDISAH_11/NODDI_processing/Case1_error_code.nii,1'
                                                   '/mnt/c/WSL2_dir/NODDISAH_11/NODDI_processing/Case1_fibredirs_xvec.nii,1'
                                                   '/mnt/c/WSL2_dir/NODDISAH_11/NODDI_processing/Case1_fibredirs_yvec.nii,1'
                                                   '/mnt/c/WSL2_dir/NODDISAH_11/NODDI_processing/Case1_fibredirs_zvec.nii,1'
                                                   '/mnt/c/WSL2_dir/NODDISAH_11/NODDI_processing/Case1_ficvf.nii,1'
                                                   '/mnt/c/WSL2_dir/NODDISAH_11/NODDI_processing/Case1_fiso.nii,1'
                                                   '/mnt/c/WSL2_dir/NODDISAH_11/NODDI_processing/Case1_fmin.nii,1'
                                                   '/mnt/c/WSL2_dir/NODDISAH_11/NODDI_processing/Case1_kappa.nii,1'
                                                   '/mnt/c/WSL2_dir/NODDISAH_11/NODDI_processing/Case1_odi.nii,1'
                                                   };
%%
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
