disconnectivitydir='/home/emo4002/colossus_shared3/pons_sfmodelling/processing/disconnectivity/NeMo2_outputs/july20_voxelwise_flippedAfterNeMo/'

clear subject
for i=1:23
    subject_nemo=read_avw(strcat(disconnectivitydir, 'SUB',num2str(i), '_lesion1mm_nemo_output_chacovol_res2mm_mean_right.nii.gz'));
    subject_nemo=reshape(subject_nemo, [1 902629]);
    subject{i}=subject_nemo;
end


subject=cell2mat(subject);

[coeff,score,latent,tsquared,explaind]=pca(subject);