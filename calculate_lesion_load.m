function [lesionload]=calculate_lesion_load(A,B)
% From Zhu et al., 2010: lesion load of the corticospinal tract predicts
% motor impairment in chronic stroke.
% Raw lesion load is calculated by overlaying the lesion map (binarized)
% on top of the probabilistic tract and summing the probabilities (intensities) of all intersecting voxels.

% this is preferable to using the dice coefficient since it takes into
% account probability of being part of the CST.
% inputs:
%     A = lesion mask (binary)
%     B = tract probabilistic ROI 
%
% output:
%     lesionload = probability 
%A=lesion;
%B=atlas;
    nmax= nnz(A & B);%nmax = total number of intersecting voxels between lesion map and tract
    overlap=A&B;
    tract_overlap=B(overlap);
    lesionload=0;
    for i=1:nmax
        lesionload=lesionload+tract_overlap(i);
        value(i)=lesionload;
    end
    if nmax==0
        lesionload=0;
    else
        lesionload=lesionload/nnz(B);
    end
end