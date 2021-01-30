function [tracts_unaltered] = loadTractNemo_rightonly_noalteration()
% Load voxelwise ChaCo scores for each tract & edit with appropriate mask.
    nemodir ='/home/emo4002/colossus_shared3/pons_sfmodelling/processing/BrainStemAtlas/NeMo/*'
    a=dir(nemodir);
    tractnames=a.name;
    for i=3:25
        tracts_unaltered{i-2}=read_avw(strcat(nemodir, a(i).name))
    end
    
end

