matlabbatch{1}.spm.util.reorient.srcfiles = {'input_file,1'};
a=load('/computation/transform.mat');
matlabbatch{1}.spm.util.reorient.transform.transM = a.M;
matlabbatch{1}.spm.util.reorient.prefix = '';