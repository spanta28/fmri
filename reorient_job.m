% List of open inputs
% Reorient Images: Reorientation Matrix - cfg_entry
nrun = 1; % enter the number of runs here
jobfile = {'/computation/reorient.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(1, nrun);
spm('defaults', 'PET');
spm_jobman('run', jobs, inputs{:});
