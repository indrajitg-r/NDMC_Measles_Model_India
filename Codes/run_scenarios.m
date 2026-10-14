function scenarios = run_scenarios(runName)
%RUN_SCENARIOS Manuscript scenarios 0--11, using a completed bootstrap run.
% In this folder, run: scenarios = run_scenarios;
% For a named run:     scenarios = run_scenarios('india_1000_repeat');
% Requires successful calibration/identifiability and fit-adequacy checks.
% Results: outputs/<runName>/scenarios/.
if nargin<1,runName='india_1000';end
assert(ischar(runName) && ~isempty(regexp(runName,'^[A-Za-z0-9][A-Za-z0-9_-]*$','once')));
folder=fileparts(mfilename('fullpath'));
assert(exist(fullfile(folder,'outputs',runName,'calibration.mat'),'file')==2, ...
    'india10:MissingCalibration','Run run_bootstrap_1000 with this name first.');
scenarios=run_india_manuscript_scenarios(runName);
end
