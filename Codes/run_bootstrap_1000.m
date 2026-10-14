function fit = run_bootstrap_1000(runName,backend)
%RUN_BOOTSTRAP_1000 Calibrate six parameters, then generate and fit 1000 samples.
% In this folder, run: fit = run_bootstrap_1000;
% Optional new name:   fit = run_bootstrap_1000('india_1000_repeat');
% Reference MATLAB:   fit = run_bootstrap_1000('india_1000_matlab','matlab');
% Results: outputs/<runName>/REPORT.txt and calibration.mat.
% Scenarios are a separate step: run_scenarios(runName).
if nargin<1,runName='india_1000';end
if nargin<2,backend='auto';end
fit=run_india_manuscript(1000,runName,'full',backend);
end
