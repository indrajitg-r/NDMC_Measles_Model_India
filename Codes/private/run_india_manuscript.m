function result = run_india_manuscript(nBootstrap,runName,mode,backend)
%RUN_INDIA_MANUSCRIPT Six-parameter OLS calibration and 20%-SD Normal bootstrap.
% Example: run_india_manuscript(1000,'india_lsq_1000');
% 'check' reduces the fit budget, permits <=3 replicates and never unlocks scenarios.
if nargin<1,nBootstrap=1000;end
if nargin<2,runName='india_lsq_1000';end
if nargin<3,mode='full';end
if nargin<4,backend='auto';end
assert(isscalar(nBootstrap) && isfinite(nBootstrap) && nBootstrap>=0 && nBootstrap==round(nBootstrap));
assert(ischar(runName) && ~isempty(regexp(runName,'^[A-Za-z0-9][A-Za-z0-9_-]*$','once')));
assert(any(strcmp(mode,{'full','check'})));
helperFolder=fileparts(mfilename('fullpath'));
folder=fileparts(helperFolder);
output=fullfile(folder,'outputs',runName);
assert(exist(output,'dir')~=7,'Run already exists. Choose a new run name to preserve results.');
mkdir(output);
options.substeps=4;options.maxIterations=300;options.maxEvaluations=4000;
options.startCount=8;options.bootstrapStarts=3;options.seed=13102026;options.mode=mode;
if strcmp(mode,'check')
    assert(nBootstrap<=3,'Check mode allows at most three bootstrap replicates.');
    options.maxIterations=60;options.maxEvaluations=900;options.startCount=2;options.bootstrapStarts=2;
end
data=india10_manuscript_data();
[data,validation]=india10_backend(data,backend,fullfile(output,'build'));
% Keep the exact implementation used for a reproducible source audit.
locations={folder,helperFolder};labels={'entrypoints','private'};
for k=1:numel(locations)
    dest=fullfile(output,'source',labels{k});mkdir(dest);
    files=[dir(fullfile(locations{k},'*.m'));dir(fullfile(locations{k},'*.cs'))];
    for j=1:numel(files),copyfile(fullfile(locations{k},files(j).name),dest);end
end
result.data=data;result.options=options;result.backendValidation=validation;
result.nBootstrap=nBootstrap;result.started=datestr(now,30);result.matlabVersion=version;
stream=RandStream('mt19937ar','Seed',options.seed);
starts=repmat(data.start,options.startCount,1).*exp(.65*randn(stream,options.startCount,6));
starts(1,:)=data.start;
result.base=india10_lsq_fit(data,options,starts);
fprintf('Base OLS fit: SSE %.10g, exitflag %d, agreeing starts %d/%d.\n', ...
    result.base.sse,result.base.exitflag,result.base.agreeingStarts,options.startCount);
result.parameters=NaN(nBootstrap,6);result.fits=cell(nBootstrap,1);
result.observations=NaN(24,nBootstrap);result.completed=false(nBootstrap,1);
result.converged=false(nBootstrap,1);
result.bootstrapStepError=NaN(nBootstrap,1);
save(fullfile(output,'checkpoint.mat'),'result','-v7');
if ~result.base.converged && nBootstrap>0
    error('india10:CalibrationFailed','Generating fit did not converge. Check checkpoint; no bootstrap or scenarios run.');
end
for b=1:nBootstrap
    stream=RandStream('mt19937ar','Seed',options.seed+b);
    synthetic=data;
    synthetic.cases=result.base.predicted.*(1+.20*randn(stream,24,1));
    result.observations(:,b)=synthetic.cases;
    starts=repmat(result.base.parameters,options.bootstrapStarts,1).* ...
        exp(.3*randn(stream,options.bootstrapStarts,6));starts(1,:)=result.base.parameters;
    fit=india10_lsq_fit(synthetic,options,starts);
    if ~fit.converged || fit.agreeingStarts<2
        retry=options;retry.maxIterations=2*options.maxIterations;retry.maxEvaluations=2*options.maxEvaluations;
        extra=repmat(result.base.parameters,3,1).*exp(.6*randn(stream,3,6));
        fit=india10_lsq_fit(synthetic,retry,[fit.parameters;starts;extra]);
    end
    refined=india10_predict(fit.parameters,data,2*options.substeps);
    result.bootstrapStepError(b)=max(abs(fit.predicted-refined)./max(1,refined));
    result.fits{b}=fit;result.parameters(b,:)=fit.parameters;
    result.converged(b)=fit.converged && fit.agreeingStarts>=2 && result.bootstrapStepError(b)<.01;
    result.completed(b)=true;
    save(fullfile(output,'checkpoint.mat'),'result','-v7');
    fprintf('OLS bootstrap %d/%d: SSE %.8g, converged/agreed %d.\n',b,nBootstrap,fit.sse,result.converged(b));
end
result.checks=india10_lsq_checks(result);
if strcmp(mode,'check'),result.checks.scenarioGenerationAllowed=false;end
result.finished=datestr(now,30);
save(fullfile(output,'calibration.mat'),'result','-v7');
india10_lsq_report(result,output);
fprintf('Read %s\n',fullfile(output,'REPORT.txt'));
end
