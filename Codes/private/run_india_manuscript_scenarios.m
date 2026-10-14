function result = run_india_manuscript_scenarios(runName)
%RUN_INDIA_MANUSCRIPT_SCENARIOS Six-parameter projections after diagnostic checks.
assert(ischar(runName) && ~isempty(regexp(runName,'^[A-Za-z0-9][A-Za-z0-9_-]*$','once')));
helperFolder=fileparts(mfilename('fullpath'));folder=fileparts(helperFolder);
input=fullfile(folder,'outputs',runName);loaded=load(fullfile(input,'calibration.mat'),'result');
fits=loaded.result;
assert(strcmp(fits.data.modelVersion,'india10-manuscript-lsq-2026-10-13-v2'), ...
    'india10:ModelVersionChanged','This release requires a new calibration run; previous outputs are preserved.');
% Do not apply old acceptance of a different, two-parameter model to this fit.
assert(fits.checks.scenarioGenerationAllowed && strcmp(fits.options.mode,'full'), ...
    'india10:ScenariosBlocked','Calibration checks have not passed. Read the calibration REPORT.txt.');
assert(all(fits.completed) && all(fits.converged) && size(fits.parameters,2)==6);
locations={folder,helperFolder};labels={'entrypoints','private'};
for k=1:numel(locations)
    files=[dir(fullfile(input,'source',labels{k},'*.m'));dir(fullfile(input,'source',labels{k},'*.cs'))];
    assert(~isempty(files),'Missing source snapshot.');
    for j=1:numel(files)
        assert(strcmp(fileread(fullfile(locations{k},files(j).name)), ...
            fileread(fullfile(input,'source',labels{k},files(j).name))), ...
            'Source changed since calibration: %s. Refit with consistent code.',files(j).name);
    end
end
output=fullfile(input,'scenarios');assert(exist(output,'dir')~=7,'Scenario output already exists.');mkdir(output);
settings=india10_scenario_settings();settings.substeps=fits.options.substeps;
settings.backend='matlab';
n=fits.nBootstrap;ny=numel(settings.years);annual=zeros(ny,12,n);totals=zeros(12,n);
annualTrue=zeros(ny,12,n);unassignedLoss=zeros(12,n);
equationBalanceError=zeros(12,n);
result.settings=settings;result.parameters=fits.parameters;result.data=fits.data;
[~,history]=india10_simulate(fits.base.parameters,fits.data,settings.substeps);
result.stepRefinementError=zeros(1,12);
for j=1:12
    s=settings.scenarios(j);initial=history.finalStatePeople/1e6;
    coarse=india10_scenario_simulate(fits.base.parameters,initial,s,settings,fits.data);
    fine=india10_scenario_simulate(fits.base.parameters,initial,s,settings,fits.data,2*settings.substeps);
    result.stepRefinementError(j)=max(abs(coarse.trueCases(:)-fine.trueCases(:))./max(1,fine.trueCases(:)));
    assert(result.stepRefinementError(j)<.01,'Forecast timestep refinement failed.');
end
for b=1:n
    p=fits.parameters(b,:);[mu,history]=india10_simulate(p,fits.data,settings.substeps);
    assert(max(abs(mu-fits.fits{b}.predicted)./max(1,mu))<1e-7,'Fitted history changed.');
    for j=1:12
        projection=india10_scenario_simulate(p,history.finalStatePeople/1e6,settings.scenarios(j),settings,fits.data);
        assert(projection.maxEquationBalanceErrorPeople<1e-3 && projection.minimumCompartmentPeople>=-1e-4);
        annual(:,j,b)=sum(projection.reportedCases,2);totals(j,b)=sum(annual(:,j,b));
        annualTrue(:,j,b)=sum(projection.trueCases,2);
        unassignedLoss(j,b)=projection.cumulativeUnassignedV1LossPeople;
        equationBalanceError(j,b)=projection.maxEquationBalanceErrorPeople;
    end
    if mod(b,25)==0 || b==n,fprintf('Six-parameter scenarios: %d/%d parameter vectors.\n',b,n);end
end
result.annualReported=annual;result.cumulativeReported=totals;
result.annualTrue=annualTrue;
result.cumulativeUnassignedV1LossPeople=unassignedLoss;
result.equationBalanceErrorPeople=equationBalanceError;
result.percentAverted=100*bsxfun(@rdivide,bsxfun(@minus,totals(1,:),totals),totals(1,:));
result.intervalMeaning='Bootstrap percentiles of fitted model means; parameter uncertainty only.';
result=india10_scenario_report(result,output);
fprintf('Six-parameter scenarios saved in %s\n',output);
end
