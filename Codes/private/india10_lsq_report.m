function india10_lsq_report(r,output)
%INDIA10_LSQ_REPORT Export labelled OLS estimates and identifiability warnings.
fid=fopen(fullfile(output,'parameters.csv'),'w');assert(fid>=0);
fprintf(fid,'parameter,point,bootstrap_p025,bootstrap_mean,bootstrap_p975,near_boundary_fraction\n');
for k=1:6
    fprintf(fid,'%s,%.12g,%.12g,%.12g,%.12g,%.12g\n',r.data.names{k}, ...
        r.base.parameters(k),r.checks.percentiles(1,k),r.checks.bootstrapMean(k), ...
        r.checks.percentiles(3,k),r.checks.boundaryFraction(k));
end
fclose(fid);
fid=fopen(fullfile(output,'historical_fit.csv'),'w');assert(fid>=0);
ok=find(r.completed & r.converged);ny=numel(r.data.years);
lower=NaN(ny,1);average=NaN(ny,1);upper=NaN(ny,1);
if numel(ok)>=2
    curves=zeros(ny,numel(ok));
    for j=1:numel(ok),curves(:,j)=r.fits{ok(j)}.predicted;end
    summary=india10_sample_summary(curves);
    lower=summary.lower;average=summary.mean;upper=summary.upper;
end
fprintf(fid,'year,observed,point_predicted,bootstrap_generation_sd,bootstrap_p025,bootstrap_mean,bootstrap_p975,n_bootstrap_in_summary\n');
fprintf(fid,'%d,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%d\n', ...
    [r.data.years r.data.cases r.base.predicted .2*r.base.predicted lower average upper numel(ok)*ones(ny,1)]');
fclose(fid);
fid=fopen(fullfile(output,'bootstrap_parameters.csv'),'w');assert(fid>=0);
fprintf(fid,'replicate,converged,alpha_1_1,alpha_1_2,alpha_1_3,alpha_1_4,rho,kappa\n');
for b=1:r.nBootstrap,fprintf(fid,'%d,%d,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n',b,r.converged(b),r.parameters(b,:));end
fclose(fid);
fid=fopen(fullfile(output,'bootstrap_fits.csv'),'w');assert(fid>=0);
fprintf(fid,'replicate,year,synthetic_observation,fitted_mean\n');
for b=1:r.nBootstrap
    fprintf(fid,'%d,%d,%.12g,%.12g\n',[b*ones(24,1) r.data.years r.observations(:,b) r.fits{b}.predicted]');
end
fclose(fid);
fid=fopen(fullfile(output,'initial_conditions.csv'),'w');assert(fid>=0);
fprintf(fid,'age_index,population,M,S,F,I,R,V1,V2,VSIA\n');
for a=1:10,fprintf(fid,'%d,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n',a,r.data.initialPopulation(a),1e6*r.data.initialState(:,a));end
fclose(fid);
fid=fopen(fullfile(output,'REPORT.txt'),'w');assert(fid>=0);
fprintf(fid,'INDIA TRC MANUSCRIPT: SIX-PARAMETER NONLINEAR LEAST SQUARES\n');
fprintf(fid,'Four transmission rates (without the legacy 0.3686 multiplier), rho and kappa estimated.\n');
fprintf(fid,'Unweighted annual residuals; Normal bootstrap SD=0.20*generating mean; no clipping.\n');
fprintf(fid,'Ten ages, historical simulation 1970-2023, observed annual incidence 2000-2023.\n');
fprintf(fid,'Maternal waning M_i -> S_i at 1/6 per biweek in every age group; ordinary ageing in all compartments.\n');
fprintf(fid,'Mode %s; point SSE %.12g; point convergence %d; agreeing starts %d.\n', ...
    r.options.mode,r.base.sse,r.base.converged,r.base.agreeingStarts);
fprintf(fid,'Bootstrap accepted fits %d/%d; all attempted fits retained in calibration.mat.\n',sum(r.converged),r.nBootstrap);
fprintf(fid,'CSV central summaries use arithmetic bootstrap means; 95%% intervals remain empirical percentiles.\n');
fprintf(fid,'Summary sample: completed, converged refits only; no forecast proceeds unless all refits pass.\n');
fprintf(fid,'Local log-sensitivity rank %d/6; condition number %.8g.\n',r.checks.rank,r.checks.conditionNumber);
fprintf(fid,'Max absolute bootstrap log-parameter correlation %.8g.\n',r.checks.maximumCorrelation);
fprintf(fid,'Standardized residual RMS %.8g; conditional bootstrap tail fraction %.8g.\n', ...
    r.checks.standardizedResidualRMS,r.checks.fitTailFraction);
fprintf(fid,'Practical screen passed %d; fit adequacy warning %d.\n',r.checks.practicalScreenPassed,r.checks.fitAdequacyWarning);
fprintf(fid,'Fit adequacy assessed with >=100 successful replicates: %d. A warning of 0 is not a pass when unassessed.\n',r.checks.fitAdequacyAssessed);
if r.nBootstrap<100
    fprintf(fid,'Fewer than 100 replicates: exported percentiles/correlations are execution diagnostics only, not inferential results.\n');
end
fprintf(fid,'scenarioGenerationAllowed=%d. No scenarios are run by this calibration command.\n',r.checks.scenarioGenerationAllowed);
fprintf(fid,'These are finite-sample practical diagnostics, not proof of structural identifiability.\n');
fprintf(fid,'Inspect rank, correlations, boundary frequencies and interval widths before scenario use.\n');
fprintf(fid,'Numerics: refinement %.8g; backend difference %.8g; mass error %.8g people.\n', ...
    r.checks.maximumRelativeStepError,r.checks.maximumRelativeBackendError,r.checks.massBalanceErrorPeople);
fprintf(fid,'Previous manuscript parameter tables and figures are not regenerated or validated by code fixes.\n');
fprintf(fid,'MATLAB %s; started %s; finished %s.\n',r.matlabVersion,r.started,r.finished);
fclose(fid);
end
