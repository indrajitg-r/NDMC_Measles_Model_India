function result = india10_scenario_report(result,output)
%INDIA10_SCENARIO_REPORT Arithmetic means, percentile intervals and elimination.
years=result.settings.years;ny=numel(years);
n=size(result.annualTrue,3);ns=size(result.annualTrue,2);
assert(n>=2 && ns==12 && isequal(size(result.annualTrue),size(result.annualReported)));
threshold=result.settings.eliminationThresholdTrueCases;
result.eliminationYear=NaN(ns,n);result.eliminationYearFromMean=NaN(ns,1);
meanYear=NaN(ns,1);loYear=NaN(ns,1);hiYear=NaN(ns,1);
for j=1:ns
    cases=reshape(result.annualTrue(:,j,:),ny,n);
    result.eliminationYear(j,:)=calculate_elimination_years(years,cases,threshold);
    result.eliminationYearFromMean(j)=calculate_elimination_years(years,mean(cases,2),threshold);
    reached=result.eliminationYear(j,isfinite(result.eliminationYear(j,:)));
    if ~isempty(reached)
        meanYear(j)=mean(reached);
        if numel(reached)==1
            loYear(j)=reached;hiYear(j)=reached;
        else
            q=india10_sample_summary(reached);loYear(j)=q.lower;hiYear(j)=q.upper;
        end
    end
end
result.eliminationReached=isfinite(result.eliminationYear);
result.eliminationFraction=mean(result.eliminationReached,2);
result.meanEliminationYearAmongReached=meanYear;
result.eliminationIntervalAmongReached=[loYear hiYear];
result.eliminationInterpretation=sprintf(['Strictly <%.12g true cases nationally in each full year from the first qualifying year ' ...
    'through %d. NaN = not reached by the horizon. No claim about resurgence after the horizon.'],threshold,years(end));
save(fullfile(output,'scenarios.mat'),'result','-v7');

fid=fopen(fullfile(output,'annual_scenarios.csv'),'w');assert(fid>=0);
fprintf(fid,'scenario,year,p025_reported,mean_reported,p975_reported,p025_true,mean_true,p975_true\n');
for j=1:ns
    r=india10_sample_summary(reshape(result.annualReported(:,j,:),ny,n));
    t=india10_sample_summary(reshape(result.annualTrue(:,j,:),ny,n));
    fprintf(fid,'%d,%d,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n', ...
        [(j-1)*ones(ny,1) years r.lower r.mean r.upper t.lower t.mean t.upper]');
end
fclose(fid);
s=india10_sample_summary(result.percentAverted);
fid=fopen(fullfile(output,'percent_cases_averted.csv'),'w');assert(fid>=0);
fprintf(fid,'scenario,p025_percent,mean_percent,p975_percent\n');
fprintf(fid,'%d,%.12g,%.12g,%.12g\n',[(0:ns-1)' s.lower s.mean s.upper]');fclose(fid);

fid=fopen(fullfile(output,'elimination_years.csv'),'w');assert(fid>=0);
fprintf(fid,['scenario,year_from_mean_true_trajectory,fraction_reached_by_horizon,' ...
    'n_reached,n_not_reached,mean_year_among_reached,p025_year_among_reached,' ...
    'p975_year_among_reached,strict_true_case_threshold,last_forecast_year\n']);
for j=1:ns
    nr=sum(result.eliminationReached(j,:));
    fprintf(fid,'%d,%.12g,%.12g,%d,%d,%.12g,%.12g,%.12g,%.12g,%d\n', ...
        j-1,result.eliminationYearFromMean(j),result.eliminationFraction(j),nr,n-nr, ...
        meanYear(j),loYear(j),hiYear(j),threshold,years(end));
end
fclose(fid);
fid=fopen(fullfile(output,'bootstrap_elimination_years.csv'),'w');assert(fid>=0);
fprintf(fid,'scenario,replicate,elimination_year,reached_by_horizon\n');
for j=1:ns
    fprintf(fid,'%d,%d,%.12g,%d\n',[(j-1)*ones(n,1) (1:n)' ...
        result.eliminationYear(j,:)' result.eliminationReached(j,:)']');
end
fclose(fid);

s=india10_sample_summary(result.cumulativeUnassignedV1LossPeople);
fid=fopen(fullfile(output,'manuscript_equation_balance.csv'),'w');assert(fid>=0);
fprintf(fid,'scenario,p025_unassigned_V1_loss_people,mean_unassigned_V1_loss_people,p975_unassigned_V1_loss_people,max_numerical_balance_residual_people\n');
fprintf(fid,'%d,%.12g,%.12g,%.12g,%.12g\n',[(0:ns-1)' s.lower s.mean s.upper ...
    max(result.equationBalanceErrorPeople,[],2)]');fclose(fid);
fid=fopen(fullfile(output,'SCENARIO_NOTES.txt'),'w');assert(fid>=0);
fprintf(fid,'Literal printed MCV0 equations: N excludes V0; V1 loss has no V2 receiving term.\n');
fprintf(fid,'Unassigned V1 loss is reported separately from numerical integration residuals and natural deaths.\n');
fprintf(fid,'Annual cases and cases-averted summaries are arithmetic means with percentile intervals.\n');
fprintf(fid,'%s\n',result.eliminationInterpretation);
fprintf(fid,'Elimination-year means/intervals are conditional on reaching the criterion; consult the fraction and count not reached.\n');
fprintf(fid,'Elimination of the mean trajectory differs from the mean of individual elimination years.\n');
fprintf(fid,'These are conditional deterministic projections, not evidence of formal elimination.\n');
fclose(fid);
end
