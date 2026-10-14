function fit = india10_lsq_fit(data,options,starts)
%INDIA10_LSQ_FIT Fit predicted annual cases, not squared residuals, to data.
% A common scale improves conditioning without changing the OLS minimizer.
scale=max(1,sqrt(mean(data.cases.^2)));
settings=optimset('Algorithm','trust-region-reflective','Display','off', ...
    'MaxIter',options.maxIterations,'MaxFunEvals',options.maxEvaluations, ...
    'TolFun',1e-9,'TolX',1e-8,'FinDiffType','central');
records=struct('parameters',{},'sse',{},'exitflag',{},'predicted',{},'output',{});
for j=1:size(starts,1)
    seed=min(max(starts(j,:),data.lower),data.upper);
    [z,~,~,flag,out]=lsqcurvefit(@prediction,log(seed),data.years, ...
        data.cases/scale,log(data.lower),log(data.upper),settings);
    mu=prediction(z,data.years)*scale;
    records(j).parameters=exp(z);records(j).sse=sum((mu-data.cases).^2);
    records(j).exitflag=flag;records(j).predicted=mu;records(j).output=out;
    fprintf('OLS start %d/%d: SSE %.8g, exitflag %d.\n',j,size(starts,1),records(j).sse,flag);
end
[~,best]=min([records.sse]);fit=records(best);fit.starts=records;
fit.converged=fit.exitflag>0 && isfinite(fit.sse) && all(fit.predicted<1e25);
fit.agreeingStarts=sum([records.exitflag]>0 & ...
    abs([records.sse]-fit.sse)<=1e-4*max(1,fit.sse));
    function mu=prediction(z,~)
        try
            mu=india10_predict(exp(z),data,options.substeps)/scale;
        catch err
            if ~any(strcmp(err.identifier,{'india10:InvalidState','india10:InvalidMean'}))
                rethrow(err)
            end
            mu=1e30*ones(size(data.cases));
        end
    end
end
