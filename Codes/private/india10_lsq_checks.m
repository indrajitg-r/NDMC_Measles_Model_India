function checks = india10_lsq_checks(result)
%INDIA10_LSQ_CHECKS Numerical, local-sensitivity and bootstrap warning checks.
% These are practical screens, not a structural-identifiability proof.
d=result.data;p=result.base.parameters;h=result.options.substeps;
[mu,state]=india10_simulate(p,d,h);fine=india10_simulate(p,d,2*h);
checks.maximumRelativeStepError=max(abs(mu-fine)./max(1,fine));
checks.maximumRelativeBackendError=max(abs(mu-result.base.predicted)./max(1,mu));
checks.massBalanceErrorPeople=state.maxMassBalanceErrorPeople;
checks.minimumCompartmentPeople=state.minCompartmentPeople;
checks.numericsPassed=checks.maximumRelativeStepError<.01 && ...
    checks.maximumRelativeBackendError<1e-7 && checks.massBalanceErrorPeople<1e-3;
J=zeros(24,6);delta=1e-4;
for k=1:6
    lo=p;hi=p;lo(k)=max(d.lower(k),p(k)*exp(-delta));
    hi(k)=min(d.upper(k),p(k)*exp(delta));
    J(:,k)=(log(india10_predict(hi,d,h))-log(india10_predict(lo,d,h)))/log(hi(k)/lo(k));
end
checks.logSensitivity=J;s=svd(J,0);checks.singularValues=s;
checks.rank=sum(s>1e-6*s(1));checks.conditionNumber=s(1)/max(s(end),realmin);
q=(log(p)-log(d.lower))./(log(d.upper)-log(d.lower));
checks.baseNearBoundary=q<=.01 | q>=.99;
checks.standardizedResidual=(d.cases-mu)./(.20*mu);
checks.standardizedResidualRMS=sqrt(mean(checks.standardizedResidual.^2));
ok=result.completed & result.converged;
checks.successfulFits=sum(ok);checks.totalFits=result.nBootstrap;
checks.bootstrapCorrelation=NaN(6);checks.boundaryFraction=NaN(1,6);
checks.percentiles=NaN(3,6);checks.maximumCorrelation=NaN;
checks.bootstrapMean=NaN(1,6);
checks.relativeIntervalWidth=NaN(1,6);checks.fitTailFraction=NaN;
if sum(ok)>=2
    values=result.parameters(ok,:);checks.bootstrapCorrelation=corrcoef(log(values));
    c=abs(checks.bootstrapCorrelation);c(1:7:end)=0;checks.maximumCorrelation=max(c(:));
    q=bsxfun(@rdivide,bsxfun(@minus,log(values),log(d.lower)),log(d.upper)-log(d.lower));
    checks.boundaryFraction=mean(q<=.01 | q>=.99,1);
    summary=india10_sample_summary(values');
    checks.percentiles=[summary.lower';summary.median';summary.upper'];
    checks.bootstrapMean=summary.mean';
    checks.relativeIntervalWidth=(summary.upper'-summary.lower')./p;
    rms=zeros(sum(ok),1);indices=find(ok);
    for j=1:numel(indices)
        b=indices(j);pred=result.fits{b}.predicted;
        rms(j)=sqrt(mean(((result.observations(:,b)-pred)./(.20*pred)).^2));
    end
    checks.fitTailFraction=(1+sum(rms>=checks.standardizedResidualRMS))/(1+numel(rms));
end
checks.practicalScreenPassed=result.nBootstrap>=100 && all(ok) && ...
    result.base.converged && result.base.agreeingStarts>=2 && checks.numericsPassed && ...
    checks.rank==6 && checks.conditionNumber<1e6 && ~any(checks.baseNearBoundary) && ...
    isfinite(checks.maximumCorrelation) && checks.maximumCorrelation<.95 && ...
    all(checks.boundaryFraction<=.10) && all(checks.relativeIntervalWidth<=1);
checks.fitAdequacyAssessed=sum(ok)>=100;
checks.fitAdequacyWarning=checks.fitAdequacyAssessed && ...
    isfinite(checks.fitTailFraction) && checks.fitTailFraction<.05;
checks.scenarioGenerationAllowed=checks.practicalScreenPassed && ~checks.fitAdequacyWarning;
end
