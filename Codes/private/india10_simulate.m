function [mu, diagnostics] = india10_simulate(p, data, substeps)
%INDIA10_SIMULATE Historical ten-age simulation only: no future scenarios.
% Fourth-order Runge-Kutta, with exact annual/biweekly intervention boundaries.
if nargin<2 || isempty(data), data=india10_data(); end
if nargin<3, substeps=2; end
assert(numel(p)==6 && all(isfinite(p)) && all(p>0), 'Expected six positive parameters.');
assert(substeps>=1 && substeps==round(substeps), 'Substeps must be a positive integer.');
p=p(:)';
scale=0.3686; % Retain the archived multiplier convention by default.
if isfield(data,'transmissionScale'), scale=data.transmissionScale; end
alpha=scale*[p(1)*ones(1,4) p(2)*ones(1,4) p(3) p(4)];
maternal=p(6)*[0.2 0.6 1 zeros(1,7)];
x=data.initialState;
h=1/substeps;
annual=zeros(54,10);
minState=min(x(:));
massError=0;
for step=1:54*26
    year=floor((step-1)/26)+1;
    birth=(40-(24/52)*(year-1))/(1000*26);
    death=(17-(8/52)*(year-1))/(1000*26);
    rates=[birth death data.mcv1(year) data.mcv2(year) data.sia(step)];
    for sub=1:substeps
        phase=step+(sub-1)*h;
        [a,ia]=measles_model(x,phase,alpha,maternal,rates,data);
        xb=x+h*a/2;
        [b,ib]=measles_model(xb,phase+h/2,alpha,maternal,rates,data);
        xc=x+h*b/2;
        [c,ic]=measles_model(xc,phase+h/2,alpha,maternal,rates,data);
        xd=x+h*c;
        [d,id]=measles_model(xd,phase+h,alpha,maternal,rates,data);
        next=x+h*(a+2*b+2*c+d)/6;
        if any(~isfinite(next(:))) || min(next(:)) < -1e-10
            error('india10:InvalidState','Nonfinite/negative compartment: increase substeps.');
        end
        % Compare total change to births minus deaths integrated by the same rule.
        expected=(birth-death)*h*(sum(x(:))+2*sum(xb(:))+2*sum(xc(:))+sum(xd(:)))/6;
        massError=max(massError,abs(sum(next(:))-sum(x(:))-expected));
        annual(year,:)=annual(year,:)+h*(ia+2*ib+2*ic+id)/6;
        x=next;
        minState=min(minState,min(x(:)));
    end
end
mu=p(5)*1e6*sum(annual(31:54,:),2);
if any(~isfinite(mu)) || any(mu<=0)
    error('india10:InvalidMean','Predicted annual cases must be finite and positive.');
end
if nargout>1
    diagnostics.minCompartmentPeople=1e6*minState;
    diagnostics.maxMassBalanceErrorPeople=1e6*massError;
    diagnostics.annualTrueCases=1e6*annual(31:54,:);
    diagnostics.finalStatePeople=1e6*x;
    diagnostics.substeps=substeps;
end
end
