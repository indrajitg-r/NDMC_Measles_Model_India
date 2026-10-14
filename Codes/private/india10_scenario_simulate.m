function projection = india10_scenario_simulate(theta,initial,scenario,settings,data,substeps)
%INDIA10_SCENARIO_SIMULATE Reference RK4 projection from end-2023 state.
if nargin<6, substeps=settings.substeps;end
assert(any(numel(theta)==[2 6]) && all(isfinite(theta)) && all(theta>0));
theta=theta(:)';
transmission=theta(1);rho=theta(2);
if numel(theta)==6,transmission=theta;rho=theta(5);end
assert(isequal(size(initial),[8 10]) || isequal(size(initial),[9 10]));
assert(all(isfinite(initial(:))) && all(initial(:)>=-1e-10));
assert(substeps>=1 && substeps==round(substeps));
assert(settings.years(1)==data.historyYears(end)+1 && all(diff(settings.years)==1));
x=initial;if size(x,1)==8,x=[x;zeros(1,10)];end
ny=numel(settings.years);h=1/substeps;
annual=zeros(ny,10);personTime=zeros(ny,10);endPopulation=zeros(ny,10);
annualUnassignedLoss=zeros(ny,1);
minimum=min(x(:));balanceError=0;
for step=1:ny*26
    year=floor((step-1)/26)+1;
    rates=india10_scenario_rates(step,scenario,settings);
    for sub=1:substeps
        phase=numel(data.historyYears)*26+step+(sub-1)*h;
        [a,ia,pa,la]=measles_scenario_model(x,phase,transmission,rates,scenario,settings,data);
        xb=x+h*a/2;
        [b,ib,pb,lb]=measles_scenario_model(xb,phase+h/2,transmission,rates,scenario,settings,data);
        xc=x+h*b/2;
        [c,ic,pc,lc]=measles_scenario_model(xc,phase+h/2,transmission,rates,scenario,settings,data);
        xd=x+h*c;
        [d,id,pd,ld]=measles_scenario_model(xd,phase+h,transmission,rates,scenario,settings,data);
        next=x+h*(a+2*b+2*c+d)/6;
        if any(~isfinite(next(:))) || min(next(:)) < -1e-10
            error('india10:InvalidProjection','Nonfinite/negative forecast state.');
        end
        integratedPopulation=h*(pa+2*pb+2*pc+pd)/6;
        % N in births excludes V0, but natural deaths apply to all nine rows.
        integratedTotal=h*(sum(x(:))+2*sum(xb(:))+2*sum(xc(:))+sum(xd(:)))/6;
        omitted=h*(la+2*lb+2*lc+ld)/6;
        expected=rates.birth*sum(integratedPopulation)-rates.death*integratedTotal-omitted;
        balanceError=max(balanceError,abs(sum(next(:))-sum(x(:))-expected));
        annualUnassignedLoss(year)=annualUnassignedLoss(year)+omitted;
        annual(year,:)=annual(year,:)+h*(ia+2*ib+2*ic+id)/6;
        personTime(year,:)=personTime(year,:)+integratedPopulation/26;
        x=next;minimum=min(minimum,min(x(:)));
    end
    if mod(step,26)==0,endPopulation(year,:)=sum(x(1:8,:),1);end
end
projection.trueCases=1e6*annual;projection.reportedCases=rho*projection.trueCases;
projection.personYears=1e6*personTime;projection.endPopulation=1e6*endPopulation;
projection.finalStatePeople=1e6*x;projection.minimumCompartmentPeople=1e6*minimum;
projection.maxEquationBalanceErrorPeople=1e6*balanceError;
projection.annualUnassignedV1LossPeople=1e6*annualUnassignedLoss;
projection.cumulativeUnassignedV1LossPeople=1e6*sum(annualUnassignedLoss);
projection.populationDefinition='Printed N = sum of eight primary compartments; excludes V0.';
projection.substeps=substeps;
end
