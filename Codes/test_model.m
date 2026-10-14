function varargout = test_model(backend)
%TEST_MODEL Check the equation changes without a production refit.
if nargin<1,backend='matlab';end
d=india10_manuscript_data();p=d.start;
assert(d.transmissionScale==1 && d.waning==1/6 && isequal(d.maternalOmega,d.omega));
assert(isequal(d.upper,[5 5 5 5 1 1]) && d.relativeSD==.20);
assert(all(d.initialState(:)>=0));
assert(max(abs(1e6*sum(d.initialState,1)-d.initialPopulation))<1e-6);
assert(abs(1e6*sum(d.initialState(2,:))-d.initialSusceptibleTotal)<1e-6);
% An isolated M cohort loses immunity into S of the SAME age and ages into M.
x=zeros(8,10);x(5,:)=1;x(1,3)=2;
dx=measles_model(x,0,zeros(1,10),zeros(1,10),zeros(1,5),d);
assert(abs(dx(2,3)-2/6)<1e-12 && dx(2,4)==0);
assert(abs(dx(1,4)-2*d.omega(3))<1e-12 && abs(sum(dx(:)))<1e-12);
% Direct force of infection: four coefficients enter without an extra scale.
x=reshape(1:80,8,10)/100;alpha=[p(1)*ones(1,4) p(2)*ones(1,4) p(3:4)];
maternal=p(6)*[.2 .6 1 zeros(1,7)];rates=[.001 .0003 .1 .03 .2];
[dx,inc]=measles_model(x,3,alpha,maternal,rates,d);
expected=alpha*sum(x(4,:)./sum(x,1))*(1-.05*cos(2*pi*3/26)).* ...
    (maternal.*x(1,:)+x(2,:)+x(3,:));
assert(max(abs(inc-expected))<1e-12);
assert(abs(sum(dx(:))-(rates(1)-rates(2))*sum(x(:)))<1e-12);
[mu,state]=india10_simulate(p,d,4);fine=india10_simulate(p,d,8);
assert(all(mu>0) && isequal(size(mu),[24 1]));
assert(max(abs(mu-fine)./fine)<.01 && state.maxMassBalanceErrorPeople<1e-3);
for k=1:6
    q=p;q(k)=p(k)*.9;
    assert(norm(india10_simulate(q,d,4)-mu)/norm(mu)>1e-7,'Inactive parameter %d',k);
end
% No intervention change: historical and scenario RHS must agree exactly.
settings=india10_scenario_settings();s=settings.scenarios(1);
r.birth=rates(1);r.death=rates(2);r.mcv1=rates(3);r.mcv2=rates(4);r.sia=rates(5);r.mcv0=0;
[future,finc]=measles_scenario_model([x;zeros(1,10)],3,p,r,s,settings,d);
assert(max(max(abs(future(1:8,:)-dx)))<1e-12 && max(abs(finc-inc))<1e-12);
% Literal extra equations: no V0 -> V1, no V2 receipt of the V1 outflow.
s=settings.scenarios(12);r.mcv0=.28;x9=reshape(1:90,9,10)/100;
[future,~,population,loss]=measles_scenario_model(x9,3,p,r,s,settings,d);
expected=measles_model(x9(1:8,:),3,alpha,maternal,rates,d);
v0age=x9(9,:).*d.omega;
expected(9,:)=[0 v0age(1:9)]-v0age-r.death*x9(9,:);
zero=.70*r.mcv0*x9(1,3);unassigned=.98*r.mcv2*x9(6,5);
expected(1,3)=expected(1,3)-zero;expected(9,3)=expected(9,3)+zero;
expected(6,5)=expected(6,5)-unassigned;
assert(max(abs(future(:)-expected(:)))<1e-12 && abs(loss-unassigned)<1e-12);
assert(max(abs(population-sum(x9(1:8,:),1)))<1e-12);
assert(abs(sum(future(:))-(r.birth*sum(population)-r.death*sum(x9(:))-loss))<1e-12);
for j=1:12
    s=settings.scenarios(j);
    projection=india10_scenario_simulate(p,state.finalStatePeople/1e6,s,settings,d,4);
    assert(all(projection.trueCases(:)>=0) && projection.maxEquationBalanceErrorPeople<1e-3);
    assert(max(abs(projection.reportedCases(:)-p(5)*projection.trueCases(:)))<1e-8);
    if j<=8
        assert(projection.cumulativeUnassignedV1LossPeople==0);
    else
        assert(projection.cumulativeUnassignedV1LossPeople>0);
    end
end
% Complete years; transient dips, equality at 1 and non-attainment are explicit.
years=(2024:2027)';counts=[.2 .2 2 .2 2; .3 2 1 .3 3; .4 .4 .2 .4 4; .5 .5 .3 2 .4];
[ey,reached]=calculate_elimination_years(years,counts);
assert(isequaln(ey,[2024 2026 2026 NaN 2027]) && isequal(reached,[true true true false true]));
assert(isnan(calculate_elimination_years(years,ones(4,1))));
assert(calculate_elimination_years(years,zeros(4,1))==2024);
% Synthetic exact data exercise lsqcurvefit's prediction/target convention.
synthetic=d;synthetic.cases=mu;options.substeps=4;options.maxIterations=30;options.maxEvaluations=300;
fit=india10_lsq_fit(synthetic,options,p);
assert(fit.converged && fit.sse<1e-8 && max(abs(fit.predicted-mu))<1e-4);
[~,validation]=india10_backend(d,backend,fullfile(fileparts(mfilename('fullpath')),'outputs','test_model_build'));
if nargout>0
    info.predicted=mu;info.parameters=p;info.backendValidation=validation;
    info.solver=functions(@india10_simulate);info.fitter=functions(@india10_lsq_fit);
    varargout{1}=info;
end
fprintf('IMPLEMENTATION TESTS PASSED: literal equations, six active parameters, OLS, 12 scenarios and elimination logic.\n');
fprintf('Printed MCV0 equations have unassigned V1 loss; its magnitude is recorded separately.\n');
end
