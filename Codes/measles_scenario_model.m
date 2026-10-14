function [dx,incidence,population,unassignedLoss] = measles_scenario_model(x,phase,alpha,rates,scenario,settings,data)
%MEASLES_SCENARIO_MODEL Ten ages; M,S,F,I,R,V1,V2,VSIA,V0 in rows 1--9.
% Literal printed manuscript: N excludes V0; V1 has an unmatched MCV2 loss.
% unassignedLoss reports that loss explicitly; it is not a natural death.
population=sum(x(1:8,:),1);
if numel(alpha)==6
    p=alpha; % Manuscript parameter order: four alpha values, rho, kappa.
    alpha=[p(1)*ones(1,4) p(2)*ones(1,4) p(3) p(4)];
    maternal=p(6)*[0.2 0.6 1 zeros(1,7)];
else
    alpha=0.3686*alpha;maternal=[0.2 0.6 1 zeros(1,7)];
end
lambda=alpha*sum(x(4,:)./population)*(1-0.05*cos(2*pi*phase/26));
mi=maternal.*x(1,:).*lambda;si=x(2,:).*lambda;fi=x(3,:).*lambda;
incidence=mi+si+fi;
outAge=bsxfun(@times,x,data.omega);outAge(1,:)=x(1,:).*data.maternalOmega;
dx=[zeros(9,1) outAge(:,1:9)]-outAge-rates.death*x;
dx(1,1)=dx(1,1)+rates.birth*sum(population);
dx(1,:)=dx(1,:)-mi;dx(2,:)=dx(2,:)-si;dx(3,:)=dx(3,:)-fi;
dx(4,:)=dx(4,:)+incidence-x(4,:);dx(5,:)=dx(5,:)+x(4,:);
if isfield(data,'sameAgeWaning') && data.sameAgeWaning
    wane=data.waning*x(1,:);dx(1,:)=dx(1,:)-wane;dx(2,:)=dx(2,:)+wane;
else
    wane=data.waning*x(1,3);dx(1,3)=dx(1,3)-wane;dx(2,4)=dx(2,4)+wane;
end
a=scenario.mcv1Age;
first=rates.mcv1*x(1:2,a);
dx(1:2,a)=dx(1:2,a)-first;
dx(3,a)=dx(3,a)+(1-scenario.mcv1Efficacy)*sum(first);
dx(6,a)=dx(6,a)+scenario.mcv1Efficacy*sum(first);
second=settings.mcv2Efficacy*rates.mcv2*x(3,5);
dx(3,5)=dx(3,5)-second;dx(7,5)=dx(7,5)+second;
campaign=settings.siaEfficacy*rates.sia*x(1:3,settings.siaAges);
dx(1:3,settings.siaAges)=dx(1:3,settings.siaAges)-campaign;
dx(8,settings.siaAges)=dx(8,settings.siaAges)+sum(campaign,1);
early=settings.mcv0Efficacy*rates.mcv0*x(1,3);
dx(1,3)=dx(1,3)-early;dx(9,3)=dx(9,3)+early;
% Eq. for V0 has only MCV0 inflow, ageing and death: no V0 -> V1 transfer.
% The modified V1 equation subtracts this flow, but the printed V2 equation
% has no receiving term. Reproduce that omission, and expose its magnitude.
unassignedLoss=0;
if scenario.mcv0Interval>0
    unassignedLoss=settings.mcv2Efficacy*rates.mcv2*x(6,5);
    dx(6,5)=dx(6,5)-unassignedLoss;
end
end
