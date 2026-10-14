function [dx, incidence] = measles_model(x, phase, alpha, maternal, rates, data)
%MEASLES_MODEL Conservative flows, all evaluated at the same current state.
% rates = [birth, death, MCV1, MCV2, SIA]; units are inverse biweeks.
population = sum(x,1);
lambda = alpha*sum(x(4,:)./population)*(1-0.05*cos(2*pi*phase/26));
mInfection = maternal.*x(1,:).*lambda;
sInfection = x(2,:).*lambda;
fInfection = x(3,:).*lambda;
incidence = mInfection+sInfection+fInfection;
% Ageing includes infected people and cancels in the total population.
outAge = bsxfun(@times,x,data.omega);
outAge(1,:) = x(1,:).*data.maternalOmega;
dx = [zeros(8,1) outAge(:,1:9)]-outAge-rates(2)*x;
dx(1,1) = dx(1,1)+rates(1)*sum(population);
dx(1,:) = dx(1,:)-mInfection;
dx(2,:) = dx(2,:)-sInfection;
dx(3,:) = dx(3,:)-fInfection;
dx(4,:) = dx(4,:)+incidence-x(4,:);
dx(5,:) = dx(5,:)+x(4,:);
if isfield(data,'sameAgeWaning') && data.sameAgeWaning
    % Manuscript: immunity loss is M_i -> S_i, separate from ageing.
    wane = data.waning*x(1,:);
    dx(1,:) = dx(1,:)-wane;
    dx(2,:) = dx(2,:)+wane;
else
    wane = data.waning*x(1,3);
    dx(1,3) = dx(1,3)-wane;
    dx(2,4) = dx(2,4)+wane;
end
mcv1M = rates(3)*x(1,4);
mcv1S = rates(3)*x(2,4);
dx(1,4) = dx(1,4)-mcv1M;
dx(2,4) = dx(2,4)-mcv1S;
dx(3,4) = dx(3,4)+0.05*(mcv1M+mcv1S);
dx(6,4) = dx(6,4)+0.95*(mcv1M+mcv1S);
mcv2 = 0.98*rates(4)*x(3,5);
dx(3,5) = dx(3,5)-mcv2;
dx(7,5) = dx(7,5)+mcv2;
sia = 0.70*rates(5)*x(1:3,4:9);
dx(1:3,4:9) = dx(1:3,4:9)-sia;
dx(8,4:9) = dx(8,4:9)+sum(sia,1);
end
