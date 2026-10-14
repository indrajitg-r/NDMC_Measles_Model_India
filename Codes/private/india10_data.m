function data = india10_data()
%INDIA10_DATA Historical inputs copied from the supplied NDMC ten-age model.
% All coverage arrays refer to 1970:2023. Cases refer to 2000:2023.
data.modelVersion = 'india10-corrected-2026-10-12-v1';
data.years = (2000:2023)';
data.historyYears = (1970:2023)';
data.cases = flipud([65150 40967 5700 5604 10430 19474 12032 17250 ...
    30168 26530 8285 3305 33634 31458 56188 44258 41144 64185 ...
    36711 55443 47147 40044 51780 38835]');
data.names = {'alpha_1_1','alpha_1_2','alpha_1_3','alpha_1_4','rho','kappa'};
data.lower = [0.3 0.001 0.01 0.001 0.0025 0.001];
data.upper = [2.4 0.75 0.25 0.05 1 1];
% Original start had alpha_1_4=0.0009 below its lower bound; corrected here.
data.start = [0.9094 0.279 0.09 0.0012 0.0045 0.9];
data.relativeSD = 0.20;
data.mcv1 = [zeros(1,10) 0 0 0 0 0 1 10 24 32 42 56 43 51 59 ...
    67 72 66 55 53 56 56 57 56 60 64 68 69 70 72 78 82 84 83 ...
    83 85 87 88 90 93 95 89 89 95 93]/(100*6);
data.mcv2 = [zeros(1,30) 0 0 0 0 0 0 0 0 0 0 0 27 36 51 60 ...
    69 76 80 82 84 81 82 90 90]/(100*26);
data.siaSteps = 30*26 + [1 26 10*26+11*2 11*26+4*2 11*26+12*2 ...
    12*26+3*2 12*26+12*2 13*26+6*2 17*26+2*2 17*26+4*2 ...
    17*26+8*2 17*26+9*2 17*26+12*2 18*26+2*2 18*26+3*2 ...
    18*26+5*2 18*26+6*2 18*26+7*2 18*26+9*2 18*26+10*2 ...
    18*26+11*2 18*26+12*2 19*26+2*2];
data.siaCoverage = [0.7 0.9 3.6 1 10.7 3.5 15.7 14.4 4.3 4.6 ...
    3 3 2.5 0.1 3 0.3 0 3.8 3.9 5.4 3.4 27.2 16.5]/100;
data.sia = zeros(1,54*26);
data.sia(data.siaSteps) = data.siaCoverage;
data.omega = [ones(1,4)/6 ones(1,4)/26 1/260 0];
% Preserve the original convention: M1 -> M2 -> M3 -> S4.
data.maternalOmega = [1/6 1/6 zeros(1,8)];
data.waning = 0.1667;
data.initialPopulation = 545864267 * [15.8/20*ones(1,4) ...
    15.8/5*ones(1,4) 25.5 58.7]/100;
x = zeros(8,10); % M,S,F,I,R,V1,V2,VSIA; counts in millions internally.
x(1,1:3) = data.initialPopulation(1:3);
ageWeights = [0 0 0 962.58 610.71 68.22 68.22 68.22 57.72 14.93];
x(4,1:5) = 1;
x(1,1:3) = x(1,1:3)-x(4,1:3); % Allocate seed infections without R<0.
% Preserve the supplied total S and relative age weights subject to capacity.
% Each susceptible requires space for the original F=0.01*S as well.
% The user selected this allocation rather than 20% within each age group.
data.initialSusceptibleTotal = 0.2*sum(data.initialPopulation(4:9));
capacity = max(0,(data.initialPopulation-sum(x,1))/1.01);
available = ageWeights>0;
remaining = data.initialSusceptibleTotal;
while any(available)
    allocation = remaining*ageWeights(available)/sum(ageWeights(available));
    indices = find(available);
    capped = allocation>capacity(available);
    if ~any(capped)
        x(2,indices) = allocation;
        break
    end
    saturated = indices(capped);
    x(2,saturated) = capacity(saturated);
    remaining = remaining-sum(capacity(saturated));
    available(saturated) = false;
end
assert(abs(sum(x(2,:))-data.initialSusceptibleTotal)<1e-6);
x(3,:) = 0.01*x(2,:);
x(5,:) = data.initialPopulation-sum(x,1);
% Roundoff at exactly saturated age groups must not create negative counts.
x(5,abs(x(5,:))<1e-6) = 0;
assert(all(x(:)>=0), 'Initial compartments must be nonnegative.');
data.initialState = x/1e6;
assert(numel(data.mcv1)==54 && numel(data.mcv2)==54);
end
