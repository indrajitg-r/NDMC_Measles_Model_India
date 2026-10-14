function data = india10_manuscript_data()
%INDIA10_MANUSCRIPT_DATA Six-parameter version of the TRC methods equations.
% Parameter order is alpha_1_1, alpha_1_2, alpha_1_3, alpha_1_4, rho, kappa.
data=india10_data();
data.modelVersion='india10-manuscript-lsq-2026-10-13-v2';
data.transmissionScale=1; % alpha is the transmission rate itself, per biweek.
data.sameAgeWaning=true;  % M_i -> S_i in the same age group, as in the equations.
data.maternalOmega=data.omega;
data.waning=1/6;         % Table 1 gives one rate; apply it in each age group.
% Positive numerical floor for log coordinates; upper limits follow the table.
data.lower=1e-8*ones(1,6);
data.upper=[5 5 5 5 1 1];
data.start=[0.8249 0.5971 0.0658 0.0012 0.0045 0.9];
data.method='unweighted nonlinear least squares (lsqcurvefit)';
end
