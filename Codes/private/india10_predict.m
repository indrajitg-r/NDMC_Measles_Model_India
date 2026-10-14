function mu = india10_predict(p,data,substeps)
%INDIA10_PREDICT Optional verified accelerator; MATLAB remains the reference.
if isfield(data,'backend') && strcmp(data.backend,'dotnet')
    if isfield(data,'sameAgeWaning') && data.sameAgeWaning
        mu=double(NDMC.India10Kernel.SimulateManuscript(p(:)',data.initialState(:)',data.mcv1, ...
            data.mcv2,data.sia,data.omega,data.maternalOmega,data.waning,int32(substeps)));
    else
        mu=double(NDMC.India10Kernel.Simulate(p(:)',data.initialState(:)',data.mcv1, ...
            data.mcv2,data.sia,data.omega,data.maternalOmega,data.waning,int32(substeps)));
    end
    mu=mu(:);
    if any(~isfinite(mu)) || any(mu<=0)
        error('india10:InvalidMean','Accelerated simulation returned invalid annual means.');
    end
else
    mu=india10_simulate(p,data,substeps);
end
end
