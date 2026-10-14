function rates = india10_scenario_rates(step,scenario,settings)
%INDIA10_SCENARIO_RATES Calendar aligned future rates, per biweek.
yearIndex=floor((step-1)/26); year=settings.years(1)+yearIndex;
rates.birth=(40-(24/52)*(year-1970))/(1000*26);
rates.death=(17-(8/52)*(year-1970))/(1000*26);
rates.mcv1=settings.mcv1Coverage/6;rates.mcv2=settings.mcv2Coverage/26;
rates.sia=0;rates.mcv0=0;
if mod(step-1,26)==0
    if scenario.siaInterval>0 && mod(yearIndex,scenario.siaInterval)==0
        rates.sia=scenario.siaCoverage;
    end
    if scenario.mcv0Interval>0 && mod(yearIndex,scenario.mcv0Interval)==0
        rates.mcv0=scenario.mcv0Coverage;
    end
end
assert(rates.birth>=0 && rates.death>=0,'Demographic extrapolation became negative.');
end
