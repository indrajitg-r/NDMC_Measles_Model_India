function settings = india10_scenario_settings()
%INDIA10_SCENARIO_SETTINGS Manuscript scenarios 0--11 for India, 2024--2050.
settings.version='india10-literal-manuscript-scenarios-2026-10-13-v2';
settings.years=(2024:2050)'; settings.substeps=2;
settings.mcv1Coverage=0.93; settings.mcv2Coverage=0.90;
settings.standardMCV1Efficacy=0.95; settings.earlyMCV1Efficacy=0.70;
% Early MCV1 has an absolute reduction of 25 percentage points: .95-.25.
settings.mcv2Efficacy=0.98; settings.siaEfficacy=0.70;
settings.siaCoverage=0.165; settings.siaAges=4:9;
settings.mcv0Coverage=0.28; settings.mcv0Efficacy=0.70;
settings.campaignConvention='coverage used as a rate during one biweek; not an instantaneous jump';
settings.eliminationThresholdTrueCases=1; % Strictly below, nationally per year.
settings.eliminationCriterion='First full calendar year with <1 true case and no later year >=1 through 2050.';
settings.equationConvention='Literal printed MCV0 equations: N excludes V0; unmatched V1 outflow is recorded.';
settings.ageNames={'0-3 months','3-6 months','6-9 months','9-12 months', ...
    '1-2 years','2-3 years','3-4 years','4-5 years','5-15 years','15-80 years'};
intervals=[0 3 2 1 0 3 2 1 0 3 2 1];
for id=0:11
    s.id=id; s.mcv1Age=4; s.mcv1Efficacy=settings.standardMCV1Efficacy;
    if id>=4 && id<=7, s.mcv1Age=3;s.mcv1Efficacy=settings.earlyMCV1Efficacy;end
    s.siaInterval=intervals(id+1);s.siaCoverage=settings.siaCoverage;
    s.mcv0Interval=0;s.mcv0Coverage=0;
    if id>=8, s.mcv0Interval=3;s.mcv0Coverage=settings.mcv0Coverage;end
    s.name=sprintf('Scenario %d',id);
    settings.scenarios(id+1)=s;
end
end
