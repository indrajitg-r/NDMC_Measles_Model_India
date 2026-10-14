function [eliminationYear,reached] = calculate_elimination_years(years,annualTrueCases,threshold)
%CALCULATE_ELIMINATION_YEARS First sustained year below a true-case threshold.
% Rows are full consecutive calendar years; columns are trajectories.
% Default: <1 true case nationally per year, with no subsequent resurgence
% within the supplied horizon. NaN means not reached within that horizon.
% This is a deterministic model criterion, not certification of elimination.
if nargin<3,threshold=1;end
years=years(:);
assert(~isempty(years) && all(isfinite(years)) && all(years==round(years)) && all(diff(years)==1), ...
    'india10:InvalidYears','Supply consecutive full calendar years.');
assert(isnumeric(annualTrueCases) && ismatrix(annualTrueCases) && size(annualTrueCases,1)==numel(years) && ...
    size(annualTrueCases,2)>=1 && all(isfinite(annualTrueCases(:))) && all(annualTrueCases(:)>=0), ...
    'india10:InvalidIncidence','Supply finite nonnegative annual TRUE case counts, one year per row.');
assert(isscalar(threshold) && isfinite(threshold) && threshold>0, ...
    'india10:InvalidThreshold','The strict case threshold must be positive.');
staysBelow=flipud(cumprod(double(flipud(annualTrueCases<threshold)),1))>0;
eliminationYear=NaN(1,size(annualTrueCases,2));
for j=1:size(annualTrueCases,2)
    first=find(staysBelow(:,j),1,'first');
    if ~isempty(first),eliminationYear(j)=years(first);end
end
reached=isfinite(eliminationYear);
end
