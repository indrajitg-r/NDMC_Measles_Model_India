function summary = india10_sample_summary(values)
%INDIA10_SAMPLE_SUMMARY Rows are quantities; columns are paired replicates.
assert(ismatrix(values) && size(values,2)>=2 && all(isfinite(values(:))));
n=size(values,2);sorted=sort(values,2);index=1+(n-1)*[.025 .5 .975];
q=zeros(size(values,1),3);
for k=1:3
    lo=floor(index(k));hi=ceil(index(k));w=index(k)-lo;
    q(:,k)=(1-w)*sorted(:,lo)+w*sorted(:,hi);
end
summary.lower=q(:,1);summary.median=q(:,2);summary.upper=q(:,3);
summary.mean=mean(values,2);summary.n=n;
end
