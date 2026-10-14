function [data,validation] = india10_backend(data,requested,buildFolder)
%INDIA10_BACKEND Build an optional Windows accelerator and verify equivalence.
% 'matlab' is portable. 'auto' falls back if .NET/compiler is unavailable.
assert(any(strcmp(requested,{'auto','matlab','dotnet'})));
data.backend='matlab';
validation.requested=requested; validation.backend='matlab';
validation.maximumRelativeDifference=0; validation.testedVectors=0;
if strcmp(requested,'matlab'), return; end
compiler=fullfile(getenv('WINDIR'),'Microsoft.NET','Framework64','v4.0.30319','csc.exe');
if ~ispc || exist(compiler,'file')~=2
    assert(~strcmp(requested,'dotnet'),'Windows .NET Framework compiler unavailable.');
    fprintf('Using reference MATLAB simulator; optional compiler unavailable.\n');
    return
end
if exist(buildFolder,'dir')~=7, mkdir(buildFolder); end
source=fullfile(fileparts(mfilename('fullpath')),'India10Kernel.cs');
dll=fullfile(buildFolder,'India10Kernel.dll');
digest=java.security.MessageDigest.getInstance('SHA-256');
digest.update(uint8(fileread(source)));
sourceHash=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2)',1,[]));
hashFile=fullfile(buildFolder,'India10KernelFingerprint.cs');
hashSource=sprintf(['namespace NDMC {public static partial class India10Kernel {' ...
    'public static string SourceHash(){return "%s";}}}'],sourceHash);
if exist(dll,'file')~=2 || exist(hashFile,'file')~=2 || ~strcmp(fileread(hashFile),hashSource)
    fid=fopen(hashFile,'w');assert(fid>=0);fprintf(fid,'%s',hashSource);fclose(fid);
    [status,message]=system(sprintf('"%s" /nologo /optimize+ /target:library /out:"%s" "%s" "%s"', ...
        compiler,dll,source,hashFile));
    if status~=0
        assert(~strcmp(requested,'dotnet'),'Accelerator build failed: %s',message);
        fprintf('Using reference MATLAB simulator; accelerator build failed.\n'); return
    end
end
try
    NET.addAssembly(dll);
catch err
    if strcmp(requested,'dotnet'), rethrow(err); end
    fprintf('Using reference MATLAB simulator; .NET load failed: %s\n',err.message); return
end
data.backend='dotnet';
assert(strcmp(char(NDMC.India10Kernel.SourceHash()),sourceHash), ...
    'A different accelerator is already loaded. Restart MATLAB after editing the C# source.');
validation.sourceHash=sourceHash;
% Use a private random stream: parity checks do not change analysis seeds.
stream=RandStream('mt19937ar','Seed',4182026);
u=rand(stream,12,6);
vectors=bsxfun(@plus,data.lower,bsxfun(@times,u,data.upper-data.lower));
vectors=[data.start;data.lower;data.upper;0.98156 0.39183 0.01 0.05 0.00523 1;vectors];
maximum=0;
bothInvalid=0;
for j=1:size(vectors,1)
    reference=[];accelerated=[];
    try
        reference=india10_simulate(vectors(j,:),data,2);
    catch err
        if ~any(strcmp(err.identifier,{'india10:InvalidMean','india10:InvalidState'})), rethrow(err); end
    end
    try
        accelerated=india10_predict(vectors(j,:),data,2);
    catch err
        if ~strcmp(err.identifier,'india10:InvalidMean'), rethrow(err); end
    end
    assert(isempty(reference)==isempty(accelerated),'Accelerator feasibility differs from reference.');
    if isempty(reference), bothInvalid=bothInvalid+1; continue; end
    maximum=max(maximum,max(abs(reference-accelerated)./reference));
end
% Verify nondefault time step and passed-in input arrays, not just defaults.
changed=data;changed.mcv1=0.91*data.mcv1;changed.mcv2=0.87*data.mcv2;
changed.sia=0.83*data.sia;changed.waning=0.95*data.waning;
changed.initialState(2:3,:)=0.9*data.initialState(2:3,:);
reference=india10_simulate(vectors(4,:),changed,4);
accelerated=india10_predict(vectors(4,:),changed,4);
maximum=max(maximum,max(abs(reference-accelerated)./reference));
assert(maximum<1e-7,'Accelerator differs from reference; aborting analysis.');
validation.backend='dotnet';validation.maximumRelativeDifference=maximum;
validation.testedVectors=size(vectors,1)+1;validation.assembly=dll;
validation.matchedInvalidVectors=bothInvalid;
fprintf('Accelerator verified against MATLAB: %d cases; maximum relative difference %.3g.\n', ...
    validation.testedVectors,maximum);
end
