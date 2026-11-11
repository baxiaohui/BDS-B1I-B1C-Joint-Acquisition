function cfg = resolveJointAcquisitionOptions(cfg)
% Separate branch durations and manual complex-amplitude weights.
% Legacy coherentMs-only calls retain their unnormalized [1 1] statistic.
hasI=isfield(cfg,'b1iCoherentMs') && ~isempty(cfg.b1iCoherentMs);
hasC=isfield(cfg,'b1cCoherentMs') && ~isempty(cfg.b1cCoherentMs);
legacy=~hasI && ~hasC && isfield(cfg,'coherentMs');
if ~hasI
    assert(isfield(cfg,'coherentMs'),'Set b1iCoherentMs.');
    cfg.b1iCoherentMs=cfg.coherentMs;
end
if ~hasC
    assert(isfield(cfg,'coherentMs'),'Set b1cCoherentMs.');
    cfg.b1cCoherentMs=cfg.coherentMs;
end
for field={'b1iCoherentMs','b1cCoherentMs'}
    value=cfg.(field{1});
    assert(isnumeric(value) && isreal(value) && isscalar(value) && isfinite(value) && value>0, ...
        '%s must be a positive finite duration in milliseconds.',field{1});
end
cfg.maxCoherentMs=max(cfg.b1iCoherentMs,cfg.b1cCoherentMs);
if ~isfield(cfg,'jointPhaseMode'), cfg.jointPhaseMode='continuous'; end
phaseMode=lower(string(cfg.jointPhaseMode));
assert(isscalar(phaseMode) && ismember(phaseMode,["continuous","grid"]), ...
    'jointPhaseMode must be continuous or grid.');
cfg.jointPhaseMode=char(phaseMode);
if ~isfield(cfg,'jointPhaseBins'), cfg.jointPhaseBins=360; end
bins=cfg.jointPhaseBins;
assert(isnumeric(bins) && isreal(bins) && isscalar(bins) && isfinite(bins) && ...
    bins>=1 && bins==fix(bins) && bins<=flintmax, ...
    'jointPhaseBins must be a positive integer (number of uniform phases over 360 degrees).');
cfg.jointPhaseBins=double(bins);
if phaseMode=="grid"
    cfg.jointPhaseEffectiveBins=cfg.jointPhaseBins;
    cfg.jointPhaseStepDeg=360/cfg.jointPhaseBins;
else
    cfg.jointPhaseEffectiveBins=0; cfg.jointPhaseStepDeg=0;
end
cfg.jointPhaseMaxQuantizationErrorDeg=cfg.jointPhaseStepDeg/2;
if ~isfield(cfg,'b1cBranch'), cfg.b1cBranch='pilot'; end
branch=lower(string(cfg.b1cBranch));
assert(isscalar(branch) && ismember(branch,["data","pilot","both"]),'b1cBranch must be data, pilot or both.');
cfg.b1cBranch=char(branch);
if ~isfield(cfg,'b1cDataPilotWeights'), cfg.b1cDataPilotWeights=sqrt([11 29]); end
v=cfg.b1cDataPilotWeights;
assert(isnumeric(v) && isreal(v) && isvector(v) && numel(v)==2 && ...
    all(isfinite(v)) && all(v>=0) && any(v>0), ...
    'b1cDataPilotWeights must be finite nonnegative [data pilot11] amplitude weights.');
cfg.b1cDataPilotWeights=double(v(:).');
if branch=="both"
    scaled=cfg.b1cDataPilotWeights/max(cfg.b1cDataPilotWeights);
    cfg.b1cEffectiveDataPilotWeights=scaled/norm(scaled);
elseif branch=="data", cfg.b1cEffectiveDataPilotWeights=[1 0];
else, cfg.b1cEffectiveDataPilotWeights=[0 1]; end
if ~isfield(cfg,'b1cPilotPhaseDeg'), cfg.b1cPilotPhaseDeg=90; end
assert(isnumeric(cfg.b1cPilotPhaseDeg) && isreal(cfg.b1cPilotPhaseDeg) && ...
    isscalar(cfg.b1cPilotPhaseDeg) && isfinite(cfg.b1cPilotPhaseDeg), ...
    'b1cPilotPhaseDeg must be a finite real SIGNAL phase in degrees.');
if ~isfield(cfg,'jointWeights'), cfg.jointWeights=[1 1]; end
w=cfg.jointWeights;
assert(isnumeric(w) && isreal(w) && isvector(w) && numel(w)==2 && ...
    all(isfinite(w)) && all(w>=0) && any(w>0), ...
    'jointWeights must be [B1I B1C], finite nonnegative amplitude weights; at least one must be positive.');
cfg.jointWeights=double(w(:).');
if ~isfield(cfg,'jointWeightNormalization')
    if legacy, cfg.jointWeightNormalization='none';
    else, cfg.jointWeightNormalization='l2'; end
end
mode=lower(string(cfg.jointWeightNormalization));
assert(isscalar(mode) && ismember(mode,["l2","none"]),'jointWeightNormalization must be l2 or none.');
cfg.jointWeightNormalization=char(mode);
if mode=="l2"
    scaled=cfg.jointWeights/max(cfg.jointWeights); % stable even for large inputs
    cfg.jointEffectiveWeights=scaled/norm(scaled);
else, cfg.jointEffectiveWeights=cfg.jointWeights; end
% Compatibility: old minBranchRatio sets both SINGLE-signal thresholds only.
if isfield(cfg,'minBranchRatio'), oldSingle=cfg.minBranchRatio; else, oldSingle=8; end
if ~isfield(cfg,'b1iMinRatio'), cfg.b1iMinRatio=oldSingle; end
if ~isfield(cfg,'b1cMinRatio'), cfg.b1cMinRatio=oldSingle; end
if ~isfield(cfg,'minJointScore')
    if mode=="l2", cfg.minJointScore=16; else, cfg.minJointScore=32; end
end
for field={'b1iMinRatio','b1cMinRatio','minJointScore'}
    threshold=cfg.(field{1});
    assert(isnumeric(threshold) && isreal(threshold) && isscalar(threshold) && ...
        isfinite(threshold) && threshold>=0,'%s must be a finite nonnegative threshold.',field{1});
end
% Obsolete joint branch guards must never reintroduce a hidden veto when an
% old options struct is reused. The joint decision is score-only.
if isfield(cfg,'jointMinBranchRatio'), cfg=rmfield(cfg,'jointMinBranchRatio'); end
end
