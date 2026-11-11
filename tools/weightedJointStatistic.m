function [score,effective] = weightedJointStatistic(powerI,powerC,weights,normalization,phaseResidual)
% Phase-optimized coherent statistic; inputs are background-normalized powers.
% Weights apply to COMPLEX AMPLITUDE, not directly to power.
weights=double(weights(:).');
assert(isreal(weights) && numel(weights)==2 && all(isfinite(weights)) && ...
    all(weights>=0) && any(weights>0),'Invalid joint amplitude weights.');
assert(ismember(lower(string(normalization)),["l2","none"]),'Unknown weight normalization.');
if strcmpi(normalization,'l2')
    scaled=weights/max(weights); effective=scaled/norm(scaled);
else, effective=weights; end
a=effective(1)*sqrt(powerI); b=effective(2)*sqrt(powerC);
score=(a+b).^2;
if nargin>=5 && all(effective>0)
    % Equivalent to a^2+b^2+2*a*b*cos(residual), but preserves the old
    % continuous result and avoids cancellation in 1-cos(small residual).
    score=max(score-4*a.*b.*sin(phaseResidual/2).^2,0);
end
end
