function [theta,residual,index] = jointPhaseHypothesis(difference,mode,bins)
% difference = arg(Z_C)-arg(Z_I); compensate I by exp(+1j*theta).
% For nonnegative amplitude weights, the nearest uniform phase maximizes
% the cosine cross term EXACTLY; no K-times loop or extra phase dimension.
% Continuous theta retains the old (possibly unwrapped) angle convention.
switch char(mode)
    case 'continuous'
        theta=double(difference);
        residual=zeros(size(theta),'like',theta);
        index=nan(size(theta),'like',theta);
    case 'grid'
        assert(isscalar(bins) && isfinite(bins) && bins>=1 && bins==fix(bins), ...
            'jointPhaseBins must be a positive integer.');
        step=2*pi/double(bins);
        delta=double(difference);
        index=mod(round(delta/step),double(bins)); % 0-based in [0,K-1]
        theta=index*step;
        residual=mod(delta-theta+pi,2*pi)-pi;
    otherwise
        error('jointPhaseMode must be continuous or grid.');
end
end
