function replica = b1cMatchedReplica(dataCode,pilotCode,nav,secondary,cfg)
% SIGNAL replica. FFT correlation conjugates it, hence +j here means -j
% rotation of the pilot correlation. Never multiply nav and secondary.
switch cfg.b1cBranch
    case 'data'
        replica=dataCode.*single(nav);
    case 'pilot'
        replica=pilotCode.*single(secondary);
    case 'both'
        w=cfg.b1cEffectiveDataPilotWeights;
        phase=single(exp(1j*pi/180*cfg.b1cPilotPhaseDeg));
        replica=single(w(1))*dataCode.*single(nav) + ...
            single(w(2))*phase*pilotCode.*single(secondary);
    otherwise
        error('Unknown b1cBranch.');
end
end
