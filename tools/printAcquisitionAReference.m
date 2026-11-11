function printAcquisitionAReference(run,cfg,results)
% Console-only strong-A tracking/acquisition reference. Never saves results
% or supplies Doppler/code-phase hypotheses to the weak-B acquisition.
root=fileparts(fileparts(mfilename('fullpath')));
oldPath=path; restore=onCleanup(@()path(oldPath));
sourceA=readLabsat3wInfo(canonical(run.cfg.A.file));
sourceB=readLabsat3wInfo(canonical(cfg.ls3wFile));
assert(isequaln(sourceB,run.B.B1I.settings.sourceInfo) && ...
    abs(cfg.fileStartSecond-run.cfg.B.startSecond)<1e-9 && ...
    abs(cfg.durationSeconds-run.cfg.B.durationSeconds)<1e-9, ...
    'A reference run belongs to another B recording/window.');
fprintf('\nPreparing console-only A reference: short acquisition uses about 200 ms of ORIGINAL strong A.\n');
reference=struct([]);
for signal={'B1I','B1C'}
    name=signal{1};
    q=load(run.A.(name).trackingFile,'key','settings','trackResults');
    assert(isequaln(q.key.source,sourceA) && q.key.start==run.cfg.A.startSecond && ...
        q.key.seconds==run.cfg.A.processSeconds, ...
        'A %s tracking cache changed or belongs to another processing interval.',name);
    if strcmp(name,'B1I'), folder='mBDS_B1I'; period=.001;
    else, folder='mBDS-3_B1C'; period=.010; end
    path(oldPath);
    addpath(root,fullfile(root,'tools'),fullfile(root,folder), ...
        fullfile(root,folder,'include'),fullfile(root,folder,'Common'));
    clear acquisitionCPU acquisitionGPU;
    for k=1:numel(results)
        rec=results(k); prn=rec.PRN;
        model=run.alignment.(name); model=model([model.PRN]==prn);
        track=q.trackResults([q.trackResults.PRN]==prn & [q.trackResults.status]=='T');
        assert(isscalar(model) && isscalar(track), ...
            'Missing unique A tracking/alignment reference for %s C%02d.',name,prn);
        settings=q.settings;
        bSecond=rec.BStartSecond;
        assert(abs(bSecond-cfg.startSecond)<1e-9,'B capture start mismatch.');
        aSecond=model.offsetAminusBSeconds+model.scaleAperB*bSecond;
        aFileSecond=run.cfg.A.startSecond+aSecond;
        span=ceil((settings.codeLength-double(track.remCodePhase(:)))./ ...
            (double(track.codeFreq(:))/settings.samplingFreq))/settings.samplingFreq;
        time=double(track.absoluteSample(:))/settings.samplingFreq+span/2 ...
            -settings.filterDelaySeconds-settings.segmentStartSecond;
        doppler=double(track.carrFreq(:))-settings.IF;
        assert(all(diff(time)>0) && all(isfinite(doppler)), ...
            'Incomplete A %s tracking reference.',name);
        atStart=interp1(time,doppler,aSecond,'linear',nan);
        if strcmp(name,'B1I')
            windowStart=rec.B1IWindowStartSecond; windowEnd=rec.B1IWindowEndSecond;
            bDoppler=rec.B1IDopplerHz; jointDoppler=rec.jointB1CDopplerHz*1561.098/1575.42;
            bCode=rec.B1ICodeMs; jointCode=rec.jointCodeMs;
        else
            windowStart=rec.B1CWindowStartSecond; windowEnd=rec.B1CWindowEndSecond;
            bDoppler=rec.B1CDopplerHz; jointDoppler=rec.jointB1CDopplerHz;
            bCode=rec.B1CCodeMs; jointCode=rec.jointB1CCodeMs;
        end
        bounds=model.offsetAminusBSeconds+model.scaleAperB*[windowStart windowEnd];
        t=[bounds(1);time(time>bounds(1) & time<bounds(2));bounds(2)];
        meanDoppler=trapz(t,interp1(time,doppler,t,'linear',nan))/diff(bounds);
        assert(isfinite(atStart) && isfinite(meanDoppler), ...
            'A tracking does not cover %s C%02d at the requested B time.',name,prn);
        % Native acquisition seeks into the unmodified A IQ cache. Its code
        % phase is local to this NEW read position, not A's initial capture.
        settings.skipNumberOfSamples=round(aFileSecond*settings.samplingFreq);
        settings.acqSatelliteList=prn;
        readSecond=settings.skipNumberOfSamples/settings.samplingFreq;
        originalIF=settings.IF;
        fid=fopen(settings.fileName,'rb'); assert(fid>0,'Missing A IQ conversion cache.');
        closeInput=onCleanup(@()fclose(fid));
        if exist('gpuDeviceCount','file') && gpuDeviceCount>0
            evalc('acquired=acquisitionGPU(fid,settings);');
        else
            evalc('acquired=acquisitionCPU(fid,settings);');
        end
        clear closeInput;
        detected=acquired.carrFreq(prn)~=0 && acquired.peakMetric(prn)>settings.acqThreshold;
        captureDoppler=nan; codeSample=nan; rawDelay=nan; bEquivalentDelay=nan;
        if detected
            carrier=acquired.carrFreq(prn);
            % Original receivers encode exact zero carrier with a 1-Hz sentinel.
            if originalIF==0 && carrier==1, carrier=0; end
            captureDoppler=carrier-originalIF;
            codeSample=acquired.codePhase(prn);
            rawDelay=(codeSample-1)/settings.samplingFreq;
            correctedDelay=mod(rawDelay-settings.filterDelaySeconds+readSecond-aFileSecond,period);
            bEquivalentDelay=correctedDelay/model.scaleAperB;
        end
        item=struct('PRN',prn,'signal',name,'BSecond',bSecond,'ASecond',aSecond, ...
            'AFileSecond',aFileSecond,'AReadSecond',readSecond,'trackingHz',atStart, ...
            'trackingMeanHz',meanDoppler,'detected',detected,'captureHz',captureDoppler, ...
            'peakMetric',acquired.peakMetric(prn),'threshold',settings.acqThreshold, ...
            'codeSample1',codeSample,'rawDelayMs',rawDelay*1000, ...
            'bEquivalentDelayMs',bEquivalentDelay*1000,'BAloneHz',bDoppler, ...
            'BJointHz',jointDoppler,'BAloneCodeMs',bCode,'BJointCodeMs',jointCode);
        if isempty(reference), reference=item; else, reference(end+1)=item; end %#ok<AGROW>
    end
end
fprintf('\n========== A / B acquisition reference (console only; no reference files saved) ==========\n');
fprintf('A source: %s\n',sourceA.fileName);
for prn=unique([reference.PRN])
    rows=reference([reference.PRN]==prn);
    fprintf('\nC%02d | B acquisition start %.9f s\n',prn,rows(1).BSecond);
    for r=rows
        fprintf('  %s: B -> A processing %.9f s, A file %.9f s; actual IQ read %.9f s.\n', ...
            r.signal,r.ASecond,r.AFileSecond,r.AReadSecond);
    end
    fprintf('  Signal | A track@start Hz | A track mean Hz | A short acq Hz | B alone Hz | B joint Hz\n');
    for r=rows
        fprintf('  %-6s | %16.3f | %15.3f | %14.3f | %10.3f | %10.3f\n', ...
            r.signal,r.trackingHz,r.trackingMeanHz,r.captureHz,r.BAloneHz,r.BJointHz);
    end
    fprintf('  Signal | A code sample(1) | A raw lag ms | A->B lag ms | B alone lag ms | B joint lag ms\n');
    for r=rows
        fprintf('  %-6s | %16.0f | %12.9f | %11.9f | %14.9f | %14.9f\n', ...
            r.signal,r.codeSample1,r.rawDelayMs,r.bEquivalentDelayMs,r.BAloneCodeMs,r.BJointCodeMs);
        fprintf('    A %s capture: detected=%d, peak ratio %.4f (threshold %.4f).\n', ...
            r.signal,r.detected,r.peakMetric,r.threshold);
    end
end
fprintf(['A->B lag removes A FIR delay and sample-rounding offset, then applies the A/B clock scale.\n' ...
    'A track mean uses the corresponding B standalone integration window. Frequencies are carrier minus IF.\n' ...
    'Reference values are diagnostics; A/B RF offsets or dynamics can cause differences.\n' ...
    '=======================================================================================\n']);
end
function file=canonical(file)
file=char(java.io.File(file).getCanonicalPath());
end
