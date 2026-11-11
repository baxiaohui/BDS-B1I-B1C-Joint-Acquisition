function status = check_public_setup(cfg,checkData)
%CHECK_PUBLIC_SETUP Read-only dependency/data check. Does not acquire/track.
% check_public_setup([],false) checks code/toolboxes without a recording.
root=fileparts(mfilename('fullpath'));
if nargin<1 || isempty(cfg), cfg=ab_stripping_settings(); end
if nargin<2, checkData=true; end
oldPath=path; restorePath=onCleanup(@()path(oldPath));
addpath(fullfile(root,'tools'));
issues={};
required={'fir1','kaiser','pwelch','gf','bchenc'};
for k=1:numel(required)
    if exist(required{k},'file')==0
        issues{end+1}=['Missing MATLAB function: ' required{k}]; %#ok<AGROW>
    end
end
if exist('dsp.FIRDecimator','class')==0
    issues{end+1}='DSP System Toolbox is required (dsp.FIRDecimator).';
end
code={'run_ab_stripping.m','convert_labsat3w.m','tools/alignAidsToB.m', ...
    'tools/aidedJointAcquisition.m','tools/resolveB1CCarrierPolarity.m', ...
    'mBDS_B1I/trkChannelsParallel.m','mBDS_B1I/include/NAVdecoding.m', ...
    'mBDS-3_B1C/trkChannelsParallel.m','mBDS-3_B1C/include/BCNAV1decoding.m', ...
    'mBDS-3_B1C/include/ldpcDecoder/getH_informationSF2.m', ...
    'LS3W_Real_B1I_B1C/Joint_B1I_B1C_PhaseComp_Acquisition.m'};
for k=1:numel(code)
    if ~isfile(fullfile(root,code{k})), issues{end+1}=['Missing source: ' code{k}]; end %#ok<AGROW>
end
if checkData
    for label={'A','B'}
        name=label{1}; file=cfg.(name).file;
        try
            info=readLabsat3wInfo(file);
            assert(info.SMP==58e6,'Requires SMP=58000000.');
            fprintf('%s: %.3f seconds, valid supported LS3W/INI layout.\n',name,info.durationSeconds);
        catch err
            issues{end+1}=sprintf('%s: %s',name,err.message); %#ok<AGROW>
        end
    end
end
status=struct('ready',isempty(issues),'issues',{issues},'dataChecked',logical(checkData), ...
    'matlabVersion',version,'release',version('-release'));
if status.ready
    fprintf('Setup passed (%s). Source and required functions available.\n',status.release);
else
    fprintf('Setup incomplete:\n');
    for k=1:numel(issues), fprintf('  - %s\n',issues{k}); end
end
end
