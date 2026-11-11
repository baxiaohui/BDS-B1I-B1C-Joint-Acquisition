function result = run_public(cfg)
%RUN_PUBLIC Strong receivers -> air-symbol alignment -> aided joint acquisition.
% result=run_public; or cfg=ab_stripping_settings(); result=run_public(cfg);
root=fileparts(mfilename('fullpath'));
oldPath=path; restorePath=onCleanup(@()path(oldPath));
addpath(root,fullfile(root,'tools'));
if nargin<1, cfg=ab_stripping_settings(); end
status=check_public_setup(cfg);
assert(status.ready,'PublicRelease:Setup','Resolve setup issues listed above before running.');
result=run_ab_stripping(cfg);
end
