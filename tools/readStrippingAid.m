function aid = readStrippingAid(file)
% Read the simple timestamped TXT format, including optional receive time.
fid=fopen(file,'rt','n','UTF-8'); assert(fid>0,'Aid file missing: %s',file);
closeFile=onCleanup(@()fclose(fid)); header={};
while true
    line=fgetl(fid); assert(ischar(line),'Missing column header in %s.',file);
    if startsWith(line,'BDT_SOW'), names=strsplit(strtrim(line)); break; end
    header{end+1}=line; %#ok<AGROW>
end
assert(numel(names)>=5,'Use the A/B export with RX_SECONDS and PHASE_INDEX: %s',file);
values=textscan(fid,repmat('%f',1,numel(names)),'CollectOutput',true);
aid.rows=values{1}; aid.header=strjoin(header,newline); aid.file=file;
assert(~isempty(aid.rows) && size(aid.rows,2)>=5,'Empty/invalid aid: %s',file);
assert(all(ismember(aid.rows(:,2),[0 1])) && ...
    all(aid.rows(:,3)==1-2*aid.rows(:,2)) && all(diff(aid.rows(:,4))>0), ...
    'Invalid bits or receive timestamps: %s',file);
tok=regexp(aid.header,'BDT_WEEK=(\d+)','tokens','once'); aid.week=nan;
if ~isempty(tok), aid.week=str2double(tok{1}); end
tok=regexp(aid.header,'bit_duration_seconds=([\d.]+)','tokens','once');
if ~isempty(tok), aid.step=str2double(tok{1});
elseif contains(file,'_NH.txt'), aid.step=0.001;
else, aid.step=0.010; end
if size(aid.rows,2)<6, aid.rows(:,6)=aid.rows(:,4)+aid.step; end
end
