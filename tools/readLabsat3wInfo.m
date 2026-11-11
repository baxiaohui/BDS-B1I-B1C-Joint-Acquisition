function info = readLabsat3wInfo(fileName)
% Read and validate the LS3W sidecar. This reader supports QUA=1, CHN=3,
% SFT=6 (the supplied File_048 recording); reject other layouts explicitly.
[folder, stem] = fileparts(fileName);
iniName = fullfile(folder, [stem '.ini']);
assert(isfile(fileName), 'LS3W data file does not exist: %s', fileName);
assert(isfile(iniName), 'LS3W sidecar does not exist: %s', iniName);
ini = fileread(iniName);
keys = {'SMP','QUA','CHN','SFT','CFA','BWA'};
for k = 1:numel(keys)
    token = regexp(ini, ['(?m)^' keys{k} '\s*=\s*(\d+)'], 'tokens', 'once');
    assert(~isempty(token), 'Missing INI key: %s', keys{k});
    info.(keys{k}) = str2double(token{1});
end
assert(info.QUA == 1 && info.CHN == 3 && info.SFT == 6, ...
    'Supported LS3W layout: QUA=1, CHN=3, SFT=6 only.');
info.fileName = char(fileName);
info.iniName = iniName;
info.samplesPerRegister = 10;
info.spareBits = 4;
d = dir(fileName);
info.fileBytes=d.bytes;
info.modifiedDatenum=d.datenum;
assert(mod(d.bytes,8) == 0, 'Incomplete LS3W 64-bit register at EOF.');
info.totalSamples = d.bytes/8*10;
info.durationSeconds = info.totalSamples/info.SMP;
end
