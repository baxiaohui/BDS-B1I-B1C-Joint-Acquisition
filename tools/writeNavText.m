function files = writeNavText(track, messages, settings, signal, outputDir)
% Export only timestamped navigation bits and NH/pilot secondary code.
if ~isfolder(outputDir), mkdir(outputDir); end
files=writeStrippingText(track,messages,settings,signal,outputDir);
end
