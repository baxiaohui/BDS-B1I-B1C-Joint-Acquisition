function [filters,rf] = labsatAcquisitionFilters(fs)
% One definition shared by the acquisition reader and noise calibration.
rf=[1561.098e6 1575.42e6]; cut=[2.3e6 5.8e6];
filters=cell(1,2);
for k=1:2
    filters{k}=fir1(320,cut(k)/(fs/2),kaiser(321,7));
end
end
