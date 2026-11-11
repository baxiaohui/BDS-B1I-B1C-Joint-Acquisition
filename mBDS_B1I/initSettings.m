function settings = initSettings()
%Functions initializes and saves settings. Settings can be edited inside of
%the function or updated from the command line.
%Edit this file to configure receiver parameters.  
%
%All settings are described inside function code.
%
%settings = initSettings()
%
%   Inputs: none
%
%   Outputs:
%       settings     - Receiver settings (a structure). 

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR  
% (C) Updated by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
% Based on the original work by Darius Plausinaitis,Peter Rinder, 
% Nicolaj Bertelsen and Dennis M. Akos
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

% CVS record:
% $Id: initSettings.m,v 1.9.2.31 2006/08/18 11:41:57 dpl Exp $

%% Processing settings ============================================
% Number of milliseconds to be processed used 36000 + any transients (see
% below - in Nav parameters) to ensure nav subframes are provided
settings.msToProcess        = 45000;        %[ms]
% Number of channels to be used for signal processing
settings.numberOfChannels   = 20;
% Move the starting point of processing. Can be used to start signal
% processing at any point in a long data record. The offset is expressed
% in CONVERTED IQ samples (5.8 Msps for the LS3W B1I profile): one I/Q
% pair counts as one sample. It is not a byte offset in the LS3W file.
settings.skipNumberOfSamples     = 0;

%% Raw signal file name and other parameter =======================
% Original LS3W recording. init.m / postProcessing.m automatically call
% prepareReceiverInput before fread; never read packed LS3W as int8/int16.
% Select 'AB' to run ab_stripping_settings.m using both original receivers.
settings.processingMode = 'single'; % 'single' or 'AB'
settings.signal             = 'B1I';
settings.inputFormat        = 'LS3W';
settings.fileName           = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'data', 'File_048.LS3W');
settings.labsatCacheDir      = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'cache');
% The following describe the CONVERTED IQ BIN, not the packed LS3W file.
settings.dataType           = 'int8';
% File Types
%1 - real samples: S0,S1,S2,...
%2 - complex samples: I0,Q0,I1,Q1,I2,Q2,...                      
settings.fileType           = 2;
% Intermediate, sampling and code frequencies
settings.IF                 = 0;             %[Hz], after digital downconversion
settings.samplingFreq       = 5.8e6;          %[Hz], 58 MHz / 10 after FIR

%% Code parameter setting =================================================
% Define number of chips in a code period
settings.codeLength         = 2046; % Beidou B1I/B2I has 2046 chip length.
settings.codeFreqBasis      = 2.046e6;           %[Hz] 

%% Acquisition settings ===================================================
% Enable use of GPU acceleration for acquisition
settings.gpuACQflag         = 0;             % 0 - Off; 1 - On
% Skips acquisition in the script postProcessing.m if set to 1
settings.skipAcquisition    = 0;
% List of satellites to look for. Some satellites can be excluded to speed
% up acquisition
settings.acqSatelliteList   = 1:63;  % [PRN numbers]   
% One-sided band around IF to search for satellite signal. Depends on the max Doppler.
settings.acqSearchBand      = 5000;              % [Hz]
% Non-coherent integration times after 1ms coherent integration
settings.acqNonCohTime      = 20;                %[ms]
% Threshold for the signal presence decision rule
settings.acqThreshold       = 3;
% Frequency search step for coarse acquisition
settings.acqSearchStep      = 500;               % [Hz]
%% Tracking loops settings ================================================
% Enable/disable use of SIMD/GPU MEX to accelerate tracking
settings.trkMode                 = 1;   % 0 - Channel-serial tracking; 
                                        % 1 - Channel-parallel tracking;
% Enable/disable use of SIMD/GPU MEX to accelerate tracking
settings.correlatorType          = 0;   % 0 - Matlab correlator; 
                                        % 1 - SIMD correlator;
                                        % 2 - GPU correlator;
% Number of right-shift bits for IF data to prevent SIMD correlator 
% overflow when ADC valid bits occupy the high bits
if settings.correlatorType == 1
    settings.rShiftBits              = 5;   % 0 - No shift;
end

% Code tracking loop parameters
settings.dllDampingRatio         = 0.7;
settings.dllNoiseBandwidth       = 1.5;       %[Hz]
settings.dllCorrelatorSpacing    = 0.5;       %[chips]
% Carrier tracking loop parameters
settings.pllDampingRatio         = 0.7;
settings.pllNoiseBandwidth       = 25;        %[Hz]
% Integration time for DLL and PLL
settings.intTime                 = 0.001;      %[s]
%% Navigation solution settings ===================================
% Period for calculating pseudoranges and position
settings.navSolPeriod       = 500;          % [ms]
% Elevation mask to exclude signals from satellites at low elevation
settings.elevationMask      = 5;            % [degrees 0 - 90]
% Enable/disable use of tropospheric correction
settings.useTropCorr        = 0;            % 0 - Off; 1 - On
% True position of the antenna in UTM system (if known). Otherwise enter
% all NaN's and mean position will be used as a reference .
settings.truePosition.E     = nan;
settings.truePosition.N     = nan;
settings.truePosition.U     = nan;

%% Plot settings ==================================================
% Enable/disable plotting of the tracking results for each channel
settings.plotTracking       = 1;            % 0 - Off; 1 - On

%% Constants ======================================================
% The speed of light 
settings.c                  = 299792458;    %[m/s]
% Initial sign. travel time. It is 68.802 ms for BDS MEO and 120 ms for
% GEO/IGSO which is 35,786 km from MSL. Here we take 94 ms as an average value. 
settings.startOffset        = 94;           %[ms] 

%% CNo Settings====================================================
% Accumulation interval in Tracking (in Sec)
settings.CNo.accTime = 0.001;
% Accumulation interval for computing VSM C/No (in ms)
settings.CNo.VSMinterval = 100;

%% BDS B1I/B2I carrier frequency ==================================
% Use 1207.14e6 Hz when this receiver is configured for B2I.
settings.carrFreqBasis    = 1561.098e6;      %[Hz]

%% Navigation TXT export
settings.saveNavMessages = true;
% Continuous on-air symbols, including partial frames; FEC is optional repair.
settings.exportAirBits = true;
settings.exportReceiveTime = true;
settings.messageOutputDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'output');

