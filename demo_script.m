% (Change or add directory: .../clust3r)

% Load demo data 
load demo_data

% Generate sorting params structure (optional)
opts.sorterID = 'AN';
opts.experimentTime = '2015-10-21_block-1';
opts.savepath = pwd; % Where to save data files 
opts.figpath = opts.savepath; % Where to save figures
opts.channel = 0;
opts.nWavSamples = size(waves,2); % number of samples per channel (for multichannel waveform data)
opts.tspkUnits = 'ms'; % 's' for seconds, 'ms' for milliseconds, 'samples' for samples
opts.fs = ''; % sample rate (for documentation purpose only) 
    
% Launch clust3r  
[ cluster_indices, waves_out, spiketimes_out, sortinfo] = clust3r( waves, spiketimes, opts );
