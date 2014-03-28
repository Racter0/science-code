%%%%%%%%%%%% Parameter and run file for the edge-method DQE program, implemented in MATLAB
%%%%% syntax:  /gware/linux64bin/matlab -nodesktop -nosplash -r parfile
%%%%% Assign values to the below variables. 
%%%%% Terminate with a semicolon to suppress output.
%%%%% Strings must be enclosed in single quotes.

%% method: 1 for edge method, 2 for FT method (minimally functional but not refined implementation)

method=1;

%% gain factor: the gain conversion factor of the camera

gain_fac=1;

%% inputfile: the input file name, enclosed in single quotes (eg, '/path/to/dir/image.mrc')

inputfile='/groups/grigorieff/home/ruskinr/DQE/emconvs/K2Count-200kv-50e-3eps.i.mrc';

%% use_noise_file: this program can calculate the noise power spectrum from a single image of a straight
%% edge, or it can use a pure noise image, taken with the exact beam conditions, in addition to the straight edge
%% choose 'y' to use a noise image, 'n' to just use the single image

use_noise_file= 'n';

%% noise_file: if you are using a noise file, this is its filename

noise_file= '/path/to/dir/noise.mrc';

%% lsf_filter_level: default is' none'. However if there is spurious high frequency in the LSF,
%% this value can be set to 'low', med', or 'high' and a modified Hamming filter will be applied to the LSF

lsf_filter_level='none';

%% nps_squash: default is 0.01. nps_squash is the fraction of points at the beginning of the noise power spectrum
%% that are replaced by the average of the subsequent nps_squash/3 points. Useful if there are low-frequency trends
%% and correlations caused by bad gain reference.

nps_squash=0.1;

%% no_align: default is 0. If you think the program's edge alignment isn't helping, (eg, by looking at the figures of pre- and post-alignment
%% edges, you can turn alignment off by setting this to 1.

no_align=0;

%% divide_by_2: default set to 'n'. Some software programs divide 16-bit data by two when saving. If your data has been saved in this way,
%% set divide_by_2 to 'y'. A good sign that you need to do this is if your DQE comes out twice as small as you expect.

divide_by_2 = 'n';

%% steep: default set to 'n'. If you are measuring an especially steep edge spread function, for example from the K2 counting mode, you may
%% want to remove some of the constraints normally applied to the spline-fitting of the ESF so that the sharp features
%% can be captured. In that case, set steep to 'y'. Normally, the constraints prevent smoothing artifacts and should be kept.

steep= 'n';

%% now run the program and save the DQE to DQE.dat

DQE=FindDQE(method, gain_fac, inputfile, use_noise_file, lsf_filter_level, nps_squash, no_align, divide_by_2, steep, noise_file);

save K2Count_200.dat DQE /ascii

