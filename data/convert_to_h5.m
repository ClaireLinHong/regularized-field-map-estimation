% convert data from proprietary .mat format to open .h5 format

% this did not work for complex arrays, so i switched to julia instead!

if 0
 in = load('Ankle_sensemap.mat');
 sensemap = in.sensemap;
 sensemap = single(sensemap); % single precision suffices
 file = 'Ankle_sensemap.h5';
 h5create(file, '/sensemap', size(sensemap))
 h5write(file, '/sensemap', sensemap) % failed :(
end
