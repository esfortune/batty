function signal = lpf( sng, Fs, hilord );
% Signal = wrenfilter( signal, Fs, [ctf order] );
% signal is the data
% Fs is the samplerate in Hz
% Optional... 
%     "low" freq cutoff for lowpass
%     "high" freq cutoff of highpass
%     order of filter (odd number, usual 3-7)
% Default is  ctf=500, order=7

signal = sng - mean(sng);

if nargin == 2 
	order = 7; low = 500; 
end

if nargin == 3
	low = hilord(1); order = hilord(2);
end

        Wnlow = low*2/Fs;
        [bb,aa] = butter(order,Wnlow,'low');
        signal = filtfilt(bb,aa,signal);
