function [onsets, onsetStrength] = onset_spectral(x, win, hop, th, gamma,fs)
% Spectral-based onset detection (spectral flux). The spectral flux calculation is based on the compressed magnitude spectrogram: spectrogram_compressed = log(1 + gamma * spectrogram). This is to reduce the dynamic range of the spectrogram.
%
% Input
% - x : audio waveform
% - win : window function
% - hop : window hop size (in samples)
% - th : threshold to determine onsets
% - gamma : parameter for spectrogram compression
% Output
% - onsets : frame indices of the onsets
% - onsetsStrength : normalized onset strength curve, one value per frame, range in [0, 1]
f = length(win);
s=spectrogram(x,win,(length(win)-hop),f,fs);
Ya = log(1+gamma*abs(s));
[a b]=size(Ya);
nframes =b,;
Spec_delta = zeros(1,nframes);
 for n = 1:nframes-1
     for k = 1:(round(length(win)/2))
        Spec_delta(n)= Spec_delta(n)+ abs(Ya(k,n+1)-Ya(k,n));
     end
 end
 % normalized onset strength
 onsetStrength = (Spec_delta-min(Spec_delta))/(max(Spec_delta) - min(Spec_delta));
 %local average
 u = zeros(1,nframes);
 M = 1000;
 onsetStrength_S = [onsetStrength, zeros(1,M)];
 for i = 1:nframes
     for m = 1:M 
        u(n) = u(n) + 1/M *  onsetStrength_S(n+m);
     end
 end
 bar_onsetStrength = zeros(1,nframes);
for i = 1:nframes
    bar_onsetStrength(i) = abs(onsetStrength(i) - u(i));
end
 onsets = find(bar_onsetStrength > th);
end