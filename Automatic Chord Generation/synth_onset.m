function y = synth_onset(x,frameLen, frameHop, onsets)
% Synthesize each onset as a 1-frame long white noise signal, and add it to
% the original audio signal.
%
% Input
% - x : input audio waveform
% - frameLen : frame length (in samples)
% - frameHop : frame hop size (in samples)
% - onsets : detected onsets (in frames)
% Output
% - y : output audio waveform which is the mixture of x and
% synthesized onset impulses.
noise = zeros(1,frameLen);
for i = 1:frameLen
    noise(i) = rand(1)/2;
end
nframes = ceil(length(x)/frameHop);
x = [x', zeros(1,nframes*frameHop-length(x))];
empty = zeros(1,frameHop*nframes);
for i = 1:length(onsets)
    if onsets(i) > 0
       empty((onsets(i)-1)*frameLen+1 : (onsets(i))*frameLen) = noise;
    end
end
y = x + empty;
end