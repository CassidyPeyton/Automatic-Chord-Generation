%Albert Peyton
function [onsets,onsetStrength,test] = onset_energy(x,win,hop, th)
% Energy-based onset detection

% Input
% - x : audio waveform
% - win : window function
% - hop : window hop size (in samples)
% - th : global threshold to determine onsets

% Output
% - onsets : frame indices of the onsets
% - onsetsStrength : normalized onset strength curve, one value per frame, range in [0, 1]

% Signal envelope
e=zeros(length(x),1);
for n=1:hop:length(x)
    for m=1:length(win)
        if((n+m-(round(length(win)/2)))<1 || (n+m-round((length(win)/2))>=length(x)))
        else
            e(n)=e(n)+(x(n+m-round((length(win)/2)))*win(m))^2;
        end
    end
end
test=e;
onsetStrength=zeros(length(e)-hop,1);
for n=1:hop:length(e)-hop
    onsetStrength(n)=abs(e(n+hop)-e(n));
end
onsetStrength=onsetStrength/max(onsetStrength);
onsets=find(onsetStrength > th);
onsets=round(onsets/length(win));
end