function beats = beat_dp(onsetStrength, tempoExpected, lambda)
% beat tracking by dynamic programming.
%
% Input
% - wavData : the audio signal
% - onsetStrength : onset strength in each audio frame
% - tempoExpected : the expected tempo (in BPM)
% - lambda : tradeoff between the onset strength objective and beat regularity objective
% Output
% - beats : the estimated beat sequence (in frame number)
fs = 44100;
winsize = 1000;
nframes = length(onsetStrength);
space = round(1/(tempoExpected/60)*fs/winsize);
D = zeros(1,nframes);
P = zeros(1,nframes);
M = zeros(1,nframes);
for n = 1:nframes
            if n == 1
                D(n) =  onsetStrength(n);
            else  
                for m = 1:n-1
                    M(m) = D(m)+lambda*(-(log2((n-m)/space))^2);          
                end 
                D(n) = onsetStrength(n) + max(max(M),0); 
                   if D(n) == onsetStrength(n)
                       P(n) = 0;
                   else
                       [~, P(n)] = max(M);
                   end
            end   
end
l = 1; 
[~, A(l)] = max(D);
while P(A(l)) ~= 0
    l = l + 1;
    A(l) = P(A(l-1));
end 
beats = A;
end