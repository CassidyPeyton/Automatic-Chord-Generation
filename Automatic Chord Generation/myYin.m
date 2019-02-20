function [pitch, ap_pwr, rms1] = myYin(wavData, fs, tHop, tW, f0Min, f0Max, dp_th)

% My YIN implementation of Step 1-5.
% The first estimate is performed at time 0 (i.e. the first integration
% window is from the first sample to time tW, which means the window is
% centered at tW/2) and the following estimates are spaced by tHop.
%
% Input
%   - wavData       : input single-channel audio wave (a column vector)
%   - fs            : sampling rate of wavData
%   - tHop          : time interval between two adjacent estimates in second (default 0.01)
%   - tW            : integration window size in second (default 0.025)
%   - f0Min         : lowest possible F0 in Hz (default 40)
%   - f0Max         : highest possible F0 in Hz (default 400)
%   - dp_th         : the threshold of dips of d_prime(default 0.1)
% Output
%   - pitch         : estimated pitches (a row vector)
%   - ap_pwr        : the corresponding d_prime, which is approximatedly
%                       the aperiodic power over the total signal power.
%                       It can be used as the salience of this pitch
%                       estimate.
%   - rms           : the RMS value of the signal in the integration
%                       window, which can also be used to determine if
%                       there is a pitch in the signal or not.
%
% Author: XXX
% Created: XXX
% Last modified: XXX
%load('violin_gt_pitch.mat')
% default parameters for speech
if nargin<7 dp_th=0.1; end
if nargin<6 f0Max=400; end
if nargin<5 f0Min=40; end
if nargin<4 tW=0.025; end
if nargin<3 tHop=0.01; end

% Part 2
%converting seconds into samples
winlen=round(tW*fs);
hopsize=round(tHop*fs);
N=length(wavData);
rmx=rms(wavData);
%number of frames
Nframes=ceil((N-winlen)/hopsize);
%zeropadding
x= [wavData', zeros(1,Nframes*winlen-N)];
frame=zeros(Nframes, winlen);
init=1;
%creating frames
rms1=zeros(1,Nframes);
for i= 1:Nframes
    frame(i,:)=x(init:(init+winlen-1));
    init=init+hopsize;
    i
    rms1(i)=rms(frame(i,:));
   % if -30>20*log10(rms1(i)/rmx)
        
  %  end
end
d=zeros(Nframes,winlen);
x_temp = [frame, zeros(Nframes,winlen)];
%differencing
for t = 0:winlen-1
    for j = 1:winlen
        d(:,t+1)=d(:,t+1)+ (x_temp(:,j)- x_temp(:,j+t)).^2;
    end
end

%Part 3
%normalizing the difference function
ap_pwr = zeros(Nframes,winlen);
ap_pwr(:,1) = 1;
for i = 1:Nframes
    for t = 1:(winlen -1)
        ap_pwr(i,t+1)=d(i,t+1)/((1/t)*sum(d(i,1:t+1)));
    end
    i
end
%Part 4 
%absoulte threshold
lag = zeros(1,Nframes);
dp_th = .1;
for i = 1:Nframes
    %finds value with lowest dip
    dip = find(ap_pwr(i,:) < dp_th,1);
    if isempty(dip) == 1
        [v,dip] = min(ap_pwr(i,:));
    end
    lag(i)=dip;
    i
end

%Part 5
period= zeros(1,Nframes);
time = zeros(Nframes,winlen);
pitch = zeros(Nframes,1);


for i=1:Nframes
    if(lag(i) > 1 && lag(i) < winlen)
        a = ap_pwr(i,lag(i)-1);
        b = ap_pwr(i,lag(i));
        gamma = ap_pwr(i,(lag(i)+1));
        peak = .5*(a-gamma)/(a-2*b+gamma);
    else
        peak=0;
    end
    period(i)=(lag(i)-1)+peak;
    pitch(i)=fs/period(i);
    time(i,:) = ((i-1)*winlen:i*winlen-1)/fs;
end
ap_pwr=mean(ap_pwr');


%zero the threshold
for i = 1:Nframes
    if pitch(i) > f0Max
        %
        pitch(i) = 0;
    end
    if pitch(i) < f0Min
        %
        pitch(i) = 0;
    end
end
end