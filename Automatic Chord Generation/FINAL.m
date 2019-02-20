%Albert Peyton Ryan Bhular Sheng Xu 2018
[wavdata, fs]= audioread('guitar1.wav');
%first 30 seconds of recording
%piano 103bpm
%guitar 106bpm
wavdata=wavdata(1:length(wavdata),1);
%call the yin function to find pitches

tHop=.0464;
tW=.0464;
[pitch,~,~]= myYin(wavdata, fs, tHop, 0.0464, 40, 2000, 0.1);
for i=1:length(pitch)
    if pitch(i)==0
    else
        pitch(i)=69+12*log2(pitch(i)/440);
        pitch(i)=round(pitch(i));
        pitch(i)=mod(pitch(i),12);
    end
end
%creating the frames x 12 array for matrix multiplication

ditch=zeros(length(pitch),12);
for i=1:length(pitch)
    if pitch(i)==0
        ditch(i,12)=1;
    else
        ditch(i,pitch(i))=1;
    end
end
prob=zeros(12,12);
%probability of scale degree in each key
%sample=rand(12,1);
%sample([2,4,7,9,11])=0;
sample=[1,0,1,0,1,1,0,1,0,1,0,1];
%creating a matrix of the
for i=1:12
    if i==1
        prob(:,i)=sample;
    else
        prob(:,i)=circshift(sample,i-1);
    end
end
test=ditch*prob;
summer=sum(test);
[M, i1]=max(summer);
summer(i1)=0;
[M, i2]=max(summer);
%% key identification
result=prob(:,i1)+prob(:,i2);
for i=1:12
   if result(i)<.5
      result(i)=0; 
   else
       result(i)=1;
   end
end
%% tempoestimation
%[bpm,~]=beat_tracking(wavdata,fs);
%% beat detection
bpm=90;
hop=round(tHop*fs);
win=hamming(tW*fs);
[onsets, onsetStrength] = onset_spectral(wavdata, win, hop, .3,100,fs );
beats=beat(onsetStrength, 106, 4);
beats=sort(beats);
%% beat detection text
z=synth_onset(wavdata,hop, hop, beats);
%%
% our classical transition matrix
tran=[.1 .05 .1 .25 .25 .2 .05; .05 .1 .1 .2 .25 .2 .15; .15 .05 .1 .15 .2 .25 .1; .2 .1 .05 .1 .2 .2 .15; .4 .2 .1 .05 .1 .05 .1; .1 .25 .25 .2 .05 .1 .05; .25 .1 .1 .1 .3 .05 .1];
%observation sequence of notes
O=pitch(beats(1:length(beats))+2);
obs=zeros(1,length(beats));
%converting pitches to scale degrees
for i=1:length(beats)
    deg=mod(O(i)+(12-i1),12);
    if deg==0
        obs(i)=1;
    elseif  deg==2
        obs(i)=2;
    elseif  deg==4
        obs(i)=3;
    elseif  deg==5
        obs(i)=4;
    elseif  deg==7
        obs(i)=5;
    elseif  deg==9
        obs(i)=6;
    elseif  deg==11
        obs(i)=7;
    end
end
%%
%major/minor detection
min=0;
maj=0;
for i=1:length(obs)
   if obs(i)==6||obs(i)==7||obs(i)==3
       min=min+1;
       
   elseif obs(i)==1||obs(i)==5||obs(i)==2
       maj=maj+1;
   end
       
end
minor=0;
if min>maj
    minor=1;
    obs=obs+2;
    obs=mod(obs,7);
    for i=1:length(obs)
       if obs(i)==0
           obs(i)=7;
       end
    end
end
%%

%emission matrix calculation
emis=zeros(7,length(obs));
for i=1:length(obs)
    if obs(i)==1
        emis(1,i)=1/2;
        emis(4,i)=1/3;
        emis(6,i)=1/6;
    elseif obs(i)==2
        emis(2,i)=1/2;
        emis(5,i)=1/3;
        emis(7,i)=1/6;        
    elseif obs(i)==3
        emis(3,i)=1/2;
        emis(6,i)=1/3;
        emis(1,i)=1/6;       
    elseif obs(i)==4
        emis(4,i)=1/6;
        emis(7,i)=1/6;
        emis(2,i)=1/6;
        emis(5,i)=1/2;
    elseif obs(i)==5
        emis(5,i)=1/2;
        emis(1,i)=1/3;
        emis(3,i)=1/6;
    elseif obs(i)==6
        emis(6,i)=1/3;
        emis(2,i)=1/2;
        emis(4,i)=1/6;
    elseif obs(i)==7
        emis(7,i)=1/6;
        emis(3,i)=1/4;
        emis(5,i)=1/4;        
    end
end
emis=emis+.00001;
%inital probabilites based on first note
if obs(1)==1||obs(1)==3
    pi=[.6,0,0,0,.4,0,0];
elseif obs(1)==5
    pi=[.5,0,0,0,.5,0,0];
elseif obs(1)==2
    pi=[0,.2,0,0,.7,0,0];
elseif obs(1)==7
    pi=[.3,0,0,0,.7,0,0];
elseif obs(1)==4||obs(1)==6
    pi=[0,.2,0,.8,0,0,0];
end
%[a,b,p]=BaumWelch_n(transMat',emis,pi,obs',200)
%% running viterbi with the different transition matricies
%our pop matrix
load('matlab.mat');
[~,~,pathclass]=ShengViterbi(tran,emis,pi);
[~,~,pathpop]=ShengViterbi(transMat',emis,pi);
%[~,~,pathbw]=ShengViterbi(a,b,p);
%%
chord=zeros(length(beats),3);
%major
for i=1:length(beats)
    if pathclass(i)==1
        chord(i,:)=[0 4 7] + i1;
    elseif pathclass(i)==2
        chord(i,:)=[2 5 9] + i1;
    elseif pathclass(i)==3
        chord(i,:)=[4 7 11] + i1;
    elseif pathclass(i)==4
        chord(i,:)=[5 9 12] + i1;
    elseif pathclass(i)==5
        chord(i,:)=[7 11 14] + i1;
    elseif pathclass(i)==6
        chord(i,:)=[9 12 16] + i1;
    elseif pathclass(i)==7
        chord(i,:)=[11 14 17] + i1;
    end
end
% %minor
if minor==1
for i=1:length(beats)
    if pathclass(i)==1
        chord(i,:)=[0 3 7] + i1-3;
    elseif pathclass(i)==2
        chord(i,:)=[2 5 8] + i1-3;
    elseif pathclass(i)==3
        chord(i,:)=[3 7 10] + i1-3;
    elseif pathclass(i)==4
        chord(i,:)=[5 8 12] + i1-3;
    elseif pathclass(i)==5
        chord(i,:)=[7 10 14] + i1-3;
    elseif pathclass(i)==6
        chord(i,:)=[8 12 15] + i1-3;
    elseif pathclass(i)==7
        chord(i,:)=[10 14 17] + i1-3;
    end
end
end
%% creating midi files
chord=chord+48;

% block chords
M = zeros(3*length(beats),6);
M(:,1)=1;
M(:,2)=1;
M(:,4)=120;
for i = 1:length(beats)
    M((1+3*(i-1)):(3+3*(i-1)),3)=chord(i,:);
    M((1+3*(i-1)):(3+3*(i-1)),5)=beats(i)*tW;
    M((1+3*(i-1)):(3+3*(i-1)),6)=M((1+3*(i-1)):(3+3*(i-1)),5)+.5;
end

midi=matrix2midi(M);
writemidi(midi,'midichords.mid');
% arpeggios
N = zeros(4*length(beats),6);
N(:,1)=1;
N(:,2)=1;
N(:,4)=120;
chord(:,4)=chord(:,2);
for i = 1:length(beats)
    N((1+4*(i-1)):(4+4*(i-1)),3)=chord(i,:);

end
%subdivisions of beat
for i=1:length(beats)-1
d(i)=(beats(i+1)-beats(i));
end
q=median(d)/4;
for i=1:4*length(beats)
    N(i,5)=beats(1)*tW+q*(i-1)*tW;
    N(i,6)=N(i,5)+.25;
end
midi2=matrix2midi(N);
writemidi(midi2,'midiarpeggio.mid');

%% importing midi for demo
[y,Fs]=midi2audio(midi,fs, 'sine');
z1=.25*y+wavdata(1:length(y),1)';
z1(length(y):length(wavdata))=wavdata(length(y):length(wavdata))';
%% importing midi for demo (arrpegios)
[y1,Fs]=midi2audio(midi2,fs, 'sine');
if length(y1)>length(wavdata)
z2=.25*y1(1,1:length(wavdata))+wavdata';
z2(length(y1):length(wavdata))=wavdata(length(y1):length(wavdata))';
else
z2=.25*y1+wavdata(1:length(y1),1)';
z2(length(y1):length(wavdata))=wavdata(length(y1):length(wavdata))';
end
%% input signal
sound(wavdata,fs)
%% added block chords
sound(z1,Fs);
%% added arpeggios
sound(z2,Fs);