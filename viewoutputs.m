% File: viewoutputs.m
% Created by: M. G. Heinz
% For use with: anmodheinz00.c 
%
% Plots 3 filterbank responses from the AN model
%    BM: Basilar membrane
%   IHC: Inner hair cell
%   IFR: AN instantaneous firing rate

clear

%%%%%%%%%%%%%%%%
% Load and organize all data
%%%%%%%%%%%%%%%%
load ifr.dat
load bm.dat
load ihc.dat
load stim.dat

[Nchs,stimpts]=size(ifr);
stimpts=stimpts-1;
Nchs=Nchs-1;

time=ifr(1,2:stimpts+1);
ANcfs=ifr(2:Nchs+1,1);

bm=bm(2:Nchs+1,2:stimpts+1);
ihc=ihc(2:Nchs+1,2:stimpts+1);
ifr=ifr(2:Nchs+1,2:stimpts+1);
stim=stim(2,2:stimpts+1);


%%%%%%%%%%%%%%%%
% Plots
%%%%%%%%%%%%%%%%

lowch=1;
highch=Nchs;
startms=0;
endms=max(time);


%%%%%%%%%%%%%%%%
%%% Plot BM bank response to tone 
%%%%%%%%%%%%%%%%
figure(1); clf
gain=.5/max(max(bm(lowch:highch,:)));
hold off
for i = lowch:highch,
  plot(time,gain*bm(i,:)+(i));
  hold on
end
plot(time,stim/max(stim)/2+lowch-1.5)
axis([startms endms lowch-2.5 highch+2])
xlabel('Time (ms)')
title(sprintf('Basilar Membrane Filter Output'))
hold off
for i=1:Nchs
  text(1.01*max(time),i,sprintf('%5.0f',ANcfs(i)))
end
text(1.01*max(time),i+2,'CF (Hz)')
text(1.01*max(time),-.5,'Stimulus')

%%%%%%%%%%%%%%%%
%%% Plot IHC bank response to tone 
%%%%%%%%%%%%%%%%
figure(2); clf
gain=.5/max(max(ihc(lowch:highch,:)));
hold off
for i = lowch:highch,
  plot(time,gain*ihc(i,:)+(i));
  hold on
end
plot(time,stim/max(stim)/2+lowch-1.5)
axis([startms endms lowch-2.5 highch+2])
xlabel('Time (ms)')
title(sprintf('Inner Hair Cell Responses'))
hold off
for i=1:Nchs
  text(1.01*max(time),i,sprintf('%5.0f',ANcfs(i)))
end
text(1.01*max(time),i+2,'CF (Hz)')
text(1.01*max(time),-.5,'Stimulus')

%%%%%%%%%%%%%%%%
%%% Plot IFR bank response to tone 
%%%%%%%%%%%%%%%%
figure(3); clf
gain=1/max(max(ifr(lowch:highch,:)));
hold off
for i = lowch:highch,
  plot(time,gain*ifr(i,:)+(i));
  hold on
end
plot(time,stim/max(stim)/2+lowch-1.5)
axis([startms endms lowch-2.5 highch+2])
xlabel('Time (ms)')
title(sprintf('Auditory Nerve Firing Rate'))
hold off
for i=1:Nchs
  text(1.01*max(time),i,sprintf('%5.0f',ANcfs(i)))
end
text(1.01*max(time),i+2,'CF (Hz)')
text(1.01*max(time),-.5,'Stimulus')

