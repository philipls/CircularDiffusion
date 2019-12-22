function [theta, thetacc, thetart] = thetaclass4(Data, nw, cond);
% =========================================================================================
% Classify RT and accuracy data by stimulus angle and condition. Medians
% nw is the number of theta bins. 
%       [theta, thetacc, thetart] = thetaclass4(Data, nw, cond);
% Mar 8, 2018: Returns medians.
% =========================================================================================
eps = 0.0001;
sz = size(Data);
if sz(2) ~= 8
   disp('Wrong size data matrix, returning...');
   return
end
condx = 1;
thetax = 2;
fixation_anglex = 4;
fixation_errorx = 6;
rtx = 8;
szrt = 500; % (Arbitrary) buffer size
mq = 0.5;
thetabound = linspace(-pi, pi, nw+1);  % Bounds of theta bins - not really needed
theta = (thetabound(1:nw) + thetabound(2:nw+1))/2; % Centres of theta bins
szt = size(theta);
% Summary structures
thetacc = zeros(size(theta));
thetart = zeros(size(theta));
nacc = zeros(size(theta));
nrt = zeros(size(theta));  % Count the RTs
rt = zeros(szrt, length(theta));

% Pull out one condition.
DataCond = Data(Data(:,condx) == cond,:);
szc = size(DataCond);
ld = szc(1);

for i = 1 : ld
     thetai = DataCond(i, thetax);  % Stimulus angle for trial i ###
     %thetai = Data(i, fixation_anglex); % Classify by response.
     j = bin(thetai, nw) + nw / 2;  % Bin index for the trial i stimulus
     %[thetai, i,j]  
     thetacc(j) = thetacc(j) + DataCond(i,fixation_errorx);
     nacc(j) = nacc(j) + 1; 
     rt(nrt(j) + 1, j) = DataCond(i, rtx);
     thetart(j) = thetart(j) + DataCond(i, rtx);
     nrt(j) = nrt(j) + 1;
end

% Accuracy
for j = 1:length(theta)
    if nacc(j) > 0 
        thetacc(j) = thetacc(j) / nacc(j);
    else
        thetacc(j) = 0;
    end
end

% RT
for j = 1:length(theta)
     rtj = sort(rt(1:nrt(j), j));
     if ~mod(nrt(j), 2) % Even
         thetart(j) = (rtj(nrt(j) / 2) + rtj(nrt(j) / 2 + 1)) / 2;
     else % Odd
         thetart(j) = rtj(floor(nrt(j) / 2) + 1);
     end
end


% end main


function j = bin(thetai, nw);
% ========================================================================================
% Assign the stimulus angle to a bin. (Don't need the bin boundaries to do this.)
% ========================================================================================
   j = ceil((thetai+eps) /(2*pi/nw));




