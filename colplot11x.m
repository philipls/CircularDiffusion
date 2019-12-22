function colplot11x(Data, Pred)
% ========================================================================
% Plot fitted values of circular diffusion model with drift variability.
%     colplot11x(Data, Pred)
% Works with either 2-column or 3-column data
% ========================================================================
name = 'COLPLOT11X: ';
errmg1 = 'Data should be a 1 x 3 cell array from <makelike>...';
if size(Data) ~= [1,3]
   disp('Wrong size data matrix, exiting...')
   return
else
  sz = size(Data{1})
  two_column = sz(2) == 2; % Canonical orientation 
end
two_column

axhandle = setfig6x;
if two_column % Canonical orientation
    Theta1 = Data{1}(:,1);
    Rt1 = Data{1}(:,2);
    Theta2 = Data{2}(:,1);
    Rt2 = Data{2}(:,2);
    Theta3 = Data{3}(:,1);
    Rt3 = Data{3}(:,2);
else
    Theta1 = Data{1}(:,2);
    Rt1 = Data{1}(:,3);
    Theta2 = Data{2}(:,2);
    Rt2 = Data{2}(:,3);
    Theta3 = Data{3}(:,2);
    Rt3 = Data{3}(:,3);
end
% Predictions are a 2 x 3 cell array, 1st row is joint density, 2nd row is accuracy 
ta = Pred{1,1}(1,:);
gtam = Pred{1,1}(2,:);

tb = Pred{1,2}(1,:);
gtbm = Pred{1,2}(2,:);

tc = Pred{1,3}(1,:);
gtcm = Pred{1,3}(2,:);

thetaa = Pred{2,1}(1,:);
pthetaa = Pred{2,1}(2,:);

thetab = Pred{2,2}(1,:);
pthetab = Pred{2,2}(2,:);

thetac = Pred{2,3}(1,:);
pthetac = Pred{2,3}(2,:);

%cvec1 = [.5, .65, .85];   % Kinda blue...
%cvec1 = [.65, .85, .95];   % Kinda blue...
%cvec2 = [.85, 0, .85]   % Dark magenta
cvec2 = [.60, 0, .60]   % Dark magenta
cvec1 = [.90, .95, .95];   % Kinda blue..
%Accuracy low
axes(axhandle(1))
histogram(Theta1, 50, 'Normalization', 'pdf', 'BinLimits', [-pi,pi]);
set(gca, 'Xlim', [-pi, pi])
set(gca, 'Ylim', [0, 1.5])
xlabel('Response Error (rad)')
ylabel('Probability density')
label(gca, .65, .85, 'High');
hold
plot(thetaa, pthetaa, 'm-', 'Linewidth', 2);
c = get(gca, 'Child');
c(1).Color = cvec2;
c(3).FaceColor  = cvec1;
c(2).Color = 'k';

% Accuracy med
axes(axhandle(3));
histogram(Theta2, 50, 'Normalization', 'pdf', 'BinLimits', [-pi,pi]);
set(gca, 'Xlim', [-pi, pi])
set(gca, 'Ylim', [0, 1.5])
xlabel('Response Error (rad)')
ylabel('Probability density')
label(gca, .65, .85, 'Med');
hold
plot(thetab, pthetab, 'm-', 'Linewidth', 2);
c = get(gca, 'Child');
c(1).Color = cvec2;
c(3).FaceColor  = cvec1;
c(2).Color = 'k';

% Accuracy hi
axes(axhandle(5));
histogram(Theta3, 50, 'Normalization', 'pdf', 'BinLimits', [-pi,pi]);
set(gca, 'Xlim', [-pi, pi])
set(gca, 'Ylim', [0, 1.5])
xlabel('Response Error (rad)')
ylabel('Probability density')
label(gca, .65, .85, 'Low');
hold
plot(thetac, pthetac, 'm-', 'Linewidth', 2);
c = get(gca, 'Child');
c(1).Color = cvec2;
c(3).FaceColor  = cvec1;
c(2).Color = 'k';

% RT low
axes(axhandle(2));
histogram(Rt1, 50, 'Normalization', 'pdf', 'BinLimits', [0,4.5]);
set(gca, 'Xlim', [0, 4.0])
xlabel('Response Time (s)')
set(gca, 'Ylim', [0, 3.75])
label(gca, .65, .85, 'High');
hold
plot(ta, gtam, 'm-', 'Linewidth', 2);
c = get(gca, 'Child');
c(1).Color = cvec2;
c(3).FaceColor = cvec1;
c(2).Color = 'k';

% RT med
axes(axhandle(4));
histogram(Rt2, 50, 'Normalization', 'pdf', 'BinLimits', [0,4.5]);
set(gca, 'Xlim', [0, 4.0])
xlabel('Response Time (s)')
set(gca, 'Ylim', [0, 3.75])
label(gca, .65, .85, 'Med');
hold
plot(tb, gtbm, 'm-', 'Linewidth', 2);
c = get(gca, 'Child');
c(1).Color = cvec2;
c(3).FaceColor  = cvec1;
c(2).Color = 'k';

% RT hi
axes(axhandle(6));
histogram(Rt3, 50, 'Normalization', 'pdf', 'BinLimits', [0,4.5]);
set(gca, 'Xlim', [0, 4.0])
xlabel('Response Time (s)')
set(gca, 'Ylim', [0, 3.75])
label(gca, .65, .85, 'Low');
hold
plot(tc, gtcm, 'm-', 'Linewidth', 2);
c = get(gca, 'Child');
c(1).Color = cvec2;
c(3).FaceColor  = cvec1;
c(2).Color = 'k';






