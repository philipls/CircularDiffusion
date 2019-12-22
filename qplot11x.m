function qplot11x(fitfunc, Pvar, Pfix, Sel, Data_in);
% ----------------------------------------------------------------------
% Empirical and RT fitted RT quantiles as a function of theta.
% Plots using <aicircle11> predictions
% Equal-mass theta bounds, RTs filtered on [minrt, maxrt]
%       qplot11x(@fitfunc, Pvar, Pfix, Sel, Data_in);
% ----------------------------------------------------------------------
    minrt = 0.15;  % Filter
    maxrt = 4.5;

    errmg2 = 'QPLOT: Wrong size data matrix, exiting...';

    if size(Data_in) ~= [1,3]
       disp('Wrong size data matrix, exiting...')
       return
    else
       sz = size(Data_in{1})
       two_column = sz(2) == 2;  
    end
    if two_column
        Data = Data_in;
    else
        disp('Three column')
        Data =  cell(1,3);
        Data{1} = Data_in{1}(:, 2:3);
        Data{2} = Data_in{2}(:, 2:3);
        Data{3} = Data_in{3}(:, 2:3);
    end; 

    co = [0,0,1; 
         0,0.5,0;
         1,0,0;
         0,0.75,0.75;
         0.75,0,0.75;
         0.75,0.75,0;
         0.25,0.25,0.25];

    set(groot,'defaultAxesColorOrder',co);
    epsx = 1e-9;
    tmax = 4.0;
    labs = {'High'; 'Med'; 'Low'};

    nmass = 10;
    massm = 1.0/nmass; % 10 bins

    % Generic fit function call.
    [ll,aic, bic,Pred,Gstuff] = fitfunc(Pvar, Pfix, Sel, Data_in)

    ta = Gstuff{1,1};
    thetaa = Gstuff{2,1};
    gtma = Gstuff{3,1};
    tb = Gstuff{1,2};
    thetab = Gstuff{2,2};
    gtmb = Gstuff{3,2};
    tc = Gstuff{1,3};
    thetac = Gstuff{2,3};
    gtmc = Gstuff{3,3};

    axhandle =setfig3;
    qploti(co, axhandle(1), Data{1}, gtma, ta, thetaa, tmax, minrt, maxrt, 'High', 0);
    qploti(co, axhandle(2), Data{2}, gtmb, tb, thetab, tmax, minrt, maxrt, 'Med', 0);
    qploti(co, axhandle(3), Data{3}, gtmc, tc, thetac, tmax, minrt, maxrt, 'Low', 1);
end

function qploti(co, axi, Datai, Gt, T, thetai, tmax, minrt, maxrt, labi, do_xlabel)
% =======================================================================
% Plot 7 empirical distribution quantiles against predictions for 3
% discriminability conditions.
% =======================================================================
  %disp('in qploti')
  %size(Datai)

   symbol = ['o', 's', 'd', 'v', '^'];

   h = tmax / 300; 
   nw = 50; 
   w = 2*pi/nw;

   lnt = length(T);
   Ft = cumsum(Gt, 2) * h  * (2 * pi / nw); % normalize circular mass
   %size(Ft)
   MaxFt = Ft(:,lnt) * ones(1, lnt);
   % Calculate normalized (conditional) distribution functions.
   NormFt = Ft./ MaxFt;
   % Calculate quantiles
   Qf5 = [.1, .3, .5, .7, .9]; % Summary quantiles (can be changed).
   Qt = zeros(nw, 5);
   for i = 1:nw
       Fti = NormFt(i,:);
       Ix = (Fti >= .025 & Fti <= .975);
       if min(diff(Fti(Ix))) <= 0 
             Qti = [0,0,0,0,0];
             disp('Cannot compute Ft quantiles.');
             i
       else
             Qti=interp1(Fti(Ix)', T(Ix)', Qf5);
       end;
       Qt(i,:) = Qti;
   end;

   axes(axi);

  % Empirical RT quantiles in accuracy bins
   [Q,ThetaCentres] = bin9(Datai, minrt, maxrt);
   bound = 1.25 *  max(abs(ThetaCentres));
   theta = thetai(1:nw)'; % Because of wrap-around.
   Ix = theta >= -bound & theta <= bound;
   plot(theta(Ix), Qt(Ix, 1), '-.k', ...
        theta(Ix), Qt(Ix, 2), '-.k', ...
        theta(Ix), Qt(Ix, 3), '-.k', ...
        theta(Ix), Qt(Ix, 4), '-.k', ...
        theta(Ix), Qt(Ix, 5), '-.k')
   c = get(gca, 'Child');
   set(c(1), 'Linewidth', 2);

   hold
   set(gca, 'XLim', [-2.5,2.5]);

   set(gca, 'YLim', [0.5, 3.5]);  % maxrt
   if do_xlabel
       xlabel('Hitting Angle (rad)')
   end
   ylabel('Distribution Quantile (s)')

 
   %ThetaCentres
   %Q
   for j = 1:5
      plot(ThetaCentres, Q(j, :), 'k-')
   end
   for j = 1:5
       plot(ThetaCentres, Q(j,:), symbol(j), 'MarkerSize', 6, ...
       'MarkerEdgeColor', co(j,:), 'MarkerFaceColor', co(j,:));
   end
   label(gca, .6, .89, labi);
end


function [Q, ThetaCentres] = bin9(Data, minrt, maxrt);
% ========================================================================================
%    [Q, ThetaCentres] = bin9(Data, minrt, maxrt)
%    Bin RTs into 9 equal-mass bins, filter out long RTs.
% ========================================================================================
eps = 0.0001;

% Equal-mass theta boundaries
ntheta = 7;
%ntheta = 9;
lnd = length(Data);
thetabin = round(lnd * [0.1429, 0.2857, 0.4286,  0.5714, 0.7143, 0.8571, 1.0000]);
%thetabin = round(lnd * [0.1111, 0.2222, 0.3333,  0.4444, 0.5556, 0.6667, 0.7778, 0.8889, 1.0000])

BinRT = zeros(110, ntheta);
BinTheta = zeros(1, ntheta);
BinCount = zeros(1, ntheta);

[thetas,I] = sort(Data(:,1));
Data(:,:) = Data(I,:);

j = 1;
k = 1;
for i = 1:lnd
    % Go to next bin (data sorted by ascending theta) and reset RT counter
    if i > thetabin(j)
         j = j + 1;
         k = 1;
    end
    thetai = Data(i,1);
   %j
    BinTheta(j) = BinTheta(j) + thetai;  % Sum thetas
    BinRT(k, j) = Data(i, 2);
    %[i, j, k, Data(i, 2)]
    %pause
    BinCount(j) = BinCount(j) + 1;
    k = k + 1;
end
ThetaCentres = BinTheta ./ BinCount;

% Filter outliers, sort RTs in each bin
Q = zeros(5,ntheta);
Qp =[.1,.3,.5,.7,.9];
for j = 1 : ntheta
    rt = BinRT(1:BinCount(j), j);
    rts = sort(rt);
    truncrt = rts(find(rts >= minrt & rts <= maxrt));
    Qx = ceil(length(truncrt) * Qp);
    Q(:,j) = rts(Qx); 
end
end




