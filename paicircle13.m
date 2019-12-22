function [ll,qaic,qbic, Pred,Gstuff] = paicircle13(Pvar, Pfix, Sel, Data, trace)
% ==========================================================================
% Circular diffusion with drift anisotropies for color circle task.
% Rectangular variability in criterion.
%   [ll,qaic,qbic,Pred,Gstuff] =  paicircle11(Pvar, Pfix, Sel, Data, trace)
%    P = [v1...v3b, eta1a....eta3b, a, Ter, b1...b3, alpha, a1...a3, sa, pi1]
%          1...6         7...12    13, 14,  15..17     18   19...21  22  23
%  'Data' is a 3-element cell array
%  5-category version (S1 and S2)
%  Overdispersion set internally.
%
% ===========================================================================

   name = 'PAICIRCLE13: ';
   errmg1 = 'Incorrect number of parameters for model, exiting...';
   errmg2 = 'Incorrect length selector vector, exiting...';
   errmg3 = 'Data should be a 1 x 2 cell array from <makelike>...';


   np = 27;
   Overdispersion = [4.01, 1.75, 1.85, 1.40]; % S2, S3, S4, S1
   tau2 = Overdispersion(2); % S3

   % Number of trials in each discriminability condition. (SB a bit unbalanced).
   nlow = length(Data{1});
   nmed = length(Data{2});
   nhi = length(Data{3});

   epsx = 1e-9;

   if nargin < 5
       trace = 0;
   end;
   lp = length(Pvar) + length(Pfix);
   if lp ~= np
        [name, errmg1], length(Pvar), length(Pfix), return;
   end
   if length(Sel) ~= np
        [name, errmg2], length(Sel), return;
   end
   if size(Data) ~= [1,2]
        [name, errmg3], size(Data), return;
   end     
   
   % Assemble parameter vector.
   P = zeros(1,np);
   P(Sel==1) = Pvar;
   P(Sel==0) = Pfix;
   Ptemp = P;
   save Ptemp Ptemp   


   % Allow for stimulus bias in nonzero second drift component.
   v1a = P(1);
   v1b = P(2);
   v2a = P(3);
   v2b = P(4);
   v3a = P(5);
   v3b = P(6);
   eta1a = P(7);
   eta1b = P(8);
   eta2a = P(9);
   eta2b = P(10);
   eta3a = P(11);
   eta3b = P(12);
   a= P(13);
   ter = P(14);
   B = P(15:19);
   alpha = P(20);
   RawBias = P(21:25);   % a
   sa = P(26);
   pi1 = P(27);
   sigma =1.0; 

   % -----------------------------------------------------------------------------------------
   %    v1a...v3b,     eta1a....eta3b,    a,   Ter,  b1...b5, alpha, a1..a5,    sa     pi1]
   % ---------------------------------------------------- ------------------------------------

    Ub= [ 7.0*ones(1,6),  4.0*ones(1,6),  5.0,  1.0,  5.0*ones(1,6),  2*pi*ones(1,5),  2*a, 1.0]; 
    Lb= [-7.0*ones(1,6),  0.0*ones(1,6),  0.5,  0.1,    0*ones(1,6)  0*ones(1,5), 0, 0 ];
    Pub=[ 6.5*ones(1,6),  3.5*ones(1,6),  4.8,  0.8   4.5*ones(1,6)  (2*pi-eps)*ones(1,5), 2*a - eps, .95];
    Plb=[-6.5*ones(1,6),  0.0*ones(1,6),  0.7,  0.15  0.02*ones(1,6)  eps*ones(1,5), 0, .02];
    Pred = cell(3,3);
    if any(P - Ub > 0) | any(Lb - P > 0)
       ll = 1e7 + ...
            1e3 * (sum(max(P - Ub, 0).^2) + sum(max(Ub - P).^2));
       bic = 0;
       if trace
          max(P - Ub, 0)
          max(Lb - P, 0)
       end
   else
      penalty =  1e3 * (sum(max(P - Pub, 0).^2) + sum(max(Plb - P, 0).^2));
      if trace
          max(P - Pub, 0)
          max(Plb - P, 0)
          penalty
      end;   
      Pa = [v1a, v1b, B, eta1a, eta1b, sigma, a, ter, alpha, sa, pi1];
      % Parameters for med
      Pb = [v2a, v2b, B, eta2a, eta2b, sigma, a, ter, alpha, sa, 0];
      % Parameters for hi
      Pc = [v3a, v3b, B, eta3a, eta3b, sigma, a, ter, alpha, sa, 0];

      [ta, gtma, ftma, thetaa, pthetaa, mthetaa, mdthetaa, ethetaa, lla] = vacircle12j(Pa, Data{1}, RawBias);
      [tb, gtmb, ftmb, thetab, pthetab, mthetab, mdthetab, ethetab, llb] = vacircle12j(Pb, Data{2}, RawBias);
      [tc, gtmc, ftmc, thetac, pthetac, mthetac, mdthetac, ethetac, llc] = vacircle12j(Pc, Data{3}, RawBias);
     % Pass out the raw densities for the quantile-probability plot.
      Gstuff = cell(3,3);
      Gstuff{1,1} = ta;
      Gstuff{2,1} = thetaa;
      Gstuff{3,1} = gtma;
      Gstuff{1,2} = tb;
      Gstuff{2,2} = thetab;
      Gstuff{3,2} = gtmb;
      Gstuff{1,3} = tc;
      Gstuff{2,3} = thetac;
      Gstuff{3,3} = gtmc;
      % Minimize sum of minus LL's across two conditions.
     ll = sum(lla) + sum(llb) + sum(llc);
     qaic = 2 * ll /tau2  + 2 * sum(Sel); 
     qbic = 2 * ll/tau2 + sum(Sel) * log(nlow + nmed + nhi);
     
     % ftm is a double marginalization across w and v.
     Pgta = [ta; ftma];
     Pgtb = [tb; ftmb];
     Pgtc = [tc; ftmc];
     
     Ptha = [thetaa; pthetaa'];
     Pthb = [thetab; pthetab'];
     Pthc = [thetac; pthetac'];
  
     Rtha = [thetaa; ethetaa; mthetaa; mdthetaa];
     Rthb = [thetab; ethetab; mthetab; mdthetab];
     Rthc = [thetac; ethetac; mthetac; mdthetac];

     Pred{1,1} = Pgta;
     Pred{1,2} = Pgtb;
     Pred{1,3} = Pgtc;
     Pred{2,1} = Ptha;
     Pred{2,2} = Pthb;
     Pred{2,3} = Pthb;
     Pred{3,1} = Rtha;
     Pred{3,2} = Rthb;
     Pred{3,3} = Rthc;
    % Penalize log-likelihood quadratically.
     ll = abs(ll) + penalty;
  end
end

function [t, gtm, ftm, theta, pthetam, mtheta, mdtheta, etheta, ll0] = vacircle12j(Pj, Dataj, RawBias)
% ===============================================================================================
% Compute predictions and log-likelihoods for one condition.
% Initially return predictions marginalized across stimulus angle.
% Now puts everything into canonical orientation to avoid the etas problem
%      [t, gtm, ftm, theta, pthetam, mtheta, etheta, ll0] = vacircle12j(Pj, Dataj, RawBias, sz, nw, nv)
%      P = [v1, v2, b1, b2, b3, eta1, eta2, sigma, a, ter, alpha, sa, pii]
%  Calculate components of likelihood for one discriminability condition.
%  b is the amplitute of the bias vectors
% ===============================================================================================
   tmax = 4.0;  % Set from histograms by eye (see FikeKey.txt)  
   nw = 50; 
   nv = 50; 
   nt = 300; 
   h = tmax / nt; 
   w = 2 * pi / nw;
   nvm = fix(nv / 2);
   np = 15;
   nw1 = nw + 1;
   nv1 = nv + 1;
   BiasAngle = sort(signed_angle(RawBias)); % S1 peaks (radians unsigned);

   epsx = 1e-9;
   contamden = 0.05;  % Contaminant density.
   tmax = 4.0;
   ter = Pj(length(Pj) - 3); % Have sa on the end now.
   ld = size(Dataj);
   if ld(2) ~= 3
      disp('AICIRCLE300J: Wrong size data matrix, returning...')
      size(Dataj)
      return
   end
   if length(Pj) ~= np
      disp('Wrong length parameter vector, returning...');
      np
      return
   end
   [t, gt, thetas, theta, ptheta, mtheta] = vjoint5density([Pj(1:11), Pj(13:15)], BiasAngle, tmax, nw1, nv1, nt);
   % Filter zeros  
   gt = max(gt, epsx); % [51, 300, 51] % with wrap-around.
   % Add nondecision times
   t = t + ter;
   mtheta = mtheta + ter;
   [anglerr, time, angles]=ndgrid(theta, t, thetas);

   % Interpolate in joint density to get likelihoods of each data point
   l0 = interpn(anglerr, time, angles, gt, Dataj(:,2), Dataj(:,3), Dataj(:,1), 'linear');
   Cx = isnan(l0);
   l0(Cx) = contamden;
   ll0 = -log(l0);

  % Marginals for accuracy, joint distribution and MRT
   gtm = zeros(nw, nt);  % Last index is stimulus phase
   pthetam = zeros(nw, 1);
   mthetam = zeros(nw, 1);
   etheta = zeros(1, nv);

   gtm = sum(gt(:,:,1:nv), 3) / nv;
   % Predictions for plot
   ftm = sum(gtm) * w;
   pthetam = sum(ptheta(:,1:nv), 2) / nv;
   mthetam = sum(mtheta(:,1:nv), 2) / nv;
   % Mean error computed in canonical orientation - sum across rows for each column (stimulus)
   etheta = sum(ptheta .* theta') * w;
   for j = 1:nv1
       ptheta(:,j) = circshift(ptheta(:, j), j - nvm);
       mtheta(:,j) = circshift(mtheta(:, j), j - nvm);
   end
   % Sum across rows (response angle) gives mean RT for a given stimulus (Nondecision time added above)
   mtheta = sum(ptheta .* mtheta') * w;

   % Do medians
   mdtheta = zeros(1, nv1);
   for j = 1:nv
       gt(:, :, j) = circshift(gt(:, :, j), j - nvm);
   end
   % Sum over response - weird syntax b/c cumsum returns 1 x 300 x 51, need to reduce dimension.
   ft = zeros(nt, nv1);
   ft(:,:) = cumsum(sum(gt, 1)) * w * h;
   for j = 1:nv
        mjlo = max(find(ft(:, j) <  0.5));
        mjhi = max(find(ft(:, j) <= 0.5));
        mdtheta(j) = (t(mjlo) + t(mjhi)) / 2; 
   end
   mdtheta(nv1) = mdtheta(1); 
end



function [t, gt, thetas, thetaerr, ptheta, mtheta] = vjoint5density(P, BiasAngle, tmax, nw1, nv1, nt)
% ===============================================================================================
%  [t,gt,thetas, thetaerr,ptheta,mtheta] = vjoint5density(P, BiasAngle, tmax, nw1, nv1, nt);
% P = [v1, v2, B, eta1, eta2, sigma, a, alpha, sa, pii]
% Circular diffusion predictions as a function of stimulus angle.
% Stimuli in canonical orientation (i.e., re 0), bias computed as an offset. 
% ===============================================================================================

   v1 = P(1);
   v2 = P(2);
   B = P(3:7);
   eta1 = P(8);
   eta2 = P(9);
   sigma = P(10);
   a = P(11);
   alpha = P(12);
   sa = P(13);
   pii = P(14);
   nt = 300; % number of time steps
   nbias = 5;
   badix = 5; % should be 25 for S4
   epsilon  = 0.0001;
   n_sa_step = 11;
   nv = nv1 - 1;
   v = 2 * pi / nv;

   thetaerr = linspace(-pi, pi, nw1); % To accommodate wrap around. 
   thetas = linspace(-pi, pi,  nv1);
   ltheta = length(thetas);
   thetav = zeros(1, ltheta);
   Vtheta = zeros(2, ltheta);   % Values of drift at stimulus angle.
   %VthetaPhase = zeros(1, ltheta)
   gt = zeros(nw1, nt, nv1);  % Last index is stimulus phase
   ptheta = zeros(nw1, nv1);
   mtheta = zeros(nw1, nv1);
   etheta = zeros(nw1, nv1);
   t = linspace(0, tmax, nt);

   gt_all = zeros(nw1, nt, nv1, n_sa_step);  % Last index is stimulus phase
   ptheta_all = zeros(nw1, nv1, n_sa_step);
   mtheta_all = zeros(nw1, nv1, n_sa_step);

   % Structures for parfor

   vnorm = sqrt(v1.^2 + v2.^2);
   % Bias
   Thetasex = ones(nbias,1) * thetas;
   %thetas
   %BiasAngle
   BiasAnglex = BiasAngle' * ones(1, nv1);
   Bias = B' * ones(1, nv1);
  % Use 1 - cos distance. 
   Distance = BiasAnglex - Thetasex;
   CircularDistance = 1 - cos(Distance);
   DecayedBias =  vnorm * Bias .* exp(-alpha * CircularDistance); 
   DistanceCos = cos(Distance);
   DistanceSin = sin(Distance);
   SumBiasCos = sum(DecayedBias .* DistanceCos);
   SumBiasSin = sum(DecayedBias .* DistanceSin); 

   if sa < epsilon
        for k = 1:nv1 
             Vtheta(1,k) = v1 + SumBiasCos(k);
             Vtheta(2,k) = v2 + SumBiasSin(k);        
             [~,gt(:,:,k), ~, ptheta(:,k), mtheta(:,k)] = vdcircle300cls([Vtheta(1,k), Vtheta(2,k), ...
                   eta1, eta2, sigma, a], tmax, badix);
        end
   else % Rectangular criterion variability
        U = ones(n_sa_step, 1); 
        Rmass = U / n_sa_step ; 
        Rstep = [-(n_sa_step-1)/2:(n_sa_step-1)/2]' / (n_sa_step-1); 
        A = a + Rstep * sa;
        Vtheta = zeros(2,nv,n_sa_step);
        for i = 1:n_sa_step
            for k = 1:nv
                Vtheta(1,k,i) = v1 + SumBiasCos(k);
                Vtheta(2,k,i) = v2 + SumBiasSin(k);        
             end
        end
        %A
        parfor i = 1:n_sa_step % Was parallelized
            for k = 1:nv
                [ti, gtik, thetaik, pthetaik, mthetaik] = vdcircle300cls([Vtheta(1,k,i), Vtheta(2,k,i), ...
                       eta1, eta2, sigma, A(i)], tmax, badix);
                gt_all(:,:,k,i) =  gtik;
                ptheta_all(:,k,i) = pthetaik';
                mtheta_all(:,k,i) = mthetaik';
            end   
        end
        for i = 1:n_sa_step
             gt = gt + gt_all(:,:,:, i) / n_sa_step;
             ptheta = ptheta + ptheta_all(:,:,i) / n_sa_step;
             mtheta = mtheta + mtheta_all(:,:,i) / n_sa_step;
        end
   end;
   if pii > 0
       [~, gt0, ~, ptheta0, mtheta0] = rcircle([vnorm, sigma, a], tmax, badix); % ####
       for k = 1:nv   
           gt(:,:,k) = (1 - pii) * gt(:,:,k) + pii * gt0;
           ptheta(:,k) = (1 - pii) * ptheta(:,k) + pii * ptheta0';
           mtheta(:,k) = (1 - pii) * mtheta(:,k) + pii * mtheta0'; 
       end
   end

   % Wrap around to close for interpolation.
   gt(:,:,nv1) = gt(:,:, 1);
   ptheta(:,nv1) = ptheta(:, 1);
   mtheta(:, nv1) = mtheta(:, 1);
end


function [T, Gt, Theta, Ptheta, Mtheta] = rcircle(P, tmax, badix);
% ====================================================================
% Randomize drift around the circle. Need wrap-around
% ====================================================================
    munorm = P(1);
    sigma = P(2);
    a = P(3);
    nw = 50; % hard wire here.
    nz = 300;
    h = tmax / nz;
    [T, Gt0]= dhamanax([a, sigma], h, tmax, badix);
    Ht = besseli(0,a * munorm/sigma^2) * exp(-munorm^2/(2*sigma^2) * T) .* Gt0; 
    mt = sum(T .* Ht) * h;
    Gt = ones(nw+1, 1) * Ht / (2 * pi);
    w = 2 * pi / nw;
    Theta = [];
    Ptheta = ones(1, nw+1) / (2 * pi);
    Mtheta = mt * ones(1, nw+1);
    %mass = sum(Ptheta * w)
end


function [T, Gt] = dhamanax(P, h, tmax, badix)
% ====================================================================
% First passage time density for a bessel process, 
% Derivative of Hamana-Matsumoto solution, for x0 = 0 (Eq. 2.7)
% [t, gt] = dhamana(P, kmax, h, tmax, badix);
% P =[a, sigma];
% b is boundary, a is starting point.
% kmax controls truncation of series
% badix is number of bad initial values in Gt to zero.
% =====================================================================
if nargin < 4
    badix = 0;
end    
a = P(1);
sigma = P(2);
sigma2 = sigma^2;
%T = [0:h:tmax];
T = linspace(0, tmax, tmax/h);
a2 = a^2;
%J0k = besselzero(0, kmax, 1) % v = 0 for drift = (2v + 1)/2x 
J0k =[2.4048    5.5201    8.6537   11.7915   14.9309   18.0711   21.2116   24.3525   27.4935   30.6346 ...
     33.7758   36.9171   40.0584   43.1998   46.3412   49.4826   52.6241   55.7655   58.9070   62.0485 ...
     65.1900   68.3315   71.4730   74.6145   77.7560   80.8976   84.0391   87.1806   90.3222   93.4637 ...
     96.6053   99.7468  102.8884  106.0299  109.1715  112.3131  115.4546  118.5962  121.7377  124.8793 ...
     128.0209  131.1624  134.3040  137.4456  140.5872  143.7287  146.8703  150.0119  153.1535  156.2950];
J0k_squared = J0k.^2;
kmax = length(J0k);
J1k = besselj(1, J0k);
Rt = zeros(size(T));
scaler = sigma2 / a2;
for k = 1:kmax
    Rt = Rt + J0k(k) * exp(-J0k_squared(k) * sigma^2 * T /(2 * a2)) ...
         / J1k(k);
end;
Gt = scaler * Rt;
if badix > 0 
   Gt(1:badix) = 0;
end
Gt = max(Gt, 0);   
end
    

  
function sangle = signed_angle(a)
% ================================================================
% Convert angles on [0 : 2 *pi] to [-pi : +pi]
% ================================================================
    a = a / pi;
    sangle = pi * (rem(a , 1) - fix(a / 1));
end


