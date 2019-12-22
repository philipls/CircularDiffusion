function angplotmd(Raw, Pred, ymax)
% ====================================================================================
% Plot the marginal mean error and median RT as a function of stimulus angle. Predictions
% are as output by <aicircle9b> etc. Raw are the raw data files (ntrials x 8 columns)
%     angplotmd(Raw, Pred, {ymax})
% ====================================================================================
   co = [0,0,1; 
         0,0.5,0;
         1,0,0;
         0,0.75,0.75;
         0.75,0,0.75;
         0.75,0.75,0;
         0.25,0.25,0.25];

  set(groot,'defaultAxesColorOrder',co);

   err1 = 'ANGPLOTMD: Wrong size raw data';
   err2 = 'ANGPLOTMD: Wrong size Pred';
   c = 0.50;
   %cvec = [.7, .7, .7];
   %cvec = [.75, .85, .85];   % Kinda blue...
   cvec1 = [.90, .95, .95];   % Kinda blue...
   cvec2 = [.60, 0, .60];   % Dark magenta
   mdx = 4; % Row index of median in data structure
   szr = size(Raw)
   if szr(2) ~= 8
      err1
      return
   end
   if ~all(size(Pred) == [3,3])
      err2 
      return
   end
   if nargin < 3
     ymax = 2.0;
   end
   axhandle = setfig6;
   
   [theta1,thetacc1,thetart1]=thetaclass4(Raw, 50, 1);
   axes(axhandle(1));
   bar(theta1, thetacc1)
   hold
   plot(Pred{3,1}(1,:), Pred{3,1}(2,:), 'm-', 'Linewidth', 2)
   ch = get(gca, 'Child');
   ch(2).FaceColor = cvec1;
   ch(1).Color = cvec2;
   set(gca, 'YLim', [-1,1])
   set(gca, 'XLim', [-pi,pi])
   xlabel('Stimulus angle (rad)')
   ylabel('Mean Error (rad)');
   label(gca, .1, .9, 'High')

   axes(axhandle(2));
   bar(theta1, thetart1 - c)
   hold
   plot(Pred{3,1}(1,:), Pred{3,1}(mdx,:) - c, '-', 'Linewidth', 2)
   ch = get(gca, 'Child');
   ch(2).FaceColor = cvec1;
   ch(1).Color = cvec2;
   set(gca, 'YLim', [0,ymax])
   set(gca, 'XLim', [-pi,pi])
   xlabel('Stimulus angle (rad)')
   ylabel('MdRT - 0.5')
   label(gca, .1, .9, 'High')

   [theta2,thetacc2,thetart2]=thetaclass4(Raw, 50, 2);
   axes(axhandle(3));
   bar(theta2, thetacc2)
   hold
   plot(Pred{3,2}(1,:), Pred{3,2}(2,:), '-', 'Linewidth', 2)
   ch = get(gca, 'Child');
   ch(2).FaceColor = cvec1;
   ch(1).Color = cvec2;
   set(gca, 'YLim', [-1,1])
   set(gca, 'XLim', [-pi,pi])
   xlabel('Stimulus angle (rad)')
   ylabel('Mean Error (rad)');
   label(gca, .1, .9, 'Med')

   axes(axhandle(4));
   bar(theta2, thetart2 - c)
   hold
   plot(Pred{3,2}(1,:), Pred{3,2}(mdx,:) - c, '-', 'Linewidth', 2)
   ch = get(gca, 'Child');
   ch(2).FaceColor = cvec1;
   ch(1).Color = cvec2;
   set(gca, 'YLim', [0,ymax])
   set(gca, 'XLim', [-pi,pi])
   xlabel('Stimulus angle (rad)')
   ylabel('MdRT - 0.5')
   label(gca, .1, .9, 'Med')

   [theta3,thetacc3,thetart3]=thetaclass4(Raw, 50, 3);
   axes(axhandle(5));
   bar(theta3, thetacc3)
   hold
   plot(Pred{3,3}(1,:), Pred{3,3}(2,:), '-', 'Linewidth', 2)
   ch = get(gca, 'Child');
   ch(2).FaceColor = cvec1;
   ch(1).Color = cvec2;
   set(gca, 'YLim', [-1,1])
   set(gca, 'XLim', [-pi,pi])
   xlabel('Stimulus angle (rad)')
   ylabel('Mean Error (rad)');
   label(gca, .1, .9, 'Low')

   axes(axhandle(6));
   bar(theta3, thetart3 - c)
   hold
   plot(Pred{3,3}(1,:), Pred{3,3}(mdx,:) - c, 'm', 'Linewidth', 2)
   ch = get(gca, 'Child');
   ch(2).FaceColor = cvec1;
   ch(1).Color = cvec2;
   set(gca, 'YLim', [0,ymax])
   set(gca, 'XLim', [-pi,pi])
   xlabel('Stimulus angle (rad)')
   ylabel('MdRT - 0.5')
   label(gca, .1, .9, 'Low')


