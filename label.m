function label(axhandle, xfrac, yfrac, labeltext, font);
% ===============================================================================
% Places text on specified axes at (xmin + xfrac, ymin + yfrac)
% Optional font parameter.
% Usage:
%        label(axhandle, xfrac, yfrac, labeltext, font);
% =============================================================================
axes(axhandle);
xl = get(gca, 'XLim');
yl = get(gca, 'YLim');
x = xl(1) + xfrac * (xl(2) - xl(1));
y = yl(1) + yfrac * (yl(2) - yl(1));  
if isstr(labeltext)
    texthandle = text(x, y, labeltext);
else
    disp('LABEL: variable "labeltext" must be text, exiting...');
end;
if nargin > 4
   if isstr(font)
        set(texthandle, 'FontName', font);
   else
        disp('LABEL: variable "font" must be text, exiting...');
   end;
end; 
