function [cy] = center_y_1d(x);

% CENTER_Y_1D calculates the y-position of the peak of the 
%	cross-correlation function x by fitting the nearest 
% 	neighbors (up and down) of the maximum value of x to a parabola.

try;
   
[mi,mj]	= find(x == max(max(x)));
mi			= mi(1);
mj			= mj(1);

[xobj, yobj]	= 	meshgrid(1 : size(x,2), 1 : size(x,1));

x_fit = yobj(mi-1 : mi+1, mj);
y_fit = x(mi-1 : mi+1, mj);

[a,b] = polyfit(x_fit,y_fit,2);

x_vals = [min(x_fit) : 0.01 : max(x_fit)];
y_vals = polyval(a, x_vals);

cy = x_vals(find(y_vals == max(y_vals))) - size(x,1)/2 - 1;

if length(cy) > 1, cy = mean(cy); end;

catch;
   cy = NaN;
end;
