function [cx] = center_x_1d(x);

% CENTER_X_1D calculates the x-position of the peak of the 
%	cross-correlation function x by fitting the nearest 
% 	neighbors (left and right) of the maximum value of x to a parabola.

try;
   
[mi,mj]	= find(x == max(max(x)));
mi			= mi(1);
mj			= mj(1); 

[xobj, yobj]	= 	meshgrid(1 : size(x,2), 1 : size(x,1));

x_fit = xobj(mi, mj-1 : mj+1);
y_fit = x(mi, mj-1 : mj+1);

[a,b] = polyfit(x_fit,y_fit,2);

x_vals = [min(x_fit) : 0.01 : max(x_fit)];
y_vals = polyval(a, x_vals);

cx = x_vals(find(y_vals == max(y_vals))) - size(x,2)/2 - 1;

if length(cx) > 1, cx = mean(cx); end;

catch;
   cx = NaN;
end;
