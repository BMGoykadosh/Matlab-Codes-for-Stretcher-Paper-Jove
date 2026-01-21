function a = myfilter_exp(x,f0);

% This is a filter for displacements. Use it on the displacement matrices,
% write e.g. "new_displacements = myfilter_exp(old_displacements,
% frequency)".
% x is a matrix of displacements (either dx or dy).
% f0 is the cut-off frequency, usually 15 or so.

nr = size(x,1);
nc = size(x,2);

xf = fftshift(fft2(x));
[x1, y1] = meshgrid([-nc/2:-1, 0:nc/2-1], [-nr/2:-1, 0:nr/2-1]);
r = sqrt(x1.^2 + y1.^2);

factor = 1./(1+exp(1*(r-f0)));

xf = xf.*factor;
a = real(ifft2(fftshift(xf))); % filtered displacements