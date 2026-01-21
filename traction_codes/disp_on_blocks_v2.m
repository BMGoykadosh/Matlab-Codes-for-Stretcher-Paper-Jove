function [bx,by,corsh] = disp_on_blocks_v2(im1,im2,blocksize,subpix)

% DISP_ON_BLOCKS calculates x- and y-displacements for each distinct block 
%	taken from the image im1. The size of the blocks is determind by the 
%	parameter blocksize. For each block in the image im1 the corresponding
%	block in the image im2 is found (i.e., the block at the same position),
%	and the cross-correlation function between the two blocks is formed.
%	Coordinates of the peak of the cross-correlation function constitute
%	the displacement vector of the block. (Displacements go from im1 to im2.) 
%	The parameter subpix determines whether the resulting values of displacements 
%	are integers or not; if subpix > 0, the values of displacements are 
%	non-integers, otherwise they are integers.

%	bx is a matrix of x-displacements, by a matrix of y-displacements, and
%	corsh the fft-shifted cross-correlation by block between im1 and im2.

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

if nargin < 4,
   subpix	= 	0;
end;

% SUBTRACT THE MEAN OF EACH BLOCK
fun = @(block_struct) block_struct.data -mean2(block_struct.data);%%-- change deprecated function blockproc to blockproc
im1 = blockproc(im1,[blocksize blocksize],fun);
im2 = blockproc(im2,[blocksize blocksize],fun);

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% CALCULATE THE CROSS-CORRELATION FUNCTION
fun = @(block_struct) ifft2((fft2(block_struct.data)).*conj(fft2(block_struct.data)));
autocor_whole1 = blockproc(im1,[blocksize blocksize],fun);
autocor_whole2 = blockproc(im2,[blocksize blocksize],fun);

fun = @(block_struct) repmat(block_struct.data(1,1),size(block_struct.data));
autocor_max1   = blockproc(autocor_whole1,[blocksize blocksize],fun);
autocor_max2   = blockproc(autocor_whole2,[blocksize blocksize],fun);

fun = @(block_struct) conj(fft2(block_struct.data));
cor_whole1     = blockproc(im1,[blocksize blocksize],fun);
fun = @(block_struct) fft2(block_struct.data); 
cor_whole2     = blockproc(im2,[blocksize blocksize],fun);

cor_whole      = cor_whole2 .* cor_whole1;
fun = @(block_struct) ifft2(block_struct.data); 
cor_whole		= blockproc(cor_whole,[blocksize blocksize],fun);
cor_whole      = real(cor_whole) ./ sqrt(autocor_max1.*autocor_max2);
cor_whole(find(~isfinite(cor_whole))) = zeros(size(find(~isfinite(cor_whole))));
fun = @(block_struct) fftshift(block_struct.data); 
corsh          = blockproc(cor_whole,[blocksize blocksize],fun);


%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% FIND THE PEAK OF THE CROSS-CORRELATION FUNCTION
funx = @(block_struct)center_x_1d(block_struct.data);
funy = @(block_struct)center_y_1d(block_struct.data);

funcorrx = @(block_struct)max_cor_bx(block_struct.data);
funcorry = @(block_struct)max_cor_by(block_struct.data);

if subpix > 0
    %disp('entered wrong loop')
	bx				= blockproc(corsh,[blocksize blocksize],funx);
	by				= blockproc(corsh,[blocksize blocksize],funy);
else   
	bx				= blockproc(corsh,[blocksize blocksize],funcorrx);
	by				= blockproc(corsh,[blocksize blocksize],funcorry);
end %(if subpix > 0)





