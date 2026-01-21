function [p1new,p2new,dx_after,dy_after] = move_blocks(p1,p2,dx,dy,suspect_i,suspect_j);

% MOVE_BLOCKS creates images p1new and p2new from the images p1 and p2
%	such that at each position in p2new there is a block from p2 that is the
%	the best match to the block in p1new at the same position.

%======================================================================================

if nargin <= 4,
   suspect_i = [];
   suspect_j = [];
end;
    
    dx(find(~isfinite(dx))) = 0;
	dy(find(~isfinite(dy))) = 0;
    dx = round(dx); 
    dy = round(dy);

blocksize = size(p2,1) / size(dy,1);

[xwhole,ywhole] = meshgrid(1 : size(p2,2), 1 : size(p2,1));

%======================================================================================
% MOVE BLOCKS IN p2

	x_beg = blkproc(xwhole,[blocksize blocksize],'min2(x)');
	x_end = blkproc(xwhole,[blocksize blocksize],'max2(x)');
	y_beg = blkproc(ywhole,[blocksize blocksize],'min2(x)');
	y_end = blkproc(ywhole,[blocksize blocksize],'max2(x)');

% CORNERS IN p1
	corners_old(:,:,1) = y_beg;
	corners_old(:,:,2) = y_end;
	corners_old(:,:,3) = x_beg;
	corners_old(:,:,4) = x_end;

% CORNERS IN p2 TO COMPARE WITH p1
	corners(:,:,1) = y_beg + dy;
	corners(:,:,2) = y_end + dy;
	corners(:,:,3) = x_beg + dx;
	corners(:,:,4) = x_end + dx;

% CORNERS THAT GO OUTSIDE OF p2 HAVE TO GET INSIDE
	[fc1i, fc1j] = find(corners(:,:,1) < 1);
	corners(fc1i, fc1j, 1) = 1;
	corners(fc1i, fc1j, 2) = blocksize;
	[fc3i, fc3j] = find(corners(:,:,3) < 1);
	corners(fc3i, fc3j, 3) = 1;
	corners(fc3i, fc3j, 4) = blocksize;

	[fc2i, fc2j] = find(corners(:,:,2) > size(p2,1));
	corners(fc2i, fc2j, 1) = corners(fc2i, fc2j, 1) - (corners(fc2i, fc2j, 2) - size(p2,1));
	corners(fc2i, fc2j, 2) = size(p2,1);
	[fc4i, fc4j] = find(corners(:,:,4) > size(p2,2));
	corners(fc4i, fc4j, 3) = corners(fc4i, fc4j, 3) - (corners(fc4i, fc4j, 4) - size(p2,2));
	corners(fc4i, fc4j, 4) = size(p2,2);

% FINAL CORNERS IN p2
	[in, jn, kn] = find(~isfinite(corners));
	corners(in,jn,:) = corners_old(in,jn,:);
   corners_new = corners;
   
%======================================================================================
% CALCULATE dx_after, dy_after

dx_after = corners_new(:,:,3) - corners_old(:,:,3);
dy_after = corners_new(:,:,1) - corners_old(:,:,1);

%======================================================================================
% MAKE NEW p2

newpict = [];
oldpict = [];

if ~isempty(suspect_i),
   corners_old1 = [];
   corners_new1 = [];
	for i = 1 : length(suspect_i),
   	corners_old1 = [corners_old1; corners_old(suspect_i(i),suspect_j(i),:)];
   	corners_new1 = [corners_new1; corners_new(suspect_i(i),suspect_j(i),:)];
   end; %(for i)
   corners_old = corners_old1;
   corners_new = corners_new1;   
end; %(isempty)

for j = 1 : size(corners_new,2),
   for i = 1 : size(corners_new,1),
 	  	oldpict = [oldpict; p1(corners_old(i,j,1):corners_old(i,j,2),...
 	        corners_old(i,j,3):corners_old(i,j,4))];
 	  	newpict = [newpict; p2(corners_new(i,j,1):corners_new(i,j,2),...
 	        corners_new(i,j,3):corners_new(i,j,4))];
   end; %(for i)
end; %(for j)

p1new = oldpict;
p2new = newpict;

%======================================================================================