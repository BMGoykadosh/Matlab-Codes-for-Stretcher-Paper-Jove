function  bx = max_cor_bx(x);

[by, bx] = find(x == max(max(x))); 

bx = bx(1) - size(x,2)/2 - 1;
%by = by(1) - size(x,1)/2 - 1;

if length(bx) > 1,
   bx = 0;
end;
