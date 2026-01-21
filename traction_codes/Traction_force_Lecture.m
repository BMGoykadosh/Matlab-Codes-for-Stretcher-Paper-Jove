close all;clear all;

pixel_ratio =0.652;%microns per pixel for 10x   %0.16125; % microns per pixel
%box =50 ;% pixels-- grid size for displacements
pois =0.445;
young =6705; %Pa
young =young*1e-12;
f0 =15;



bead_cell_img =imread(fullfile(pwd,'psanv2','scene00200.png'));
% bead_cell_img =imread(fullfile(pwd,'Data','BF','BF0000.tif'));
bead_trp_img =imread(fullfile(pwd,'psanv2','scene00400.png'));
% bead_trp_img =imread(fullfile(pwd,'Data','BF','BF9999.tif'));
figure; imshow(bead_cell_img); colormap(gray);axis image;  title('Baseline');
figure; imshow(bead_trp_img); colormap(gray); axis image; title('A');
%%
%%--- Crop a small border around the edges ----%%
shx =0; shy =0; bead_cell_img2 =  bead_cell_img(50+shy:200+100+shy, 50+shx:200+100+shx); % old numbers 200 711 530 1041
shx =0; shy =0; bead_trp_img2  =  bead_trp_img(50+shy:200+100+shy,  50+shx:200+100+shx); % old numbers 200 711 530 1041
% figure; imagesc(bead_cell_img2); colormap(gray); title('Baseline with cell-CROPPED');
% figure; imagesc(bead_trp_img2); colormap(gray); title('After trypsin-CROPPED');

%%
%---Calculate rigid body displacement---%
[shx,shy]=disp_on_blocks_v2(bead_cell_img2, bead_trp_img2,size(bead_cell_img2,1),0);
bead_trp_img2  =  bead_trp_img(50+shy:200+100+shy,  50+shx:200+100+shx); % old numbers 200 711 530 1041
figure; imshow(bead_cell_img2); colormap(gray); title('Baseline-RIGID BODY SHIFT CORRECTED');
figure; imshow(bead_trp_img2); colormap(gray); title('A-RIGID BODY SHIFT CORRECTED');

%%
%----Calculate displacement_map------%
[xv,yv,uxv,uyv]=beads_imcorr_v2(bead_trp_img2, bead_cell_img2);
% convert from pixels to microns
xv = xv.*pixel_ratio;
yv = yv.*pixel_ratio;
uxv = uxv.*pixel_ratio;
uyv = uyv.*pixel_ratio;


N =sqrt(length(xv));box =size(bead_cell_img2,1)/N;
spacing = box*pixel_ratio; % box =32;

x = reshape(xv,N,N); %x = reshape(x_temp,minsize^2,1);
y = reshape(yv,N,N); %y = reshape(y_temp,minsize^2,1);
ux = reshape(uxv,N,N); %dx = reshape(dx_temp,minsize^2,1);
uy = reshape(uyv,N,N); %dy = reshape(dy_temp,minsize^2,1);


figure;
surf(x,y,sqrt(ux.^2+uy.^2)*pixel_ratio); set(gcf,'Renderer','zbuffer');hold on;
view(2); shading interp; colorbar; colormap(fire)%caxis([0 6]);
h=quiver3(x,y,100*ones(size(x)),ux,uy,zeros(size(ux)),1);
% plot3(xb,yb,100*ones(length(xb),1),'w-');
set(h,'Color','w','LineWidth',1);
title('Representative Displacement (\mum)','FontSize',18,'FontWeight','bold');axis ij;
axis image;

% Calculate strain
E=strain(ux,uy);
AvgStrain=mean(E(:))


%%
wf_cell_img =imread(fullfile(pwd,'Data2','FL','Image041.tif'));
wf_cell_img  =  wf_cell_img(1000+shy:1511+shy,  1000+shx:1511+shx); % old numbers 200 711 530 1041
figure; imshow(wf_cell_img,[],'Xdata',xv,'Ydata',yv);  hold on
h=quiver3(x,y,100*ones(size(x)),ux,uy,zeros(size(ux)),1);
% plot3(xb,yb,100*ones(length(xb),1),'w-');
set(h,'Color','r','LineWidth',1.5);
title('Displacements','FontSize',18,'FontWeight','bold');axis ij;axis image;
%%
%======================================================================

% Filter for the displacement matrix

ux  = myfilter_exp(ux, f0); % f0 is the cut-off frequency.
uy  = myfilter_exp(uy, f0);

% UNCONSTRAINED TRACTIONS
a1 = (1.0 + pois) * (1.0 - pois) / (pi * young);
b1 = (1.0 + pois) * pois / (pi * young);
c1 = (1.0 + pois) * pois / (pi * young);
% constants for unconstrained tractions
clear sx sy;
for i = 1:(N/2)
    sx(i,:) = 0:((N/2)-1);
    sy(i,1:(N/2)) = (N/2)-i;
end;

kx = [ sx  sx-(N/2);  sx  sx-(N/2) ];
ky = [ sy-(N/2)  sy-(N/2);  sy  sy ];
ky = flipud(ky);
kx(:,(N/2+1)) =  kx(:,(N/2+1));
ky((N/2+1),:) =  ky((N/2+1),:);
k_abs = sqrt(kx.^2 + ky.^2);

alpha = atan2(ky,kx);
if kx(1,1) == 0 && ky(1,1) == 0,
    alpha(1,1) = pi/2;
end;

Cx = ((k_abs * young) / (2 * (1 - pois^2))) .* (1 - pois + pois .* (cos(alpha)).^2);
Cy = ((k_abs * young) / (2 * (1 - pois^2))) .* (1 - pois + pois .* (sin(alpha)).^2);
D  = ((k_abs * young) / (2 * (1 - pois^2))) .* (pois .* sin(alpha) .* cos(alpha));

D(:,(N/2+1)) = zeros(N,1);
D((N/2+1),:) = zeros(1,N);

% Calculate Unconstrained Tractions
Dx = fft2(2*pi*ux / (N*spacing));
Dy = fft2(2*pi*uy / (N*spacing));

Tx = Cx.*Dx + D.*Dy;
Ty = D.*Dx + Cy.*Dy;

tx_uc = real(ifft2(Tx));
ty_uc = real(ifft2(Ty));

% Calculate Contractile Moment
Dxx = imag(Tx(1,2)) * (N * spacing / (2*pi)) * spacing^2 * 1e-6;
Dyy = imag(Ty(2,1)) * (N * spacing / (2*pi)) * spacing^2 * 1e-6;
Dxy = imag(Ty(1,2)) * (N * spacing / (2*pi)) * spacing^2 * 1e-6;
Dyx = imag(Tx(2,1)) * (N * spacing / (2*pi)) * spacing^2 * 1e-6;

CM_uc = [Dxx Dxy; Dyx Dyy];


% Plot Unconstrained Tractions
figure; hold on;
surf(x,y,sqrt(tx_uc.^2+ty_uc.^2)); set(gcf,'Renderer','zbuffer');
view(2); shading interp; colorbar;
h=quiver3(x,y,100*ones(size(x)),tx_uc,ty_uc,zeros(size(tx_uc)),1);
%plot3(xb,yb,100*ones(length(xb),1),'w-');
set(h,'Color','w','LineWidth',1);
axis ij; axis([min2(x) max2(x) min2(y) max2(y)]); axis square;
title('Unconstrained Tractions','FontSize',18,'FontWeight','bold');

disp(' Contractile moment');


