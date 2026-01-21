clear vars; close all; clc;
%Rigid body shift the images first
%This code can be used to calculate the discplacement and strain. You need
%to set up the

pixel_ratio = 0.32594;  % microns per pixel 

basepath='D:\Benjie_storage\Fall_2024\stretcher_experiments\11_1_second_test_stretcher\RBS_Corrected_Uniaxial_stretch_test';
process_folder = '\hello_world_1\'; %%Name and create folder to store your data. 
mkdir(fullfile(strcat(basepath,process_folder)));
processed_path = fullfile(basepath,process_folder,'\');
for i=1 %7:160

stretch_frame=sprintf('Uniaxial_1_step_10_A.tif',i-1); %Stretched image
stretch_file = fullfile(basepath,stretch_frame);
stretch_img =double(imread(stretch_file));


cell_no=i-1;

base_frame='Uniaxial_1_step_10_B.tif'; %Unstretched image
base_file = fullfile(basepath,base_frame);
base_img =double(imread(base_file));

[xv,yv,uxv,uyv]=Sam_beads_imcorr_v2_stretcher(base_img, stretch_img); % UPDATE block size, bead intensity in the function. Explanation in function

xv = xv.*pixel_ratio;        yv = yv.*pixel_ratio;
uxv = uxv.*pixel_ratio;    uyv = uyv.*pixel_ratio;
N =sqrt(length(xv));         box =size(stretch_img,1)/N;
spacing = box*pixel_ratio; 
x = reshape(xv,N,N); 
y = reshape(yv,N,N);
ux = reshape(uxv,N,N); 
uy = reshape(uyv,N,N); 

figure;
surf(x,y,sqrt(ux.^2+uy.^2)*pixel_ratio); set(gcf,'Renderer','zbuffer');hold on;
view(2); shading interp; colorbar; % caxis([0 0.5])
h=quiver3(x,y,100*ones(size(x)),ux,uy,zeros(size(ux)),1.5);
% plot3(xb,yb,100*ones(length(xb),1),'w-');
set(h,'Color','w','LineWidth',1);
title('Displacements','FontSize',18,'FontWeight','bold');axis ij;axis image;
%axis([10 600 10 600])
opfile =fullfile(strcat(processed_path),sprintf('disp_map_Frame_%03d.png',cell_no));
print(opfile,'-dpng','-r0')

E=strain(ux,uy);
EXY =squeeze(E(:,:,2,1)); 
EXX =squeeze(E(:,:,1,1)); 
EYY =squeeze(E(:,:,2,2));
E1 =sqrt(EXX.^2+EYY.^2);

midE1=E1(4:12,4:12); % selects the middle of the strain field, need to change indices depending on block size
meanE=mean(E1(:));
meanMidE=mean(midE1(:));
medianE=median(E1(:));
medianMidE=median(midE1(:));
stdE=std(E1(:));
stdMidE=std(midE1(:));

figure;
surf(x,y,E1*pixel_ratio); set(gcf,'Renderer','zbuffer');%hold on;
view(2); shading interp; colorbar; caxis([0 0.2]);
% h=quiver3(x,y,100*ones(size(x)),EXX,EYY,zeros(size(EXX)),1.5);
% plot3(xb,yb,100*ones(length(xb),1),'w-');
set(h,'Color','w','LineWidth',1);
title('Strain','FontSize',18,'FontWeight','bold');axis ij;axis image;
opfile =fullfile(strcat(basepath,process_folder),sprintf('Strain_Frame_%03d.png',cell_no));
print(opfile,'-dpng','-r0')

opfile =fullfile(strcat(basepath,process_folder),sprintf('data_frame_%03d.mat',cell_no));
save(opfile);

close all

end
     