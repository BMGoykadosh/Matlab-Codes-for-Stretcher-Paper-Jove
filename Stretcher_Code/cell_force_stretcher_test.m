close all; clear
%Code to get the average traction forces
%This would only give you an average of the traction forces. It can be
%useful for stretcher experiments but does not serve to get strain
%However, it is used to get the tractions and change in tractions of the
%cells using traction force microscopy 
dir_eval = 'C:\Users\Caroline McCormick\Desktop\Benjie\Fall_2024\9_24_first_test_stretcher\2024_09_24_16_19_38--uniaxial_test_n1_.2um_10x_13kpa\uniaxial_stretch\processed_32block_1000_thresh_to_Tryp\';
dir_save = [dir_eval,'results\']; % folder to save results

path = 'C:\Users\Caroline McCormick\Desktop\Benjie\Fall_2024\9_24_first_test_stretcher\2024_09_24_16_19_38--uniaxial_test_n1_.2um_10x_13kpa\uniaxial_stretch\processed_32block_1000_thresh_to_Tryp\';%'C:\Users\sstas\Downloads\MSM Project Files\MSM Project Files\Dataset 3\Traction Data\Traction3';
ii=1;
jj=1;
n_frames = 1; 
numblocks = 60; %This relates to the blocksize chosen in the "sams beads" function. It is pixel size / bin number. Ex:1920x1920 pixel size, 32 bins -> 1920/32 = 60
xloc=zeros(numblocks^2,n_frames);
yloc=zeros(numblocks^2,n_frames);
ux_final = zeros(numblocks^2,n_frames);
uy_final = zeros(numblocks^2,n_frames);
tx_final = zeros(numblocks^2,n_frames);
ty_final = zeros(numblocks^2,n_frames);
steps_time = 1;
meanmag = zeros(steps_time,1);

for iii = 1:steps_time
    load([path,sprintf('cell_forces_t%02d.mat',iii)]);
    tempxloc = x(:);
    tempyloc = y(:);
    tempux = ux(:);
    tempuy = uy(:);
    temptx = tx_uc(:);
    tempty = ty_uc(:);
    xloc(:,iii) = tempxloc;
    yloc(:,iii) = tempyloc;
    uxf1 = ux(:);
    uxf1 = ux(:);
    uyf1 = uy(:);
    txf1 = tx_uc(:);
    tyf1 = ty_uc(:);
    ux_final(:,iii) = uxf1;%+tempux;
    uy_final(:,iii) = uyf1;%+tempuy;
    tx_final(:,iii) = txf1;%+temptx;
    ty_final(:,iii) = tyf1;%+tempty;
    %figure; surf(sqrt((tx_uc+txmat1).^2+(ty_uc+tymat1).^2),'EdgeColor','none');view(2);colorbar;
    %pause(1);
    %close all;   
    %figure; surf(sqrt((ux+uxmat1).^2+(uy+uymat1).^2),'EdgeColor','none');view(2);colorbar;
    %pause(1);
    %close all;   
    tempmag = sqrt((tx_uc).^2+(ty_uc).^2);
    meanmag(iii,1) = mean(mean(tempmag(10:50,10:50))); %this takes an average of the traction in the center of the image. 
    %tempmagdisp = sqrt((ux+uxmat1).^2 +(uy+uymat1).^2;
    %meanmagdisp(iii,1) = mean(mean(tempmagdisp(13:20, 8:26)));
end


time_step_label_1 = [1:steps_time];
time_step_label = [time_step_label_1 * 2];
time_2 = time_step_label';
figure; plot(time_step_label,meanmag);title('Mean Traction N2')

%figure; plot(meanmagdisp);


