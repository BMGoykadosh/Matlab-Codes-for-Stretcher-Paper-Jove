 clc
clear
close all
% This is the RBS correction. You should do it for the time series as you
% wil need it. Also correct the tryp image and phase image as
% you will need them for:
%A) Tryp-comparison for the RBS 
%B) Phase-for checking the cell. If doing dag or calcium images this
%will also be useful

pixel_ratio = 0.32594; % um dimension/pixels ex: 665 um x 665 um, 2048 pix = 0.32594
%This depends on the magnification of the microscope-unit length over pixel
%size   


basepath='C:\Users\Caroline McCormick\Desktop\Benjie\Fall_2024\stretcher_experiments\11_22strretcher\'; %UPDATE for basepath. Telling the code to this respective folder 
processedpath2 = 'Uni_1'; %naming a new folder to store the new images within the basepath folder
seriespath = '\2024_11_22_09_23_25--old_uni_1\'; %Folder where you find all your series Images 
tryppath = '\2024_11_22_09_23_25--old_uni_1\';%%Folder where you find all your Tryp Images 

mkdir(fullfile(strcat(basepath,'\',processedpath2)));
%creates the path to save data. Makes the RGB Shifted folder
time_frames = 4; %UPDATE: for timeframe number (remember base zero vs standardcount) Number of images in the time frames
positions = 2; %how many postions you have, number of places within a gel you want to get data for. ie., ROIs  


xshift = 50; yshift = 50; xwidth = 1919; ywidth = xwidth;  
% Define maximum expected rigid body shift (xshift, yshift) and final image dimensions (xwidth, ywidth)
% xshift/yshift: maximum anticipated pixel displacement in x/y directions due to rigid body motion
% xwidth/ywidth: crop dimensions (1919×1919) to exclude boundary artifacts after shift correction
% Original images are 2048×2048; cropping removes edge pixels that may shift out of the imaging plane
% Boundary pixels contain unreliable data after alignment and are discarded
% N = (xwidth+1)/resolution;

% check tryp img name and series img name
 

     
for j=1:positions
    
    disp(j)
    %Here we set up the TRYP image we want to use for the comparison  
    bead_name_b = sprintf('pre_1mumbead_red.tif',j); %UPDATE: Tryp Image file name
    bead_trp_file_b = fullfile(basepath,tryppath,bead_name_b);
    bead_trp_img_b = (imread(bead_trp_file_b));
    
    bead_name = sprintf('step05_1mumbead_red.tif',j); %UPDATE: Series file name. Update the proper iterative loop in the "comment_rbs_correction_function_stretcher_test" function
    bead_trp_file = fullfile(basepath,seriespath,bead_name);
    bead_trp_img = (imread(bead_trp_file)); 
    %This is image 1 of the series-it is useful when we want to correct 
    % the image series-tryp to 1, 1 to rest
 
    rbs_correction_function_stretcher(basepath, seriespath, processedpath2, bead_trp_img, bead_trp_img_b, j, time_frames, xshift, yshift, xwidth, ywidth);

end