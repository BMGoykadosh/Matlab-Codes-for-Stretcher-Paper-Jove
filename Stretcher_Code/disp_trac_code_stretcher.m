clear
close all
% Displacement and traction calculation using block-matching method
% Compares bead positions between reference and deformed states to compute displacement field
% 
% Required inputs:
%   - Original image series path (undeformed reference state)
%   - Corrected trypsin image (stress-free reference for traction calculation)
%   - Corrected image series path (deformed states with rigid body shift correction applied)
%
% Outputs: displacement field (u) and traction stress field (t)

pixel_ratio = 0.32594;  %665.60/2048; % 40X 0.1613um/pixel,   10x  0.651855,   20x 0.32594
crop_folder_beads = 'C:\Users\Caroline McCormick\Desktop\Benjie\Fall_2024\stretcher_experiments\11_1_second_test_stretcher\RBS_Corrected_Uniaxial_stretch_test'; %UPDATE for rbs data. Corrected images 
 
f0 = 8; pois = 0.47; young = 13e3; %13e3; %0.3e3;pois = 0.47= poisson's ratio doesn't change if material is always NuSil
xshift = 50; yshift = 50; xwidth = 1920; ywidth = 1920; %Match these to the numbers in the RBS code, xwidth here is +1 to the xwidth in RBS correction code
process_folder = 'processed_64block_900_thresh_to_prev_image\'; %%Name and create folder to store your data. UPDATE 64 block and threshold size depending on parameters in sams beads file. This is just a naming normenclature and you can name it however you see fit
mkdir(fullfile(strcat(crop_folder_beads,'\',process_folder)));

for jj = 1:29   % time frames %UPDATE for image number. If you have 29 images, then this changes to j=1:29 (':' means to)
%Load the proper files to recover the displacments 
        
 
         cell_name = sprintf('Uniaxial_1_step_%02d_A.tif',jj);%,iiii,j); % RBS corrected series file %UPDATE for series image
         bead_cell_file = fullfile(crop_folder_beads,cell_name);
         bead_cell_img = double(imread(bead_cell_file));

         tryp_name = sprintf('Uniaxial_1_step_%02d_B.tif',jj);%,iiii,j); %RBS corrected tryp file%UPDATE
         bead_tryp_file = fullfile(crop_folder_beads,tryp_name);
         bead_tryp_img = double(imread(bead_tryp_file));

         [xv,yv,uxv,uyv] = comment_Sam_beads_imcorr_v2_stretcher_test(bead_tryp_img, bead_cell_img); % UPDATE block size, bead intensity in the function. Explanation in function
         %xv and yv are pos, uxv and uyv are disp, pix units
         xv = xv.*pixel_ratio;      yv = yv.*pixel_ratio; %check units now in microns
         uxv = uxv.*pixel_ratio;    uyv = uyv.*pixel_ratio;
         N = sqrt(length(xv));      box = size(bead_cell_img,1)/N; 
         spacing = box*pixel_ratio; 
         x = reshape(xv,N,N); 
         y = reshape(yv,N,N); 
         ux = reshape(uxv,N,N); %-mean(mean(ux) try to clean a bit more
         uy = reshape(uyv,N,N); 
         displacements = figure;
         surf(x,y,sqrt(ux.^2+uy.^2)*pixel_ratio); set(gcf,'Renderer','zbuffer');hold on; %surf is magnitude of displacement, qiver is direction vector %Height is mag of disp at that point
         view(2); shading interp; hh = colorbar; caxis([0 0.15]);%We took out caxis stuff but can use later 
         h = quiver3(x,y,100*ones(size(x)),ux,uy,zeros(size(ux)),1.5); %increase index 3 to raise direction arrows above "overlay", 1.5 is scale
         ylabel(hh, 'Microns', 'FontSize',12,'FontWeight','bold')
         set(h,'Color','w','LineWidth',1);
         title('Displacements','FontSize',18,'FontWeight','bold');axis ij;axis image;
         displacements.Position = [100 100 500 500];
         opfile = fullfile(strcat(crop_folder_beads,'\'),process_folder,sprintf('disp_map_t%02d.png',j));
         print(opfile,'-dpng','-r0')
    %this math calculates tractions based on variables and other things. 
    % u = K t, where K is the elasticity matrix. Since we have u and K, 
    % T = K^-1 u. We calculate this with FTTC method 
        % UNCONSTRAINED TRACTIONS
        a1 = (1.0 + pois) * (1.0 - pois) / (pi * young);
        b1 = (1.0 + pois) * pois / (pi * young);
        c1 = (1.0 + pois) * pois / (pi * young);

        clear sx sy;
        for ii = 1:(N/2)
            sx(ii,:) = 0:((N/2)-1);
            sy(ii,1:(N/2)) = (N/2)-ii;
        end

        kx = [ sx  sx-(N/2);  sx  sx-(N/2) ];
        ky = [ sy-(N/2)  sy-(N/2);  sy  sy ];
        ky = flipud(ky);
        kx(:,(N/2+1)) =  kx(:,(N/2+1));
        ky((N/2+1),:) =  ky((N/2+1),:);
        k_abs = sqrt(kx.^2 + ky.^2);

        alpha = atan2(ky,kx);
        if kx(1,1) == 0 && ky(1,1) == 0
            alpha(1,1) = pi/2;
        end

        Cx = ((k_abs * young) / (2 * (1 - pois^2))) .* (1 - pois + pois .* (cos(alpha)).^2);
        Cy = ((k_abs * young) / (2 * (1 - pois^2))) .* (1 - pois + pois .* (sin(alpha)).^2);
        D  = ((k_abs * young) / (2 * (1 - pois^2))) .* (pois .* sin(alpha) .* cos(alpha));

        D(:,(N/2+1)) = zeros(N,1);
        D((N/2+1),:) = zeros(1,N);

        % Calculate Unconstrained Tractions-fast fourier transform
        Dx = fft2(2*pi*ux / (N*spacing)); 
        Dy = fft2(2*pi*uy / (N*spacing));

        Tx = Cx.*Dx + D.*Dy; %TX and TY are real and imaginary components
        Ty = D.*Dx + Cy.*Dy;

        tx_uc = real(ifft2(Tx)); %this just takes the real parts and is what we use
        ty_uc = real(ifft2(Ty));


        trac = figure;
        surf(x,y,sqrt(tx_uc.^2+ty_uc.^2)); set(gcf,'Renderer','zbuffer');hold on;%get mag of total trac
        view(2); shading interp; hh = colorbar; %caxis([0 50]); 
        h = quiver3(x,y,100*ones(size(x)),tx_uc,ty_uc,zeros(size(tx_uc)),1.5);
        ylabel(hh, 'Pa', 'FontSize',12,'FontWeight','bold')
        set(h,'Color','w','LineWidth',1);
        title('Tractions','FontSize',18,'FontWeight','bold');axis ij;axis image;
        trac.Position = [100 100 500 500];
        opfile = fullfile(strcat(crop_folder_beads,'\'),process_folder,sprintf('trac_map_t%02d.png',j));
        print(opfile,'-dpng','-r0') 

        close all
        opfile = fullfile(crop_folder_beads,process_folder,sprintf('cell_forces_t%02d.mat',j));%.mat is variable storage in matlab %This file is created into the RBS folder
        %running this code-it is topology maps of forces
        save(opfile)

        count = count+1;
    %end

 end
 