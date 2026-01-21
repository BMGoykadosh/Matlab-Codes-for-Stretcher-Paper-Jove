function rbs_correction_function_stretcher(basepath, seriespath, processedpath2, bead_trp_img, bead_trp_img_b, j, time_frames, xshift, yshift, xwidth, ywidth)

processedpath2folder = fullfile(basepath,processedpath2,'\');

    for jj =  1:time_frames
            cell_name_1 = sprintf('step%02d_1mumbead_red.tif',jj);%j,jj-1); %This calls series files iterativly. To iterate on matlab use: %02d
            bead_cell_file_1 = fullfile(basepath,seriespath,cell_name_1);
            bead_cell_img_1 = (imread(bead_cell_file_1));
%{
            phase_name_1 = sprintf('Position%03d--t%02d--C01.tif',j,jj-1); 
            phase_cell_file_1 = fullfile(basepath,seriespath,phase_name_1);
            phase_cell_img_1 = (imread(phase_cell_file_1));
%}
            shx = 0; shy = 0; 

            bead_cell_img_1_2 = bead_cell_img_1(yshift+shy:yshift+ywidth+shy, xshift+shx:xshift+xwidth+shx);
%Compares all images to the tryp image 
            bead_trp_img_2b = bead_trp_img_b(yshift+shy:yshift+ywidth+shy, xshift+shx:xshift+xwidth+shx);
            %if you ever want to know by how much the new image shifted,the
            %below line tells you          
            [shx,shy] = disp_on_blocks_v2(bead_cell_img_1_2,bead_trp_img_2b,size(bead_cell_img_1_2,1),0);%tryp image shift correction for RBS-where the work happens
            bead_cell_img_1_2 = bead_trp_img_b(yshift+shy:yshift+ywidth+shy, xshift+shx:xshift+xwidth+shx);%here we shift by the amount calculated in disp_on_blocks_v2
            bead_cell_img_2 = bead_cell_img_1(yshift+shy:yshift+ywidth+shy, xshift+shx:xshift+xwidth+shx);
            phase_cell_img_2 = phase_cell_img_1(yshift+shy:yshift+ywidth+shy, xshift+shx:xshift+xwidth+shx);

            imwrite(bead_cell_img_1_2,[processedpath2folder, sprintf('Uni_1_step_%02d_A.tif',jj)]);%j,jj for time series
            imwrite(phase_cell_img_2,[processedpath2folder, sprintf('phase_position%03d_run_%04dA.tif',j,jj)]);
            imwrite(bead_trp_img_2b,[processedpath2folder, sprintf('Uni_1_step_%02d_B.tif',jj)]);%j,jjsame tryp image for all timeframes

    end
end