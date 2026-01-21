function rbs_correction_function_stretcher(basepath, processedpath2,  bead_trp_img,  time_frames, xshift, yshift, xwidth, ywidth)
processedpath2folder = fullfile(basepath,processedpath2,'\');
    for jj =  1:time_frames-1
            cell_name_1 = sprintf('Image%03d.tif',jj+1); % 'Series001--t%2d'
            bead_cell_file_1 = fullfile(basepath,cell_name_1);
            bead_cell_img_1 = (imread(bead_cell_file_1));
            shx = 0; shy = 0; 
            bead_cell_img_1_2 = bead_cell_img_1(yshift+shy:yshift+ywidth+shy, xshift+shx:xshift+xwidth+shx);  
            bead_trp_img_2b = bead_trp_img(yshift+shy:yshift+ywidth+shy, xshift+shx:xshift+xwidth+shx);
            [shx,shy] = disp_on_blocks_v2(bead_trp_img_2b,bead_cell_img_1_2,size(bead_cell_img_1_2,1),0);
            disp(jj)
            disp(shx)
            disp(shy)
            bead_cell_img_1_2 = bead_cell_img_1(yshift+shy:yshift+ywidth+shy, xshift+shx:xshift+xwidth+shx);

            imwrite(bead_cell_img_1_2,[processedpath2folder, sprintf('run_%04dA.tif',jj)]);
            imwrite(bead_trp_img_2b,[processedpath2folder, sprintf('run_%04dB.tif',jj)]);
    end