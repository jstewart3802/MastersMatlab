function [amidenormMaster1 namespectra1] = file_read(num_files, subfolders1)
    ROIspectraMaster = {};  % this is the raw spectra from impoly
    amidenormMaster = {}; % this has been normalised and concatenated
    ROIspectraMasterNUM = {}; % this is the number of specrtra in the ROI


    for i = 1:num_files%%% CHANGE THIS TO 1:2 IF YOU ONLY WANT TO SAY 2 length(subfolders1)
        dirName = [char(subfolders1(i)) '/']; 
        imgfilename = dir([dirName '*.dms']); 
        img = [imgfilename.folder , '/',  imgfilename.name ]; 
        %img = [filename(i).folder , '/',  filename(i).name ];

        [wavenumbers, data, width, height, filename1, acqdate]=readvarianmosaic_v4_1(img); 
        [pathstr, name, ext] = fileparts(filename1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This is David's impoly bit which allows user to draw a polygon around
    % image

        figure('Name', name); 
        imagesc(sum(data,3));

        set(gcf, 'units','normalized','outerposition',[0 0 1 1]); % make it full screen
        ROI  = impoly(gca); %imfreehand imrect(gca);

        wait(ROI); 
        ROIBW = ROI.createMask; 
        [yCoordinates, xCoordinates] = find(ROIBW);
        saveas(gcf,sprintf('fullImage%d.png',i));
        close

        x_min = min(xCoordinates); x_max = max(xCoordinates); x_range = x_max-x_min+1;
        y_min = min(yCoordinates); y_max = max(yCoordinates); y_range = y_max-y_min+1;
        dataROI2 = -1000*ones(y_range, x_range, size(data,3));
        for k = 1:numel(yCoordinates)
            dataROI2(yCoordinates(k)-y_min+1, xCoordinates(k)-x_min+1, :) = data(yCoordinates(k), xCoordinates(k), :);
        end
        dataROI2_im = sum(dataROI2,3);
        min_val = min(dataROI2_im(dataROI2_im>=-1000));
        for k1 = 1:y_range
            for k2 = 1:x_range
                if dataROI2_im(k1,k2)<=-1000
                    dataROI2_im(k1,k2) = min_val;
                end
            end
        end

        ROIspectra = zeros(numel(yCoordinates),size(data,3));
        for k = 1:numel(yCoordinates)
            ROIspectra(k, :) = data(yCoordinates(k), xCoordinates(k), :);
        end

        ROIspectraMasterNUM1 = size(ROIspectra);
        ROIspectraMasterNUM2 = ROIspectraMasterNUM1(1); 
        ROIspectraMasterNUM{i} = ROIspectraMasterNUM2; 
        ROIspectraMaster{i} = ROIspectra; 
        ROIspectraMaster1 = cell2mat(ROIspectraMaster'); % create matrix of spectra

        figure('Name', name);
        imagesc(sum(dataROI2_im,3));
        colorbar
        saveas(gcf,sprintf('magnified%d.png',i));
        close

    %%%%%%%%% end of David's impoly %%%%%%%%%%%%%%%

        meanSpectra = mean(ROIspectra);    
        Amidenorm = meanSpectra;
        % this prints out the file you are working on
        i

        [g, gN] = grp2idx(imgfilename.name);
        lstNames{i} = gN;

        % get rid of empty in lenAmidenorm and lstnames
        % lenAmidenorm(lenAmidenorm==0) = []; 
        lstNames = lstNames(~cellfun('isempty',lstNames)); 
        lstNames';
        namespectra1 = vertcat(lstNames{:}); 
        char(namespectra1);

        amidenormMaster{i} = Amidenorm; 
        amidenormMaster1 = cell2mat(amidenormMaster'); % create matrix of spectra

        save('trainingData.mat', 'amidenormMaster1', 'ROIspectraMaster1', 'lstNames', 'ROIspectraMasterNUM'); 
    end 
    %% this separates into just WT or TG
    namespectra1 = vertcat(lstNames{:}); %Do this if just 1 master file
end