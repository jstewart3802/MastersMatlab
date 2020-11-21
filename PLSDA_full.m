clear
clc  
close all
%%
% filename = dir('/Users/jcgj201/Documents/MATLAB/FTIR data/allresults/*.dmt');
%% Unhash for mean spectra of an entire image image
% imgy = '/Users/jcgj201/Documents/MATLAB/FTIR data/17-10-2016/25-02-16-WT/25-02-16-wt.dms'; 
% [wavenumbers, data]=readvarianmosaic_v4_1(imgy); 
% dataSpectra = reshape(data, 344064, 1506); 
% mu = mean(dataSpectra); 
% figure, 
% plot(wavenumbers,mu)
% set(gca, 'XDir','reverse')


%% Put in your folder directory here - it should contain subfolders with all the files from images
% [subfolders] = subdir('/Volumes/Charlie_mac_backup/TAU Project');%subdir('/Users/jcgj201/Documents/MATLAB/FTIR data/Testing_moveResultsTau_data'); % get all the subdirectories
% [subfolders] = subdir('D:\Documents\Uni\Masters\FTIR_data_for_students\Tau_data'); 
[subfolders] = subdir('D:\Documents\Uni\Masters\FTIR_data_for_students\Full_data'); 
pattern = ["tg" , "wt"]; 
logicalSub = contains(subfolders, pattern); 
subfolders1 = subfolders(logicalSub); 
subfoldersT = subfolders1'; 

%% This calls the file and gets the ROI and names of the file 
% AND saves the spectrum from that region 
ROIspectraMaster = {};  % this is the raw spectra from impoly
amidenormMaster = {}; % this has been normalised and concatenated
ROIspectraMasterNUM = {}; % this is the number of specrtra in the ROI

for i = 1:length(subfolders1)%%% CHANGE THIS TO 1:2 IF YOU ONLY WANT TO SAY 2 
    dirName = [char(subfolders1(i)) '/']; 
    imgfilename = dir([dirName '*.dms']); 
    img = [imgfilename.folder , '/',  imgfilename.name ]; 
    %img = [filename(i).folder , '/',  filename(i).name ];
    
    [wavenumbers, data, width, height, filename1, acqdate]=readvarianmosaic_v4_1(img); 
    [pathstr, name, ext] = fileparts(filename1);
    
 %%%%%%%%%%% This is imrect which works great%%%%%%%%%%%%%%%%%
%     figure('Name',filename1,'NumberTitle','off');
%     imagesc(sum(data,3));
%     ROI  = imrect(gca); 
%     positionROI = wait(ROI);
%     close
%     positionROI = round(positionROI); 
%     %%%%%%%%%% take the ROI position
%     xCoordinates = positionROI(1):positionROI(1)+positionROI(3); % postionROI is a [4 number] 
%     yCoordinates = positionROI(2):positionROI(2)+positionROI(4); % and this is the way to put togther the square
%     

%     dataROI2 = squeeze(data(y,x,:));
% 
%     figure,
%     imagesc(sum(dataROI2,3));
%     colorbar
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
    
    % pick 2 parts of spectra
%     amide = meanSpectra(:, 1:400); 
%     lipid = meanSpectra(:, 984:1040); 
    
    %figure, plot(1:457, horzcat(amide, lipid))
    
    % Do SG 2nd diff
    
%     amideSG = deriv(amide,2,5,2); %deriv(x,der,window,order); 
%     lipidSG = deriv(lipid,2,5,2);
%     SGamidelipid = horzcat(amideSG, lipidSG);
    % smooth the spectra
%     amide = smooth(amide,5); 
%     lipid = smooth(lipid,5); 
%     
%     plot(1:457,[amide lipid])
    
    % baseline amide 
%       p = amide; 
%       n = 1:400; 
%       [z,a,it,ord,s,fct] = backcor(n,p, 2, 0.01, 'atq');
%       z = z'; 
%       amide_RB = z;
%       subtractedA = amide - amide_RB; 
      
      % baseline lipid 
%       p = lipid; 
%       n = 1:57; 
%       [z,a,it,ord,s,fct] = backcor(n,p, 2, 0.01, 'atq');
%       z = z'; 
%       lipid_RB = z;
%       subtractedB = lipid - lipid_RB; 

    % SNV normalise the AmideLipid spectra and plot
%     amideN = snv(amide); 
%     lipidN = snv(lipid); 

%     amideN = snv(subtractedA); 
%     lipidN = snv(subtractedB);
%     Amidenorm = horzcat(amideN, lipidN);
    
    Amidenorm = meanSpectra;
 
    % this prints out the file you are working on
    i
    
    [g, gN] = grp2idx(imgfilename.name);
    lstNames{i} = gN;
    
    % get rid of empty in lenAmidenorm and lstnames
%     lenAmidenorm(lenAmidenorm==0) = []; 
    lstNames = lstNames(~cellfun('isempty',lstNames)); 
    lstNames';
    namespectra1 = vertcat(lstNames{:}); 
    char(namespectra1);

    amidenormMaster{i} = Amidenorm; 
    amidenormMaster1 = cell2mat(amidenormMaster'); % create matrix of spectra
    
    save('trainingData.mat', 'amidenormMaster1', 'ROIspectraMaster1', 'lstNames', 'ROIspectraMasterNUM'); 
end   
 
%% ONly need this if have to concatenate 2 master files 
% amidenormMaster1 = vertcat(amidenormMaster2, amidenormMaster3); 
% lstNames = vertcat(lstNames2{:}, lstNames3{:});

%% plot all spectra 
% figure, plot(1:457, amidenormMaster1); 
%figure, plot(1:728, amidenormMaster1); 
% figure, plot(amidenormMaster{1, 1}); 
% legend(lstNames); 

%% dont need this i don't think
%charNamespectra1 = char(namespectra1); 

%% this separates into just WT or TG
% for i=1:length(namespectra1)
namespectra1 = vertcat(lstNames{:}); %Do this if just 1 master file

%% Make normalised spectra from amidenormmaster
cropped_spectra = amidenormMaster1(1:end,1:410);
normalised_spectra=snv(cropped_spectra);


%% test plot of normalised_spectra
%figure, plot(1:728, normalised_spectra);



%%
% namespectra1 = lstNames; 
group2wtTG = contains(namespectra1, 'wt' ); 
group2wtTGnorm = contains(namespectra1, 'wt' );
%% Sorts the amidenormMaster1 spectra into WT or TG
%groupIndx = repmat(group2wtTG, 1, length(normalised_spectra)); %logical index to the same size as the Spectra
%NumWT = numel(group2wtTG(group2wtTG == 1)); % gets how many WT 
%WTb = normalised_spectra(groupIndx); 
%WT = reshape(WTb, NumWT, 1506); % have to put in the correct number
%TGb = normalised_spectra(groupIndx == 0); 
%NumTG = numel(group2wtTG(group2wtTG == 0));
%TG = reshape(TGb, NumTG, 1506); % have to put in the correct number

%% Sorts the normalised_spectra into WT or TG
groupIndxnorm = repmat(group2wtTGnorm, 1, length(normalised_spectra)); %logical index to the same size as the Spectra
NumWTnorm = numel(group2wtTGnorm(group2wtTGnorm == 1)); % gets how many WT 
WTbnorm = normalised_spectra(groupIndxnorm); 
WTnorm = reshape(WTbnorm, NumWTnorm, 410); % have to put in the correct number
TGbnorm = normalised_spectra(groupIndxnorm == 0); 
NumTGnorm = numel(group2wtTGnorm(group2wtTGnorm == 0));
TGnorm = reshape(TGbnorm, NumTGnorm, 410); % have to put in the correct number


%% this plot concat spectrum
x = 1:1506;
x = rescale(x,1000,3900);%input range of real wavenumbers here
x_crop = x(1:410);


%figure,
%h1 = plot(x, WT, 'b'); 
%hold on
%h2 = plot(x, TG, 'r'); 
%h3 = plot(x, WTnorm, 'b--');
%h4 = plot(x, TGnorm, 'r--');
%set(gca,'XDir','reverse');
%legend([h1(1), h2(1), h3(1), h4(1)], 'WT', 'TG', 'WTnorm', 'TGnorm')
%xlabel('wavenumbers (needs converting to real units)')
%ylabel('absorbance (arbitary units)')
%hold off

%% plotnorm
%figure,
%h1 = plot(x, WT, 'b'); 
%hold on
%h2 = plot(x, TG, 'r'); 
%set(gca,'XDir','reverse');
%legend([h1(1), h2(1)], 'WTnorm', 'TGnorm')
%xlabel('wavenumbers (cm^-1)')
%ylabel('absorbance (arbitary units)')
%hold off

%% plotcropped
figure,
h1 = plot(x_crop, WTnorm, 'b'); 
hold on
h2 = plot(x_crop, TGnorm, 'r'); 
set(gca,'XDir','reverse');
legend([h1(1), h2(1)], 'WTnorm', 'TGnorm')
xlabel('wavenumbers (cm^-1)')
ylabel('absorbance (arbitary units)')
hold off


%% this puts it in standard format wavennumbers from 4500 to 1000 with middle bit taken out
% gapWT = repmat(NaN, size(WT, 1), 1049); 
% gapTG = repmat(NaN, size(TG, 1), 1049);
% WTgap = [WT(:, 1:400) gapWT WT(:, 401:457)];
% TGgap = [TG(:, 1:400) gapTG TG(:, 401:457)];
% 
% figure, 
% h1 = plot(wavenumbers, WTgap, 'b'); 
% hold on
% h2 = plot(wavenumbers, TGgap, 'r'); 
% legend([h1(1), h2(1)], 'WT', 'TG')
% set(gca, 'XDir','reverse')
% xlim([950,4500]); 
% ylim([-1,5]);
% xlabel('wavenumbers')
% ylabel('normalised (with standard variance) absorbances ')
%% perform PLS-DA on the 2 spectrum % MSE = Mean squares prediction error -  Can't do on 2 spectrum 
% [g1,gN1] = grp2idx(namespectra1); 
% Yaxis = g1; 
%% This is model for indivdal mouse brains
% [XL,yl,XS,YS,beta,PCTVAR, stats, mse] = plsregress(amidenormMaster1, group2wtTG, 10, 'CV',10);

%% This is model for just 2 - WT or TG
% folds = 1;  
% [XL,yl,XS,YS,beta,PCTVAR, stats, mse] = plsregress(amidenormMaster1, group2wtTG, folds, 'CV', folds);
% 
% %% Number of components
% figure
% plot(1:folds,cumsum(100*PCTVAR(2,:)),'-bo');
% xlabel('Number of PLS components');
% ylabel('Percent Variance Explained in y');
% 
% %% plot latent variable LV1 vs LV2 (FOR XScores - this is the SAMPLES) 
% figure
% % plot(XS(:,1),XS(:,2), 'o');
% % gscatter(XS(:,1),XS(:,2),g1); 
% gscatter(XS(:,1),XS(:,2),group2wtTG); 
% 
% X = XS(:,1); %LV1
% Y = XS(:,2); %LV2
% 
% xlabel('Latent Variable 1');
% ylabel('Latent Varaible 2');
% 
% %% plot latent variable LV3 vs wavelength (for Weights - this is the VARIABLES so 476 wavelengths) 
% figure, plot(1:457, mse.W(:, 1), 1:457, mse.W(:, 2), 1:457, mse.W(:, 3)) ;
% legend('show'); 
% ylim([-0.3, 0.3]); 
% 
% %% plot latent variable LV3 vs wavelength BUT PADDED OUT
% W = mse.W; 
% Wgap = repmat(NaN, 1049, 20);
% WwithGap = vertcat(W(1:400,:), Wgap, W(401:457,:));
% figure, 
% plot(wavenumbers, WwithGap(:, 1:3)); 
% set(gca, 'XDir','reverse')
% legend('LV1', 'LV2', 'LV3'); 
% ylim([-0.15, 0.15]);
% xlabel('wavenumbers (cm-1)')
% ylabel('Latent variables Variance')
% 
% %% Plot latent varibles with a break in the axis
% start = 1800; 
% stop = 3800; 
% figure, 
% h1=BreakXAxis(wavenumbers,WwithGap(:, 1),start,stop,100); 
% hold on
% h2=BreakXAxis(wavenumbers,WwithGap(:, 2),start,stop,100); 
% hold on
% h3=BreakXAxis(wavenumbers,WwithGap(:, 3),start,stop,100); 
% set(gca, 'XDir','reverse')
% ylim([-0.2, 0.2]);
% % xlim([1000, 3500]);
% legend('LV1', 'LV2', 'LV3'); 
% xlabel('wavenumbers (cm-1)')
% ylabel('Latent variables Variance')
% 
% %% test the number of components needed THIS IS A NEW ONE
% figure, plot(0:4,PLSmsep(2,:),'b-o')%0:10,PCRmsep,'r-^');
% xlabel('Number of components');
% ylabel('Estimated Mean Squared Prediction Error');
% legend('PLSR','location','NE');%legend({'PLSR' 'PCR'},'location','NE');
% 
% 
% %% try again 1==WT 0==TG this works 
% indices = crossvalind('Kfold',group2wtTG,67);
% cp = classperf(group2wtTG); % initializes the CP object
% for i = 1:67
%     test = (indices == i); train = ~test;
%     class = classify(X(test),Y(train),group2wtTG(train)); % X=LV1 and Y=LV2
%     % updates the CP object with the current classification results
%     classperf(cp,class,test)  
% end
% cp.CorrectRate % queries for the correct classification rate
% % numel(train)
% % numel(test)
% % numel(indices)
% 
% %%
% %%%%%%%%%%%%%% try with PCA %%%%%%%%%%%%%%%%%%%%%%%%%
% [PCALoadings,PCAScores,PCAVar] = pca(amidenormMaster1,'Economy',false);
% %% plot PCA 
% figure
% plot(100*cumsum(PCAVar(1:10))/sum(PCAVar(1:10)),'r-^');
% xlabel('Number of Principal Components');
% ylabel('Percent Variance Explained in data');
% legend('PCA','location','SE');
% %% Plot 1st 2 compenents against each other 
% %this need to be projected on onto LD but not sure how
% figure
% gscatter(PCAScores(:,10),PCAScores(:,15),group2wtTG); 
% % X = XS(:,1); %LV1
% % Y = XS(:,2); %LV2
% %%
% figure
% xlabel('PCA score 1 ');
% ylabel('PCA score 2');
% plot(1:457,PCALoadings(:,1:2),'-');
% xlabel('Variable');
% ylabel('PCA Loading');
% legend({'1st Component' '2nd Component'},'location','NW');
% %%
% P = PCALoadings; 
% Pgap = repmat(NaN, 1049, 457);
% PwithGap = vertcat(P(1:400,:), Pgap, P(401:457,:));
% figure, 
% plot(wavenumbers, PwithGap(:, 1:2)); 
% set(gca, 'XDir','reverse')
% legend('PCA 1', 'PCA 2'); 
% ylim([-0.25, 0.25]);
% xlabel('wavenumbers (cm-1)')
% ylabel('PCA loadings')
% 
% 
% %% testing other PLS-DA
% [Xl,Yl,Xs,Ys,beta,pctVar,mse,stats] = plsregress(amidenormMaster1,group2wtTG,20, 'CV', 20);
% figure, plot(1:457,stats.W(:, 1:3),'-');
% xlabel('Variable');
% ylabel('PLS Weight');
% legend({'1st Component' '2nd Component' '3rd Component'},  ...
% 	'location','NW');
% ylim([-0.25, 0.25])
% 


