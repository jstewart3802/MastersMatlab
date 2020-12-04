clear
clc  
close all
%% Put in your folder directory here - it should contain subfolders with all the files from images
[subfolders] = subdir('C:\Users\James\Documents\Uni\Masters\FTIR_data_for_students\Full_data'); 

pattern = ["tg" , "wt"]; 
logicalSub = contains(subfolders, pattern); 
subfolders1 = subfolders(logicalSub); 

%% Set to either do all files in folder or a certain number
num_files = length(subfolders1);  
%num_files = 15;

%% This calls the file and gets the ROI and names of the file 
%%AND saves the spectrum from that region 
[amidenormMaster1, namespectra1] = file_read(num_files, subfolders1);


%% Crop spectra to fingerprint region and then normalise with snv normalistaion
cropped_spectra = amidenormMaster1(1:end,1:410);
normalised_spectra=snv(cropped_spectra);


%% Sorts the normalised_spectra into WT or TG
wavenumbers = 410; %needs to be right number of wavenumbers for data

[group2wtTG,WT,TG] = sort_spectra(namespectra1, normalised_spectra, wavenumbers);


%% Create x axis for cropped spectra
x = 1:1506;
x = rescale(x,1000,3900);%input range of real wavenumbers here
x_crop = x(1:410);


%% plot the cropped spectra
plot_spectra(x_crop,WT,TG);