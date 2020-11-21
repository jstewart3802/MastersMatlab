[PCAloadings, PCAscore, PCAvariances] = pca(normalised_spectra);


%% Plots the %variance explained by each PC
maxPC = input('Enter how many PCs to display; ');
figure('Name','Percent variance explained by First 10 PCs');
plot(1:maxPC,100*cumsum(PCAvariances(1:maxPC))/sum(PCAvariances(1:maxPC)),'-ob');
xlabel('Number of PCs')
ylabel('Percentage of variance explained in x')
clear maxPC;