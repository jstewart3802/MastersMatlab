[PCAloadings, PCAscore, PCAvariances] = pca(normalised_spectra);

PCi = input("Select first of 4 PC loadings to display; ");
PCj = input("Select second PC loading to display; ");
PCk = input("Select third PC loading to display; ");
PCl = input("Select last PC loading to display; ");
PCI = int2str(PCi);
PCJ = int2str(PCj);
PCK = int2str(PCk);
PCL = int2str(PCl);

%% plot loadings
figure('Name',['PCA loadings for PC',PCI,', PC',PCJ,', PC',PCK,' and PC',PCL]);
load1 = plot(x_crop,PCAloadings(:,PCi));
hold on
load2 = plot(x_crop,PCAloadings(:,PCj));
load3 = plot(x_crop,PCAloadings(:,PCk));
load4 = plot(x_crop,PCAloadings(:,PCl));
legend(['PC',PCI],['PC',PCJ],['PC',PCK],['PC',PCL]);
xlabel('Wavenumber (cm^{-1})');
ylabel('PCA Loading');
hold off
clear PCi PCI PCJ PCj PCk PCK PCl PCL;