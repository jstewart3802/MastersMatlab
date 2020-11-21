[PCAloadings, PCAscore, PCAvariances] = pca(normalised_spectra);

%use this to have user input values for which PC to use
PCi = input("Enter PC for x axis; ");
PCj = input("Enter PC for y axis; ");

% uncomment these to have this instead of user input
%PCi = 1;    %PC for x axis
PCI = int2str(PCi);
%PCj = 3;    %PC for y axis
PCJ = int2str(PCj);

%% plot loadings
figure('Name',['PCA loadings for PC',PCI,' and PC',PCJ]);
load1 = plot(x_crop,PCAloadings(:,PCi));
hold on
load2 = plot(x_crop,PCAloadings(:,PCj));
legend([load1, load2], ['PC ',PCI] ,['PC ', PCJ]);
xlabel('Wavenumber (cm^{-1})');
ylabel('PCA Loading');

%% Plot scores for WT or TG for PCs as set below
PCAWT = PCAscore(group2wtTG,1:end);
PCATG = PCAscore(group2wtTG==0,1:end);
figure('Name',['PCA scores for PC',PCI,' and PC',PCJ]);
sc1 = scatter(PCAWT(:,PCi),PCAWT(:,PCj),'b');
hold on
sc2 = scatter(PCATG(:,PCi),PCATG(:,PCj),'r');
legend([sc1(1), sc2(1)], 'WT', 'TG');
xlabel(['PC ',PCI]);
ylabel(['PC ',PCJ]);

%%Clear variables to avoid clutter
clear PCI PCJ PCi PCj;