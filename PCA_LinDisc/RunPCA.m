[PCAloadings, PCAscore, PCAvariances] = pca(normalised_spectra);
%use this to have user input values for which PC to use
PCi = input("Enter PC for x axis; ");
PCj = input("Enter PC for y axis; ");
maxPC = input('Enter how many PCs to display for variance explained; ');
% uncomment these to have this instead of user input
%PCi = 1;    %PC for x axis
PCI = int2str(PCi);
%PCj = 3;    %PC for y axis
PCJ = int2str(PCj);
PCAWT = PCAscore(group2wtTG,1:end);
PCATG = PCAscore(group2wtTG==0,1:end);

PCA_S(PCi, PCj, PCI, PCJ, PCAWT, PCATG);
PCA_L(PCi, PCj, PCI, PCJ, x_crop, PCAloadings);
PCA_E(PCAvariances, maxPC);

clear PCi PCj PCI PCJ maxPC;