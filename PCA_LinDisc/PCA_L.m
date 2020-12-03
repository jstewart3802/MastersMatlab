%%plot loadings for two PCs
function PCA_L(PCi, PCj, PCI, PCJ, x, PCAloadings)
    figure('Name',['PCA loadings for PC',PCI,' and PC',PCJ]);
    load1 = plot(x,PCAloadings(:,PCi));
    hold on
    load2 = plot(x,PCAloadings(:,PCj));
    legend(['PC',PCI],['PC',PCJ]);
    xlabel('Wavenumber (cm^{-1})');
    ylabel('PCA Loading');
    hold off
end