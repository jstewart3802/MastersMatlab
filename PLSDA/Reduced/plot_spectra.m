function plot_spectra(x,WT,TG)
    figure,
    h1 = plot(x, WT, 'b'); 
    hold on
    h2 = plot(x, TG, 'r'); 
    %set(gca,'XDir','reverse');
    legend([h1(1), h2(1)], 'WT', 'TG')
    xlabel('Wavenumbers (cm^{-1})')
    ylabel('Absorbance (arbitary units)')
    hold off
end