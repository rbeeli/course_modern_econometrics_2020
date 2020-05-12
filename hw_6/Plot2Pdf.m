function Plot2Pdf(path)
    % Saves the current plot as PDF to disk.
    % Preserves original size and aspect
    % ratio of the plot in the PDF.
    set(gcf, 'Units','inches')
    pos = get(gcf, 'Position');   
    set(gcf, 'PaperUnits','inches');
    set(gcf, 'PaperSize', [pos(3) pos(4)]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3) pos(4)]);
    print(gcf, '-dpdf', path);
end
