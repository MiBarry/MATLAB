function plot_hist_and_generalchi2pdf(dF_set, d, n_exp)
    dF_set = sort(dF_set);
    min_x = min(dF_set);
    max_x = max(dF_set);
    d_x = (max_x - min_x)/100;
    x_p = min_x:d_x:max_x;
    p = chi2pdf_general_bogdanov(x_p,d);
    figure
    h = histogram(dF_set,50);
    p = p * h.BinWidth * n_exp;
    hold on
    plot(x_p, p, 'r', 'LineWidth', 2)
    hold off

    
    n_i = h.Values / n_exp;
    p_i = zeros(1,h.NumBins);
%     h.NumBins
    for i=1:h.NumBins
        x_start = h.BinEdges(i);
        x_end = h.BinEdges(i+1);
        d_x_tmp = (x_end - x_start)/1000;
        x_tmp = x_start:d_x_tmp:x_end;
        y_tmp = chi2pdf_general_bogdanov(x_tmp,d);
        p_i(i) = trapz(x_tmp,y_tmp);
    end
    [res, p, stat] = chi2gof(1:length(p_i), 'Frequency', n_i*n_exp, 'Expected', p_i*n_exp, 'Alpha', 0.05, 'NParams', 1)

end
