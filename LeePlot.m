function LeePlot(data)

column_names = {'C1', 'C2'};

pops = {'supRS', 'supSI', 'supFS', 'L4RS', 'L4FS', 'deepRSaxon', 'deepIBaxon', 'deepFS', 'deepSI'};

pop_colors = [1 0 0; 0 .7 0; 0 0 1; 1 0 0; 0 0 1; 1 0 0; 0 0 0; 0 0 1; 0 .7 0];

no_pops = length(pops);

ha = tight_subplot(no_pops,2);

for p = 1:no_pops
   
    for c = 1:2
    
        axes(ha(2*(p - 1) + c))
    
        fig_handle = dsPlot(data, 'variable', [column_names{c}, pops{p}, '*V'],...
            'plot_type', 'raster', 'lock_gca', 1, 'suppress_textstring', 1);
        
        ax_handle = get(fig_handle, 'Children');
        
        line_handles = get(ax_handle, 'Children');
        
        for l = 1:length(line_handles{1})
            
            set(line_handles{1}(l), 'Color', pop_colors(p, :)) % 'Marker', '.', 
            
        end
        
        if mod(c, 2) == 0, ylabel(''), end
        
        if p < no_pops, set(gca, 'XTickLabel', ''), xlabel(''), end
        
    end
    
end
    
end