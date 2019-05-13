function LeePlot(data)

column_names = {'C1', 'C2'};

pops = {'supRS', 'supSI', 'supFS', 'L4RS', 'L4FS', 'deepRSaxon', 'deepIBaxon', 'deepFS', 'deepSI'};

no_pops = length(pops);

ha = tight_subplot(no_pops,2);

for p = 1:no_pops
   
    for c = 1:2
    
        axes(ha(2*(p - 1) + c))
    
        dsPlot(data, 'variable', [column_names{c}, pops{p}, '*V'], 'plot_type', 'raster', 'lock_gca', 1)
        
    end
    
end
    
end