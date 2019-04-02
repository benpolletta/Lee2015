function [data, name, sim_spec] = Lee2015functionMode(sim_struct) % , simulation) %, varargin)

if nargin < 1; sim_struct = []; end
if isempty(sim_struct); sim_struct = Lee2015initSimStruct; end

Today = datestr(datenum(date),'yy-mm-dd');

start_dir = pwd;

try
    cd /projectnb/crc-nak/brpp/Lee_2015
catch
end

try
    cd /Users/benjaminpittman-polletta/Documents/Science/Research_Projects/Lee_2015
catch
end

savepath = fullfile(pwd, 'Sims', Today);
mkdir(savepath);

Now = clock;
name = sprintf('Lee2015_%g_%g_%.4g', Now(4), Now(5), Now(6));

unpack_sim_struct

sim_spec = Lee2015simSpec(column, ach_flag, cluster_flag, excluded);

save(fullfile(savepath, [name, '_sim_spec.mat']), 'sim_spec', 'sim_struct', 'vary', 'name');

solver

if cluster_flag
    
    data = dsSimulate(sim_spec,'tspan',tspan,'downsample_factor',dsfact,'solver',solver,'coder',0,...
        'vary',vary,'verbose_flag',1,'cluster_flag',cluster_flag,...
        'debug_flag',debug_flag,'compile_flag',compile_flag,...
        'analysis_functions',analysis_functions,'analysis_options',analysis_options,...
        'overwrite_flag',1,'one_solve_file_flag',1,'qsub_mode',qsub_mode,...
        'save_data_flag',save_data_flag,'study_dir',fullfile(savepath, name));
    
    cd (start_dir)
    
    return

else
    
    tic;
    
    data = dsSimulate(sim_spec,'tspan',tspan,'downsample_factor',dsfact,'solver',solver,'coder',0,...
        'vary',vary,'verbose_flag',1,'parallel_flag',parallel_flag,'num_cores',num_cores,...
        'debug_flag',debug_flag,'compile_flag',compile_flag,...
        'analysis_functions',analysis_functions,'analysis_options',analysis_options,...
        'save_data_flag',save_data_flag,'study_dir',fullfile(savepath, name));
    
    toc;
        
end

dsPlot(data)
    
saveas(gcf, fullfile(savepath, [name, '.fig']))

save_as_pdf(gcf, fullfile(savepath, name))

dsPlot(data, 'plot_type', 'rastergram')
    
saveas(gcf, fullfile(savepath, [name, '_raster.fig']))

save_as_pdf(gcf, fullfile(savepath, [name, '_raster']))

cd (start_dir)

end