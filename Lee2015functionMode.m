function [data, name, sim_spec] = Lee2015functionMode(sim_struct) % , simulation) %, varargin)

if nargin < 1; sim_struct = []; end
if isempty(sim_struct); sim_struct = Lee2015initSimStruct; end

% Setup path
restoredefaultpath

Today = datestr(datenum(date),'yy-mm-dd');

start_dir = pwd;

% Try-catch loops for changing directories
CD_ben;
CD_dave;
addpath(genpath(pwd)); % Need to add working directory (i.e., Lee2015) to path.

savepath = fullfile(pwd, 'Sims', Today);
mkdir(savepath);

unpack_sim_struct
    
[sim_spec, sim_label] = Lee2015simSpec(column, ach_flag, bottom_up_flag, top_down_flag, excluded, column_name);

if exist('simulation', 'var')
    
    sim_struct_in = sim_struct;
    
    [vary, sim_struct] = get_sim_vary(simulation);
    
    unpack_sim_struct
    
    sim_struct = sim_struct_in;
   
    sim_struct.vary = vary;
    
    [sim_spec, sim_label] = Lee2015simSpec(column, ach_flag, bottom_up_flag, top_down_flag, excluded, column_name);
    
    sim_label = [sim_label, '_', simulation];

end

Now = clock;
name = sprintf('%s_%g_%g_%.4g', sim_label, Now(4), Now(5), Now(6));

save(fullfile(savepath, [name, '_sim_spec.mat']), 'sim_spec', 'sim_struct', 'vary', 'name');

solver

if cluster_flag
    
    if ~isempty(vary)
        
        data = dsSimulate(sim_spec,'tspan',tspan,'downsample_factor',dsfact,'solver',solver,'coder',0,...
            'vary',vary,'verbose_flag',1,'cluster_flag',cluster_flag,...
            'debug_flag',debug_flag,'compile_flag',compile_flag,...
            'analysis_functions',analysis_functions,'analysis_options',analysis_options,...
            'overwrite_flag',1,'one_solve_file_flag',1,'qsub_mode',qsub_mode,...
            'save_data_flag',save_data_flag,'study_dir',fullfile(savepath, name));
        
        cd (start_dir)
        
        return
        
    else
        
        data = dsSimulate(sim_spec,'tspan',tspan,'downsample_factor',dsfact,'solver',solver,'coder',0,...
            'verbose_flag',1,'cluster_flag',cluster_flag,...
            'debug_flag',debug_flag,'compile_flag',compile_flag,...
            'analysis_functions',analysis_functions,'analysis_options',analysis_options,...
            'overwrite_flag',1,...
            'save_data_flag',save_data_flag,'study_dir',fullfile(savepath, name));
        
        cd (start_dir)
        
        return
        
        
    end
    
else
    
    if isempty(vary), parallel_flag = 0; end
    
    tic;
    
    data = dsSimulate(sim_spec,'tspan',tspan,'downsample_factor',dsfact,'solver',solver,'coder',0,...
        'vary',vary,'verbose_flag',1,'parallel_flag',parallel_flag,'num_cores',num_cores,...
        'debug_flag',debug_flag,'compile_flag',compile_flag,...
        'analysis_functions',analysis_functions,'analysis_options',analysis_options,...
        'save_data_flag',save_data_flag,'study_dir',fullfile(savepath, name));
    
    toc;
    
end

plotdata(data, column_name, savepath, name)

cd (start_dir)

end

function [sim_vary, sim_struct] = get_sim_vary(simulation)
    
if strcmp(simulation, 'LIP')
    
    sim_vary = {'(C2L4RS,C2L4FS)', 'rate', [50 100];...
        '(C1deepSI->C1deepIBdendrite,C2deepSI->C2deepIBdendrite)', 'gSYN', [0 0.5];...
        '(C1deepIBdendrite,C1deepRSdendrite,C2deepIBdendrite,C2deepRSdendrite)', 'Sfreq', 25;...
        '(C1deepIBaxon->C1deepIBdendrite,C2deepIBaxon->C2deepIBdendrite)', 'mechanism_list', '+iNMDA';...
        '(C1deepIBaxon->C1deepIBdendrite,C2deepIBaxon->C2deepIBdendrite)', 'gNMDA', [0 .05 .1];...
        '(C1deepIBdendrite,C1deepIBsoma,C1deepIBaxon,C2deepIBdendrite,C2deepIBsoma,C2deepIBaxon)',...
        'Iapp', permute([2 1 1 2 1 1; 5 4 4 5 4 4], [3 1 2]);...
        };
    
    sim_struct = struct();
    
elseif strcmp(simulation, 'LIP_theta')
    
    sim_vary = {'(C1L4RS,C1L4FS,C2L4RS,C2L4FS)', 'rate', [50 100];...
        '(C1deepIBdendrite,C1deepRSdendrite)', 'Soffset', 125;...
        '(C1deepIBdendrite,C1deepRSdendrite)', 'Strial', 250;...
        '(C2deepIBdendrite,C2deepRSdendrite)', 'Sonset', 125;...
        '(C2deepIBdendrite,C2deepRSdendrite)', 'Strial', 250;...
        '(C1deepSI->C1deepIBdendrite,C2deepSI->C2deepIBdendrite)', 'gSYN', 0;...
        '(C1deepIBdendrite,C1deepRSdendrite,C2deepIBdendrite,C2deepRSdendrite)', 'Sfreq', 25;...
        };
    
    sim_struct = struct();

elseif strcmp(simulation, 'Pul_Poisson')
    
    sim_vary = {'(C2L4RS,C2L4FS)', 'rate', [50 100];...
        '(C1deepIBdendrite,C1supRS,C1supFS)', 'frequencyPul', [0 15];...
        '(C1deepIBdendrite,C1supRS,C1supFS)', 'gExtPul', permute([[3 .2 .02]/2; 3 .2 .02], [3 1 2]);...
        '(C1deepSI->C1deepIBdendrite,C2deepSI->C2deepIBdendrite)', 'gSYN', 0;...
        '(C1deepIBdendrite,C1deepRSdendrite,C2deepIBdendrite,C2deepRSdendrite)', 'Sfreq', 25;...
        };
    
    sim_struct = struct();
    
elseif strcmp(simulation, 'Pul_Poisson_theta')
    
    sim_vary = {'(C1L4RS,C1L4FS,C2L4RS,C2L4FS)', 'rate', [50 100];...
        '(C1deepIBdendrite,C1supRS,C1supFS,C2deepIBdendrite,C2supRS,C2supFS)', 'frequencyPul', [0 15];...
        '(C1deepIBdendrite,C1supRS,C1supFS,C2deepIBdendrite,C2supRS,C2supFS)', 'gExtPul', permute([[3 .2 .02]/2; 3 .2 .02], [3 1 2]);...
        '(C1deepIBdendrite,C1supRS,C1supFS)', 'PulOff', 125;...
        '(C1deepIBdendrite,C1supRS,C1supFS)', 'PulTrial', 250;...
        '(C1deepIBdendrite,C1deepRSdendrite)', 'Soffset', 125;...
        '(C1deepIBdendrite,C1deepRSdendrite)', 'Strial', 250;...
        '(C2deepIBdendrite,C2supRS,C2supFS)', 'PulOn', 125;...
        '(C2deepIBdendrite,C2supRS,C2supFS)', 'PulTrial', 250;...
        '(C2deepIBdendrite,C2deepRSdendrite)', 'Sonset', 125;...
        '(C2deepIBdendrite,C2deepRSdendrite)', 'Strial', 250;...
        '(C1deepSI->C1deepIBdendrite,C2deepSI->C2deepIBdendrite)', 'gSYN', 0;...
        '(C1deepIBdendrite,C1deepRSdendrite,C2deepIBdendrite,C2deepRSdendrite)', 'Sfreq', 25;...
        };
    
    sim_struct = struct();
    
elseif strcmp(simulation, 'Pul_Spikes')
    
    sim_vary = {'(C2L4RS,C2L4FS)', 'rate', [50 100];...
        '(C1deepIBdendrite,C1supRS,C1supFS)', 'gSpikesPul', permute([[3 .2 .02]/2; 3 .2 .02], [3 1 2]);...
        '(C1deepIBdendrite,C1supRS,C1supFS,C2deepIBdendrite,C2supRS,C2supFS)', 'SfreqPul', 15;...
        '(C1deepSI->C1deepIBdendrite,C2deepSI->C2deepIBdendrite)', 'gSYN', 0;...
        '(C1deepIBdendrite,C1deepRSdendrite,C2deepIBdendrite,C2deepRSdendrite)', 'Sfreq', 25;...
        };
    
    sim_struct = struct();
    
elseif strcmp(simulation, 'Pul_Spikes_theta')
    
    sim_vary = {'(C1L4RS,C1L4FS,C2L4RS,C2L4FS)', 'rate', [50 100];...
        '(C1deepIBdendrite,C1supRS,C1supFS,C2deepIBdendrite,C2supRS,C2supFS)', 'gSpikesPul', permute([3 .2 .02; zeros(1,6)], [3 1 2]);...
        '(C1deepIBdendrite,C1supRS,C1supFS,C2deepIBdendrite,C2supRS,C2supFS)', 'SfreqPul', 15;...
        '(C1deepIBdendrite,C1supRS,C1supFS)', 'StrialPul', 250;...
        '(C1deepIBdendrite,C1supRS,C1supFS)', 'SoffsetPul', 125;...
        '(C1deepIBdendrite,C1deepRSdendrite)', 'Strial', 250;...
        '(C1deepIBdendrite,C1deepRSdendrite)', 'Soffset', 125;...
        '(C1deepIBdendrite,C1supRS,C1supFS)', 'StrialPul', 250;...
        '(C2deepIBdendrite,C2supRS,C2supFS)', 'SonsetPul', 125;...
        '(C2deepIBdendrite,C2deepRSdendrite)', 'Strial', 250;...
        '(C2deepIBdendrite,C2deepRSdendrite)', 'Sonset', 125;...
        '(C1deepSI->C1deepIBdendrite,C2deepSI->C2deepIBdendrite)', 'gSYN', 0;...
        '(C1deepIBdendrite,C1deepRSdendrite,C2deepIBdendrite,C2deepRSdendrite)', 'Sfreq', 25;...
        };
    
    sim_struct = struct();
    
elseif strcmp(simulation, 'L4alpha')
    
    sim_vary = {'L4RS->L4RS', 'mechanism_list', '+iNMDA';...
        'L4RS->L4RS', 'gNMDA', 0:.01:.05;...
        'L4RS->L4RS', 'gMg', 0;... 0:.2:1;...
        'L4RS', 'gKCNH', 0:1:10;...
        'L4RS', 'gExt', 0:.02:.2;...
        };

    keys = {'column', 'ach_flag', 'bottom_up_flag', 'top_down_flag',...
        'excluded', 'column_name',...
        };
    
    values = {'par_2015', 0, 1, 0,...
        {'deep', 'sup', 'FS'}, '',...
        };
    
    sim_struct = init_struct(keys, values);
    
end           
    
end

function plotdata(data, column_name, savepath, name)

if iscell(column_name) && length(column_name) > 1
    
    no_columns = length(column_name);
    
elseif isstr(column_name)
    
    column_name = {column_name};
    
    no_columns = 1;
    
end

no_pops = sum(contains(fields(data), '_V') & contains(fields(data), column_name{1}));

if no_columns == 1 & no_pops <= 4
        
    fig1 = figure;
    
    dsPlot(data, 'suppress_textstring', 1)
    
    saveas(fig1, fullfile(savepath, [name, '.fig']))
    
    save_as_pdf(fig1, fullfile(savepath, name))
    
    fig2 = figure;
    
    dsPlot(data, 'suppress_textstring', 1, 'plot_type', 'raster')
    
    saveas(fig2, fullfile(savepath, [name, '_raster.fig']))
    
    save_as_pdf(fig2, fullfile(savepath, [name, '_raster']))
        
else
    
    no_sims = length(data);
    
    for s = 1:no_sims
        
        fig1 = figure;
        
        for c = 1:no_columns
            
            axis = subplot(1, no_columns, c);
            
            if c == 1
                
                dsPlot(data(s), 'variable', [column_name{c}, '*V'], 'lock_gca', 1)
                
            else
                
                dsPlot(data(s), 'variable', [column_name{c}, '*V'], 'lock_gca', 1, 'suppress_textstring', 1)
                
            end
            
        end
        
        saveas(fig1, fullfile(savepath, [name, '_fig', num2str(s), '.fig']))
        
        save_as_pdf(fig1, fullfile(savepath, [name, '_fig', num2str(s)]))
        
        fig2 = figure;
        
        for c = 1:no_columns
            
            axis = subplot(1, no_columns, c);
            
            if c == 1
                
                dsPlot(data(s), 'variable', [column_name{c}, '*V'], 'lock_gca', 1, 'plot_type', 'rastergram')
                
            else
                
                dsPlot(data(s), 'variable', [column_name{c}, '*V'], 'lock_gca', 1, 'plot_type', 'rastergram', 'suppress_textstring', 1)
                
            end
            
        end
        
        saveas(fig2, fullfile(savepath, [name, '_fig', num2str(s), '_raster.fig']))
        
        save_as_pdf(fig2, fullfile(savepath, [name, '_fig', num2str(s), '_raster']))
        
    end
    
end

end


function CD_ben

    try
        cd /projectnb/crc-nak/brpp/Lee2015/
    catch err
        display(err)
    end

    try
        cd /Users/benjaminpittman-polletta/Documents/Science/Research_Projects/Lee2015/
    catch err
        display(err)
    end

end

function CD_dave

    try
        cd ~/src/Lee2015/
    catch err
        display(err)
    end

    % try
    %     cd /Users/benjaminpittman-polletta/Documents/Science/Research_Projects/Lee2015
    % catch err
    %     display(err)
    % end
    
end