function [sim_spec, label] = Lee2015simSpec(column, ach_flag, bottom_up_flag, top_down_flag, excluded, column_name)
% INPUTS:
% column (string or cell of strings): one of 'par_2015', 'a1_2015', 'a1_2013'.
% ach_flag (Boolean or vector of Booleans): determines whether conductances & connectivity 
%   reflect high or low cholinergic tone.
% excluded (1D cell of strings): populations to be left out of simulation.
% column_name (string or cell of strings): prefix(es) to be prepended to
%   population names, in case of multiple columns.

if nargin < 1, column = ''; end
if isempty(column), column = 'par_2015'; end
if nargin < 2, ach_flag = []; end
if isempty(ach_flag), ach_flag = 0; end
if nargin < 3, bottom_up_flag = []; end
if isempty(bottom_up_flag), bottom_up_flag = 0; end
if nargin < 3, top_down_flag = []; end
if isempty(top_down_flag), top_down_flag = 0; end
if nargin < 4, excluded = {}; end
if nargin < 5, column_name = ''; end

label = '';
if iscell(column)
    for c = 1:length(column), label = [label, column{c}]; end
else
    label = column;
end
label = [label, '_ach'];
for a = 1:length(ach_flag), label = [label, num2str(ach_flag(a), '%d')]; end
label = [label, '_td'];
for a = 1:length(top_down_flag), label = [label, num2str(top_down_flag(a), '%d')]; end
label = [label, '_bu'];
for a = 1:length(bottom_up_flag), label = [label, num2str(bottom_up_flag(a), '%d')]; end

multicolumn_flag = iscell(column_name) & (length(column_name) > 1);

%% In case more than one column.

if multicolumn_flag
    
    sim_spec = struct();
   
    for c = 1:length(column_name)
        
        if iscell(column)
            
            column_index = min(length(column), c);
            column_type = column{column_index};
            
        elseif isstr(column)
            
            column_type = column;
            
        end
       
        ach_index = min(length(ach_flag), c);
        top_down_index = min(length(top_down_flag), c);
        bottom_up_index = min(length(bottom_up_flag), c);
       
        [column_sim_spec, ex_label] = Lee2015simSpec(column_type, ach_flag(ach_index),...
            bottom_up_flag(bottom_up_index), top_down_flag(top_down_index),...
            excluded, column_name{c});
        
        label = [label, ex_label];
        
        if isfield(sim_spec, 'populations')
        
            sim_spec.populations = [sim_spec.populations column_sim_spec.populations];
            sim_spec.connections = [sim_spec.connections column_sim_spec.connections];
            
        else
            
            sim_spec = column_sim_spec;
            
        end
        
    end
    
    [pop_list{1:length(sim_spec.populations)}] = sim_spec.populations.name;
    
    if any(contains(column, '2013'))
       
        subcategories = {'deepIBaxon', 'deepRSaxon', 'supFS', 'supSI', column_name{:}};
        
        for s = 1:length(subcategories)
            
            eval([subcategories{s}, '_index = contains(pop_list, "', subcategories{s}, '");'])
            
        end
        
        for c = 1:length(column_name)
            
            eval(['deepE_index = (deepIBaxon_index | deepRSaxon_index) & ', column_name{c}, '_index;'])
            
            for d = [1:(c - 1) (c + 1):length(column_name)]
                
                eval(['column_supFS_index = supFS_index & ', column_name{d}, '_index;'])
                eval(['column_supSI_index = supSI_index & ', column_name{d}, '_index;'])
                
                deepEtosupFS = double(deepE_index)'*double(column_supFS_index);
                deepEtosupSI = double(deepE_index)'*double(column_supSI_index);
                
                fanout = 2*(deepEtosupFS + deepEtosupSI);
                gSYN = .2*deepEtosupFS + .3*deepEtosupSI;
                tauRx = .25*deepEtosupFS + 2.5*deepEtosupSI;
                tauDx = deepEtosupFS + 50*deepEtosupSI;
                
                [I, J] = find(gSYN > 0);
                
                C_index = length(sim_spec.connections);
                
                for con = 1:length(I)
                    
                    C_index = C_index + 1;
                    
                    sim_spec.connections(C_index).direction = [pop_list{I(con)}, '->', pop_list{J(con)}];
                    
                    sim_spec.connections(C_index).mechanism_list = {'iSYN'};
                    
                    sim_spec.connections(C_index).parameters = {'gSYN', gSYN(I(con), J(con)),...
                        'tauDx', tauDx(I(con), J(con)), 'tauRx', tauRx(I(con), J(con)),...
                        'ESYN', 0, 'fanout', fanout(I(con), J(con))};
                    
                end
            
            end
            
        end
        
    end
    
    return
    
%% In case single column.
    
else
    
    sim_spec = struct();
    
    pop_list = {'supRS', 'supFS', 'supSI',...
        'L4RS', 'L4FS',...
        'deepIBdendrite', 'deepIBsoma', 'deepIBaxon',...
        'deepRSdendrite', 'deepRSsoma', 'deepRSaxon',...
        'deepFS', 'deepSI'};
    
    no_pops = length(pop_list);
    
    multicomp_pops = {'deepIB', 'deepRS'};
    
    compartments = {'dendrite', 'soma', 'axon'};
    
    included = true(no_pops, 1);
    mc_pops_included = true(length(multicomp_pops), 1);
    comps_included = true(length(compartments), 1);
    
    for e = 1:length(excluded)
        
        included(contains(pop_list, excluded{e})) = false;
        
        mc_pops_included(contains(multicomp_pops, excluded{e})) = false;
        
        comps_included(contains(compartments, excluded{e})) = false;
        
    end
    
    pop_list(~included) = [];
    multicomp_pops(~mc_pops_included) = [];
    compartments(~comps_included) = [];
    
    if multicolumn_flag
        
        label = '';
        
    else
        
        label = sprintf('%s_ach%d_td%d_bu%d', column, ach_flag, top_down_flag, bottom_up_flag);
        
    end
    
    if ~isempty(excluded)
        
        excluded_label = ['_NO_', strcat(excluded{:})];
        
        included_label = ['_', strcat(pop_list{:})];
        
        if length(excluded_label) < length(included_label)
            
            label = [label, excluded_label];
            
        else
            
            label = [label, included_label];
            
        end
        
    end
    
    no_pops = length(pop_list); % pop_names = cellfun(@(x) [column_name, x], pop_list, 'UniformOutput', 0);
    
    param_list = {'gLeak', 'gM', 'gCaH', 'gCaL', 'gAR', 'Iapp', 'IappSTD', 'gExt', 'rate', 'gSpike', 'Sfreq'};
    
    no_params = length(param_list);
    
    conductance = get_conductance(column, ach_flag, bottom_up_flag, top_down_flag);
    
    conductance = conductance(:, included);
    
    for pop = 1:no_pops
        
        sim_spec.populations(pop).name = [column_name, pop_list{pop}];
        sim_spec.populations(pop).equations = pop_list{pop};
        sim_spec.populations(pop).size = 20;
        
        parameters = cell(1, 2*no_params);
        
        for param = 1:no_params
            
            parameters{2*param - 1} = param_list{param};
            parameters{2*param} = conductance(param, pop);
            
        end
        
        sim_spec.populations(pop).parameters = parameters;
        
    end
    
    mechanisms = {'iSYN', 'iNMDA', 'iGAP'};
    
    [fanout, gSYN, GJ, gNMDA, no_mechanisms, ESYN, tauRx, tauDx] = get_connectivity(column, ach_flag, pop_list, included);
    
    C_index = 0;
    
    for p = 1:no_pops
        
        for q = 1:no_pops
            
            if gSYN(p, q) > 0
                
                C_index = C_index + 1;
                
                sim_spec.connections(C_index).direction = [column_name, pop_list{p}, '->', column_name, pop_list{q}];
                
                sim_spec.connections(C_index).mechanism_list = mechanisms(1:no_mechanisms(p,q));
                
                sim_spec.connections(C_index).parameters = {'gSYN', gSYN(p, q),...
                    'tauDx', tauDx(p, q), 'tauRx', tauRx(p, q),...
                    'ESYN', ESYN(p), 'fanout', fanout(p, q),...
                    'fanoutNMDA', 10, 'gNMDA', gNMDA(p,q)}; % supEtosupI(p,q)*gNMDA(q),...
                
                if GJ(p, q) > 0
                    
                    sim_spec.connections(C_index).mechanism_list = [sim_spec.connections(C_index).mechanism_list {'iGAP'}];
                    
                    sim_spec.connections(C_index).parameters = [sim_spec.connections(C_index).parameters {'fanoutGAP', 7, 'gGAP', GJ(p,q)}];
                    
                end
                
            elseif GJ(p, q) > 0
                
                C_index = C_index + 1;
                
                sim_spec.connections(C_index).direction = [column_name, pop_list{p}, '->', column_name, pop_list{q}];
                
                sim_spec.connections(C_index).mechanism_list = {'iGAP'};
                
                sim_spec.connections(C_index).parameters = {'fanoutGAP', 7, 'gGAP', GJ(p,q)};
                
            end
            
        end
        
    end
    
    coupling = [0, 0.4, 0;... from dendrite
        0.2, 0, 0.3;... from soma
        0, 0.3, 0]; % from axon
    
    [I,J] = find(coupling > 0);
    
    for p = 1:length(multicomp_pops)
        
        pop_name = multicomp_pops{p};
        
        for c = 1:length(I)
            
            pre_name = [column_name, pop_name, compartments{I(c)}];
            
            post_name = [column_name, pop_name, compartments{J(c)}];
            
            if sum(contains(excluded, pre_name)) + sum(contains(excluded, post_name)) < 1
                
                C_index = C_index + 1;
                
                sim_spec.connections(C_index).direction = [pre_name, '->', post_name];
                
                sim_spec.connections(C_index).mechanism_list = {'iCOM'};
                
                sim_spec.connections(C_index).parameters = {'gCOM', coupling(I(c),J(c))};
                
            end
            
        end
        
    end
    
end

end

%% Getting conductances.

function conductance = get_conductance(column, ach_flag, bottom_up_flag, top_down_flag)

switch column
    
    case 'a1_2013' %, 'L4alpha'}
        
        conductance = [0.1*ones(1, 13);... % g_L
            0.5, 0, 8, 0.3, 0, 4, 0, 2, 4, 0, 2, 0, 4;... % g_M
            zeros(1, 5), 4, 0, 0, 1.6, zeros(1, 4);... % g_CaH
            zeros(1, 13);... % g_CaL
            zeros(1, 13);... % g_h
            0, 0, -1, 1, -2.5, 3, 2, 2, 3, 2, 2, -4, -1;... % Iapp 0, 0, -1, -1, 2, 2, 1, 1, 2, 1, 1, 0, -1;... 
            0.5*ones(1,3), 0, .5, .3, .1, .1, .3, .1, .1, .5, .8;... zeros(1,13);... % IappSTD
            .2, .02, zeros(1, 11);... % gExt % 0, 1, .03, 3, 0, 0, 3, zeros(1,4);...
            50, 50, zeros(1, 11);... % rate % 0, 100, 100, 250, 0, 0, 250, zeros(1,4);... 
            zeros(1, 13);... % gSpike
            zeros(1, 13)];... % Sfreq
            
        if bottom_up_flag
           
            conductance(end - 3, [4 5]) = [1 .03];
            conductance(end - 2, [4 5]) = 100;
            
        end
        
        if top_down_flag
           
            conductance(end - 1, [6 9]) = 3;
            conductance(end, [6 9]) = 20;
            
        end
        
    case 'a1_2015'
        
        conductance = [0.1*ones(1, 13);... % g_L
            0.3, 0, 6, 0.9, 0, 0.6, 0, 0.3, 0.6, 0, 0.3, 0, 3;... % g_M
            zeros(1, 5), 4, 0, 0, 2, zeros(1, 4);... % g_CaH
            0, 0, 1, zeros(1, 9), 0.3;... % g_CaL
            zeros(1, 13);... % g_h
            -2, 0, 0, -10, -5, -5, -1.4, 0.8, -5, -1.4, 0.8, 4, -4;... % Iapp
            0.5*ones(1,3), 0, .5, .3, .1, .1, .3, .1, .1, .5, .8;... % IappSTD
            .5, 0, 0, 0.1, .03, 0.1, 0, 0, 0.1, zeros(1, 4);... % gExt
            100*ones(1, 13); % rate
            zeros(1, 13);... % gSpike
            zeros(1, 13)];... % Sfreq
        
        if ach_flag
            
            conductance(2, :) = 2.*conductance(2, :);
            conductance(6, 3) = 6.6;
            conductance(6, end - 1) = -1.3;
            conductance(6, end) = 2.6;
            
        end
        
    case {'par_2015', 'L4alpha'}
        
        conductance = [0.1*ones(1, 4), 0.2, 0.1*ones(1, 8);... % g_L
            0.6, 0, 3, 1, 0, 2, 0, 8, 3, 0, 8, 0, 3;... % g_M
            zeros(1, 5), 4, 0, 0, 2, zeros(1, 4);... % g_CaH
            0, 0, 1, zeros(1, 9), 0.5;... % g_CaL
            zeros(1, 5), 0.1, 0, 0, 0.1, zeros(1,4);... % g_h
            -4, 2, 2, 0, -4, -5, -5, -5, -1, -3, -2, 0, -3;... % Iapp
            0.5*ones(1,3), 0, .5, .3, .1, .1, .3, .1, .1, .5, .8;... % IappSTD
            0.6, 0, 0, 0.2, 0.1, zeros(1, 8);... % gExt
            100*ones(1,13); % rate
            zeros(1, 13);... % gSpike
            zeros(1, 13)];... % Sfreq
        
        
        if ach_flag
            
            conductance(2, 1) = 0.3;
            conductance(2, 3) = 6;
            conductance(2, end) = 6;
            conductance(6, end) = 3;
            
        end
        
end

end

%% Getting connectivity matrices.

function [fanout, gSYN, GJ, gNMDA, no_mechanisms, ESYN, tauRx, tauDx] = get_connectivity(column, ach_flag, pop_list, included)

switch column
    
    case 'a1_2013'
        
        fanout = [5, 10, 10, 0, 0, 20, 20, 0, 0;... % from supRS
            5, 8, 5, 0, 0, 0, 0, 0, 0;... % from supFS
            5, 5, 0, 0, 0, 0, 0, 0, 0;... % from supSI
            5, 0, 0, 10, 10, 10, 10, 20, 0;... % from L4RS
            5, 0, 0, 10, 10, 0, 0, 0, 0;... % from L4FS
            0, 2, 2, 0, 0, 10, 10, 10, 10;... % from deepIB
            0, 2, 2, 0, 0, 10, 10, 10, 10;... % from deepRS
            0, 0, 0, 0, 0, 20, 20, 20, 10;... % from deepFS
            0, 0, 0, 0, 10, 20, 10, 10, 20]; % from deepSI
        
        gSYN = [.22, .3, .03, 0, 0, .212, .212, 0, 0;... % from supRS
            .4, .6, .1, 0, 0, 0, 0, 0, 0;... % from supFS
            .1, .2, 0, 0, 0, 0, 0, 0, 0;... % from supSI
            .2, 0, 0, .4, .2, .212, .212, .3, 0;... % from L4RS
            .02, 0, 0, 1, .3, 0, 0, 0, 0;... % from L4FS
            0, .2, .2, 0, 0, .02, .02, .12, .12;... % from deepIB
            0, .2, .2, 0, 0, .02, .02, .05, .15;... % from deepRS
            0, 0, 0, 0, 0, .1, .1, .5, .3;... % from deepFS
            0, 0, 0, 0, .4, .3, .3, .6, .4]; % from deepSI
        
%         if ach_flag
%             
%             gSYN(end, :) = 2.*gSYN(end, :);
%             
%         end
        
    case 'a1_2015'
        
        fanout = [10, 20, 10, 0, 0, 20, 20, 20, 0;... % from supRS
            10, 10, 10, 8, 0, 0, 0, 0, 0;... % from supFS
            10, 10, 2, 0, 0, 10, 10, 0, 0;... % from supSI
            10, 0, 0, 10, 20, 10, 10, 0, 0;... % from L4RS
            10, 0, 0, 20, 20, 0, 0, 0, 0;... % from L4FS
            0, 2, 2, 0, 0, 5, 5, 10, 10;... % from deepIB
            0, 2, 2, 0, 0, 5, 5, 10, 10;... % from deepRS
            0, 0, 0, 0, 0, 20, 20, 20, 10;... % from deepFS
            0, 10, 0, 0, 10, 20, 20, 10, 20]; % from deepSI
        
        gSYN = [.02, .3, .03, 0, 0, .05, .05, .4, 0;... % from supRS
            1.2, .2, .6, .03, 0, 0, 0, 0, 0;... % from supFS
            .02, 1, .2, 0, 0, .01, .01, 0, 0;... % from supSI
            .3, 0, 0, .03, .25, .01, .01, 0, 0;... % from L4RS
            .1, 0, 0, .3, .3, 0, 0, 0, 0;... % from L4FS
            0, .01, .25, 0, 0, .02, .01, .3, .2;... % from deepIB
            0, .01, .25, 0, 0, .01, .01, .3, .2;... % from deepRS
            0, 0, 0, 0, 0, .4, .4, .8, .8;... % from deepFS
            0, .02, 0, 0, .6, .5, .5, .03, .4]; % from deepSI
        
        if ach_flag
            
            gSYN(5, 4) = 2.*gSYN(5, 4);
            
        end
        
    case {'par_2015', 'L4alpha'}
        
        fanout = [10, 20, 10, 0, 0, 20, 20, 20, 0;... % from supRS
            10, 10, 10, 8, 0, 0, 0, 0, 0;... % from supFS
            10, 10, 0, 0, 0, 10, 10, 0, 0;... % from supSI
            10, 0, 0, 10, 20, 10, 10, 20, 0;... % from L4RS
            10, 0, 0, 20, 20, 0, 0, 0, 0;... % from L4FS
            0, 2, 2, 0, 0, 5, 5, 10, 10;... % from deepIB
            0, 2, 2, 0, 0, 0, 5, 10, 10;... % from deepRS
            0, 0, 0, 0, 0, 20, 20, 20, 10;... % from deepFS
            0, 10, 10, 0, 10, 0, 20, 10, 20]; % from deepSI
        
        gSYN = [.3, .3, .02, 0, 0, .01, .01, .01, 0;... % from supRS
            1.4, .5, .4, .2, 0, 0, 0, 0, 0;... % from supFS
            .05, .2, 0, 0, 0, .1, .1, 0, 0;... % from supSI
            .05, 0, 0, .02, .1, .01, .01, .02, 0;... % from L4RS
            .05, 0, 0, .5, .2, 0, 0, 0, 0;... % from L4FS
            0, .01, .1, 0, 0, .25, .02, .2, .2;... % from deepIB
            0, .01, .1, 0, 0, 0, .02, .2, .2;... % from deepRS
            0, 0, 0, 0, 0, .01, .1, .6, .05;... % from deepFS
            0, .5, .1, 0, .34, 0, .4, .4, .6]; % from deepSI
        
        if ach_flag
            
            gSYN(5, 4) = 2.*gSYN(5, 4);
            
        end
        
end
    
params = {fanout, gSYN};

for param = 1:length(params)
    
    param_mat = params{param};
    
    zero_pad = zeros(1, size(param_mat, 1));
    
    param_mat = [param_mat(1:5, :); zero_pad; zero_pad; param_mat(6, :);...
        zero_pad; zero_pad; param_mat(7:end, :)];
    
    zero_pad = zeros(size(param_mat, 1), 1);
    
    param_mat = [param_mat(:, 1:6), zero_pad, zero_pad, param_mat(:, 7), zero_pad, zero_pad, param_mat(:, 8:end)];
    
    params{param} = param_mat;
    
end

fanout = params{1}; gSYN = params{2};

fanout = fanout(included, included);
gSYN = gSYN(included, included);

subcategories = {'FS', 'SI', 'sup', 'deep', 'IBaxon', 'RSaxon'};

for s = 1:length(subcategories)
    
    eval([subcategories{s}, '_index = contains(pop_list, "', subcategories{s}, '");'])

end

E_index = ~(FS_index | SI_index);

deepEtosupSI = double(deep_index & E_index)'*contains(pop_list, 'supSI');

supEtosupI = double(sup_index & E_index)'*double(sup_index & (FS_index | SI_index));

no_pops = length(pop_list);

GJ = zeros(no_pops, no_pops);

switch column
        
    case 'a1_2013'
        
        GJ = .002*double(IBaxon_index)'*IBaxon_index + .002*double(RSaxon_index)'*RSaxon_index;
        gNMDA = FS_index*.04 + SI_index*.03;
        
        gNMDA = supEtosupI*diag(gNMDA);
        
    case 'a1_2015'
        
        gNMDA = (FS_index + SI_index)*.01;
        
        gNMDA = supEtosupI*diag(gNMDA);
    
    case 'par_2015'
        
        GJ = .04*double(IBaxon_index)'*IBaxon_index;
        gNMDA = FS_index*.01 + SI_index*.05;
        
        gNMDA = supEtosupI*diag(gNMDA);
        
    case 'L4alpha'
        
        L4RS_index = E_index & ~deep_index & ~sup_index;
        gNMDA = .05*double(L4RS_index)'*double(L4RS_index);
        
end

no_mechanisms = ones(no_pops, no_pops) + double(gNMDA > 0); %  + double(GJ > 0);

ESYN = (FS_index + SI_index)*(-80);

tauRx = repmat((E_index*.25 + (FS_index + SI_index)*.5)', 1, no_pops);
tauRx(logical(deepEtosupSI)) = 2.5;

tauDx = repmat((E_index + FS_index*8 + SI_index*20)', 1, no_pops);
tauDx(logical(deepEtosupSI)) = 50;

end