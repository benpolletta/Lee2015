function sim_struct = Lee2015initSimStruct(varargin)

keys = {'solver', 'tspan', 'dsfact',...
    'save_data_flag', 'parallel_flag', 'num_cores', 'compile_flag', 'cluster_flag', 'qsub_mode', 'one_solve_file_flag',...
    'debug_flag', 'analysis_functions', 'analysis_options',...
    'column', 'ach_flag', 'vary',...
    };

values = {'euler', [0 3000], 100,...
    1, 0, 16, 0, 0, 'array', 1,...
    0,{},{},...
    'a1_2013', 0, {},...
    };

sim_struct = init_struct(keys, values, varargin{:});