function results = LFP_metrics(data, varargin)

time = data.time;

data_fields = fieldnames(data);

v_fields = data_fields(contains(data_fields, '_V'));

layers = {'', 'sup', 'L4', 'deep'};

column_names = {'C1', 'C2'};

li = 1;

LFP = nan(length(time), length(column_names)*length(layers));
Spec = nan(length(time), 100, length(column_names)*length(layers));

for c = 1:length(column_names)
    
    for l = 1:length(layers)
        
        field_index = contains(v_fields, layers{l}) & contains(v_fields, column_names{c});
        
        if strcmp(layers{l}, 'deep'), field_index = field_index & ~contains(v_fields, 'dendrite'); end
       
        LFP_fields = v_fields(field_index);
        
        LFP_pops = cellfun(@(x) x{1}, regexp(LFP_fields, '_V', 'split'), 'unif', 0);
        
        LFP_vars = {};
        
        for p = 1:length(LFP_pops)
            
            LFP_vars = [LFP_vars; data_fields(contains(data_fields, [LFP_pops{p}, '_C']) & contains(data_fields, '_iSYN_sSYNpre'))];
            
        end
        
        thisLFP = zeros(length(time), 1);
        
        for v = 1:length(LFP_vars)
           
            thisLFP = thisLFP + detrend(nanmean(data.(LFP_vars{v}), 2)); 
            
        end
        
        LFP_labels{li} = [column_names{c}, layers{l}];
        
        LFP(:, li) = thisLFP;
        
        sampling_freq = 1000*10;
        
        [Pow(:, li), f] = pmtm(thisLFP, [], [], sampling_freq);
        
        Spec(:, :, li) = wavelet_spectrogram(thisLFP, sampling_freq, 1:100, linspace(3,14,100));
        
        li = li + 1;
        
    end
    
end

results = struct('LFP', LFP, 'Spec', Spec, 'Pow', Pow, 'f', f);
results.LFP_labels = LFP_labels;