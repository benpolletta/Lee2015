function psps = multi_Poisson_matlab(no_cells, inputs_per_cell, rate, tau_1, tau_d, T, dt)

t = 0:dt:T;

% EPSP for spikes at time t = 0.
psp = exp(-max(t - tau_1, 0)/tau_d); % tau_i*(exp(-max(t - tau_1, 0)/tau_d) - exp(-max(t - tau_1, 0)/tau_r)); %/(tau_d - tau_r);
psp = psp(psp > eps);    %?
psp = [zeros(1,length(psp)) psp]; %?

no_inputs = inputs_per_cell*no_cells;

if isscalar(rate)

    spike_arrivals = poissrnd((rate*dt/1000)*inputs_per_cell, no_cells, length(t));

else
    
    spike_arrivals = poissrnd(repmat((rate*dt/1000)*inputs_per_cell, 1, no_cells)');
    
end
    
psps = nan(size(spike_arrivals)); % Calculating EPSP experienced by each cell.

for c = 1:no_cells

    psps(c, :) = conv(spike_arrivals(c, :), psp, 'same');

end