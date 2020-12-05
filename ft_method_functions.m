clearvars();
tic;
sub_id = 1;
files_sub1 = [3,4,15,16,18,21,26];
seizure_times = zeros(7,2);
seizure_times(1,:) = [2996,3036];
seizure_times(2,:) = [1467,1494];
seizure_times(3,:) = [1732,1772];
seizure_times(4,:) = [1015,1066];
seizure_times(5,:) = [1720,1810];
seizure_times(6,:) = [327,420];
seizure_times(7,:) = [1862,1963];
window_type = "hanning";
F_sampling = 256; %Hz
window_length = 1; %in s
window_shift = F_sampling*window_length/2;
gamma = 2;
filename = 'chb01_15.edf';
l2 = get_lambda(filename, window_type, window_shift, window_length, gamma);
plot(l2)
xlabel("Time (s)");
ylabel("\lambda_2");
title("\lambda_2 for Sub1 Session 15, Hanning Window 50% overlap");
toc;


% for k=1:7
%     if files_sub1(k)<10
%         filename = strcat('chb0',int2str(sub_id), '_0',int2str(files_sub1(k)),'.edf');
%     else
%         filename = strcat('chb0',int2str(sub_id), '_',int2str(files_sub1(k)),'.edf');
%     end
%     lambda2 = get_lambda(filename, window_type, window_shift, window_length, gamma);
%     lam2_thr = 1;
%     figure();
%     indices = find(lambda2>lam2_thr);
%     plot(lambda2)
%     
%     title(int2str(files_sub1(k)));
%     
% end










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%function space%%
%%%%%%%%%%%%%%%%%%
function lam2t = get_lambda(filename, window_type, window_shift, window_length, gamma)
    [head, data1] = get_data(filename);
    [n_channels, n_samples] = size(data1);
    F_sampling = head.frequency(1);
    assert(F_sampling==256);
    [data_final,data_norm, data_cen, baseline, peak] = normalize_data(data1); % peak of data_cen, baseline of data1
    N = F_sampling*window_length; % no. of nonzero samples corresponding to window length in s
    w = get_window(N, window_type);
    [stft_data, t_len, f_len] = get_stft(data_final, w, window_shift);
    stft_norm = stft_data;
    del_f = 1/window_length ; %frequency resolution
    freq_bpf = [2,10]*del_f;
    [filtered_stft, nonzero_indices] = bpf(stft_norm, del_f, freq_bpf); % del_f = 1Hz
    norm_const = 1000; % best for 1000
    P_t = total_power_from_stft(filtered_stft, norm_const);
    sigma_P = power_coupling_coeff(P_t);
    [sigma_phi, phase_stft] = phase_coupling_coeff(filtered_stft);
    [adj_mat, coupling_product] = adjacency_matrix(sigma_P, sigma_phi,gamma);
    degree_mat = degree_matrix(adj_mat);
    laplac_mat = laplacian_matrix(degree_mat, adj_mat);
    [lam2t, eig_matrix] = fiedler_eig(laplac_mat);
    P_t_all = sum(P_t, 1);
    P_t_all_norm = P_t_all/max(P_t_all(:));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [header,data] = get_data(filename)
    [header, data] = edfread(filename);

end
%%%%%%%%%%%%%%%%%
function [data_final, data_norm,data_centred, baseline, peak_val] = normalize_data(data_raw)
    
    %data_raw size = (n_channels, n_samples)
    [n_ch, n_samp] = size(data_raw);
    baseline = (1/n_samp)*sum(data_raw, 2);
    assert(isequal(size(baseline),[n_ch,1]));
    
    data_centred = data_raw-baseline;
    peak_val = max(abs(data_centred), [], 2);
    data_norm = data_centred./peak_val;
    reference_avg = sum(data_norm, 1)/n_ch;
    assert(isequal(size(reference_avg), [1, n_samp]));
    data_final = data_norm - reference_avg;
end
%%%%%%%%%%%%%%%%%%
function plot_timeseries(data, n_channels_to_plot, F_sampling, session_name)
    [~, n_samples] = size(data);
    hold on;
    t_axis = (1/F_sampling)*(1:n_samples);
    for ch=1:n_channels_to_plot
        plot(t_axis, (ch)+data(ch, :)); % plot shifted versions of data
    end
    title(strcat("Time series plot (ref to average) for channels 1-",int2str(n_channels_to_plot) , ' subject ', session_name));
    xlabel("Time (s)");
    ylabel("Channel number");
    yticks(1:n_channels_to_plot);

end
%%%%%%%%%%%%%%%%%%%
function w=get_window(window_length, type)

    if strcmp(type, 'rectangular')
        w = ones(1,window_length);
    elseif strcmp(type, 'hanning')
        w = 0.5*(1-cos(2*pi*(0:window_length-1)/window_length));
    elseif strcmp(type, 'kaiser')
        w = kaiser(window_length);
        w = w';
        assert (isequal(size(w), [1, window_length]));
    end
end
%%%%%%%%%%%%%%%%%%
function [stft, n_samp_stft, N]=get_stft(data,window, window_shift) % to undo normalization effect on transform
    N = length(window);
    assert(isequal(size(window), [1,N]));
    [n_chan, n_samp_data] = size(data);
%     window_shift = N; % window shift for adjacent segment: N->non overlapping
    n_samp_stft = 1+((n_samp_data-N)/window_shift);
    stft = zeros(n_chan, n_samp_stft, N); % 3rd dim = num of freq from FFT
    for frame=1:n_samp_stft
        data_segment = data(:, 1+ (window_shift*(frame-1)):(window_shift*(frame-1))+N);
        windowed_data = data_segment.*window;
        fft_segment = fft(windowed_data, N, 2); % fft along row
        stft(:,frame,:) = fft_segment;
    end  

end
%%%%%%%%%%%%%%%%%%%%
function stft_norm = normalize_stft(stft_data)
    % stft_data has size (n_channels, n_timepoints, n_frequencies)
    % n_timepoints is the number of points obtained by sliding the window
    % this function removes the zero peak in computed fft
    max_val = max(abs(stft_data(:)));
    stft_norm = stft_data/max_val;

end
%%%%%%%%%%%%%%%%%%%%
function [filtered_stft, nonzero_indices]= bpf(stft, del_f, freq_range)
    nonzero_indices = 1 + (freq_range(1)/del_f:freq_range(2)/del_f);
    filtered_stft = zeros(size(stft));
    filtered_stft(:,:,nonzero_indices) = stft(:,:,nonzero_indices);
    
end
%%%%%%%%%%%%%%%%%%%%
function Pow_spect_sum = total_power_from_stft(stft_filtered, normalization_const)
    mag_spec = abs(stft_filtered);
    pow_spec = mag_spec.^2;
    n_freq = sum(stft_filtered(1,1,:)~=0);
    Pow_spect_sum = sum(pow_spec, 3)/normalization_const; % sum over frequencies
    
end
%%%%%%%%%%%%%%%%%%%%%%
function sigma_P = power_coupling_coeff(pow_spect_sum)
    [n_channels, n_timepoints] = size(pow_spect_sum);
    sigma_P = zeros(n_channels, n_channels, n_timepoints);
    
    for n_i=1:n_channels
        for n_j=1:n_i-1
            temp_n = pow_spect_sum(n_i,:)+pow_spect_sum(n_j,:);
            sigma_P(n_i, n_j,:) = temp_n; 
            sigma_P(n_j,n_i,:) = temp_n;
        end
    end
                
end
%%%%%%%%%%%%%%%%%%%%%%%
function [sigma_phi, phase_stft] = phase_coupling_coeff(stft_filtered)
    [n_chan, n_timepoints, ~] = size(stft_filtered);
    n_freq = sum(abs(stft_filtered(1,1,:))~=0); % number of nonzero components (after filtering)
    assert(n_freq==9);
    phase_stft = angle(stft_filtered); % gives phase in radians
    sigma_phi = zeros(n_chan, n_chan, n_timepoints);
    for n_i=1:n_chan
        for n_j=1:n_i-1
            temp = (1/n_freq)*sum(sin(phase_stft(n_j, :, :)-phase_stft(n_i,:,:)),3);
            assert(isequal(temp, real(temp)),"sin giving imag");
            sigma_phi(n_i, n_j,:) = temp;
            sigma_phi(n_j,n_i,:) = -temp;
        end
    end
            
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [adj, prod_pow_phas]=adjacency_matrix(sigma_Pow, sigma_phas, gamma) % gamma= hyperparameter
    assert(isequal(size(sigma_Pow), size(sigma_phas)),"Error- Power and phase coupling coefficient matrices do not match in dimension");
    prod_pow_phas = sigma_Pow.*sigma_phas;
    adj = (1- exp(-(abs(prod_pow_phas).^gamma)));
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function deg_mat = degree_matrix(adj_mat)
%     [n_chan, ~, n_timepoints] = size(adj_mat);
%     degree_vals = sum( adj_mat~=0,2);
%     deg_mat = zeros(n_chan, n_chan, n_timepoints);
%     for t=1:n_timepoints
%         deg_mat(:,:,t) = diag(degree_vals(:,:,t));
%     end
%     
    [n_chan, ~, n_timepoints] = size(adj_mat);
    degree_vals = sum( adj_mat,2);
    deg_mat = zeros(n_chan, n_chan, n_timepoints);
    for t=1:n_timepoints
        deg_mat(:,:,t) = diag(degree_vals(:,:,t));
    end
    



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lap=laplacian_matrix(degree_mat, adjacency_mat)
    assert(isequal(size(degree_mat),size(adjacency_mat)));
    lap = degree_mat- adjacency_mat;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lam2_t, eigen_mat] = fiedler_eig(laplac_mat)
    [n_chan, ~, n_timepoints] = size(laplac_mat);
    lam2_t = zeros(1,n_timepoints);
    eigen_mat = zeros(n_chan, n_timepoints);
    for t=1:n_timepoints
        eigen_lap = eig(squeeze(laplac_mat(:,:,t)));
        eig_lap_sorted = sort(eigen_lap, 'ascend');
        eigen_mat(:,t) = eig_lap_sorted;
        lam2_t(1,t) = eig_lap_sorted(2); % second smallest eigenvalue  
    end
end

% code used on 1 subject
% one 1-hour epoch is being analysed
% filename = 'chb01_15.edf';
% [head, data1] = get_data(filename);
% [n_channels, n_samples] = size(data1);
% F_sampling = head.frequency(1);
% [data_final,data_norm, data_cen, baseline, peak] = normalize_data(data1); % peak of data_cen, baseline of data1
% n_chan_plot = n_channels;
% % plot_timeseries(data_final,n_chan_plot, F_sampling, ' chb01\_15'); 
% window_length = 1; % in s
% N = F_sampling*window_length; % no. of nonzero samples corresponding to window length in s
% w = get_window(N, "rectangular");
% % plot(w);
% shift = N;
% [stft_data, t_len, f_len] = get_stft(data_final, w, shift);
% % stft_norm = normalize_stft(stft_data);
% stft_norm = stft_data;
% del_f = 1/window_length ; %frequency resolution
% freq_bpf = [2,10]*del_f;
% [filtered_stft, nonzero_indices] = bpf(stft_norm, del_f, freq_bpf); % del_f = 1Hz
% % s = size(squeeze(stft_data(2,:,:)));
% % [x,y] = ndgrid(1:15,1:s(1));
% % figure(1);
% % surf(x,y,abs(squeeze(stft_data(2,:,1:15))).');
% % figure(2);
% % surf(x,y,abs(squeeze(filtered_stft(2,:,1:15))).');
% norm_const = 1000; % best for 1000
% P_t = total_power_from_stft(filtered_stft, norm_const);
% % ch_id = 2;
% % plot(P_t(ch_id,:));
% % title("Sum of squared power spectrum in [2,10] Hz");
% % xlabel("Time (s)");
% % ylabel("SUM(|H(w)|^2");
% % 
% sigma_P = power_coupling_coeff(P_t);
% chan = 1:n_chan_plot;
% % imagesc(chan, chan, sigma_P(:,:,1));
% % xlabel("channel number");
% % ylabel("channel number");
% % title("Power related coupling coefficient at t=1s");
% [sigma_phi, phase_stft] = phase_coupling_coeff(filtered_stft);
% % chan = 1:n_chan_plot;
% % imagesc(chan, chan, sigma_phi(:,:,1));
% % xlabel("channel number");
% % ylabel("channel number");
% % title("Phase related coupling coefficient at t=1s");
% gamma = 2;
% [adj_mat, coupling_product] = adjacency_matrix(sigma_P, sigma_phi,gamma);
% % imagesc(chan, chan, adj_mat(:,:,3));
% % xlabel("channel number");
% % ylabel("channel number");
% % title("Adjacency Matrix at t=1s");
% degree_mat = degree_matrix(adj_mat);
% % degree_mat = degree_mat/max(degree_mat(:));
% % imagesc(chan, chan, degree_mat(:,:,3));
% % xlabel("channel number");
% % ylabel("channel number");
% % title("Degree Matrix at t=1s");
% % degree_mat = eye(n_channels); % same at all time points
% laplac_mat = laplacian_matrix(degree_mat, adj_mat);
% [lambda2_t, eig_matrix] = fiedler_eig(laplac_mat);
% % plot(lambda2_t);
% % xlabel("Time (s)");
% % ylabel("\lambda_2");
% % title("\lambda_2 vs time for Subject chb01\_15");
% P_t_all = sum(P_t, 1);
% P_t_all_norm = P_t_all/max(P_t_all(:));
% % plot Total power
% % plot(P_t_all_norm, 'b');
% % xlabel("Time (s)");
% % ylabel("Normalized values");
% % hold on;
% % % plot seizure limits
% % % plot lambda_2
% % plot(lambda2_t/max(lambda2_t(:)),'r');
% % line([1732,1732],[0,1],'Color','black','LineStyle','--' );
% % line([1772,1772],[0,1],'Color','black','LineStyle','--' );
% % title("Comparison Between Total Power and \lambda_2");
% % legend("Sum of Power","\lambda_2","Seizure limits")