%code used on 1 subject
%one 1-hour epoch is being analysed
filename = 'chb01_15.edf';
[head, data1] = get_data(filename);
[n_channels, n_samples] = size(data1);
F_sampling = head.frequency(1);
[data_final,data_norm, data_cen, baseline, peak] = normalize_data(data1); % peak of data_cen, baseline of data1
n_chan_plot = n_channels;
% plot_timeseries(data_final,n_chan_plot, F_sampling, ' chb01\_15'); 
window_length = 1; % in s
N = F_sampling*window_length; % no. of nonzero samples corresponding to window length in s
w = get_window(N, "rectangular");
% plot(w);
shift = N;
[stft_data, t_len, f_len] = get_stft(data_final, w, shift);
% stft_norm = normalize_stft(stft_data);
stft_norm = stft_data;
del_f = 1/window_length ; %frequency resolution
freq_bpf = [2,10]*del_f;
[filtered_stft, nonzero_indices] = bpf(stft_norm, del_f, freq_bpf); % del_f = 1Hz
% s = size(squeeze(stft_data(2,:,:)));
% [x,y] = ndgrid(1:15,1:s(1));
% figure(1);
% surf(x,y,abs(squeeze(stft_data(2,:,1:15))).');
% figure(2);
% surf(x,y,abs(squeeze(filtered_stft(2,:,1:15))).');
norm_const = 1000; % best for 1000
P_t = total_power_from_stft(filtered_stft, norm_const);
% ch_id = 2;
% plot(P_t(ch_id,:));
% title("Sum of squared power spectrum in [2,10] Hz");
% xlabel("Time (s)");
% ylabel("SUM(|H(w)|^2");
% 
sigma_P = power_coupling_coeff(P_t);
chan = 1:n_chan_plot;
% imagesc(chan, chan, sigma_P(:,:,1));
% xlabel("channel number");
% ylabel("channel number");
% title("Power related coupling coefficient at t=1s");
[sigma_phi, phase_stft] = phase_coupling_coeff(filtered_stft);
% chan = 1:n_chan_plot;
% imagesc(chan, chan, sigma_phi(:,:,1));
% xlabel("channel number");
% ylabel("channel number");
% title("Phase related coupling coefficient at t=1s");
gamma = 2;
[adj_mat, coupling_product] = adjacency_matrix(sigma_P, sigma_phi,gamma);
% imagesc(chan, chan, adj_mat(:,:,3));
% xlabel("channel number");
% ylabel("channel number");
% title("Adjacency Matrix at t=1s");
degree_mat = degree_matrix(adj_mat);
% degree_mat = degree_mat/max(degree_mat(:));
% imagesc(chan, chan, degree_mat(:,:,3));
% xlabel("channel number");
% ylabel("channel number");
% title("Degree Matrix at t=1s");
% degree_mat = eye(n_channels); % same at all time points
laplac_mat = laplacian_matrix(degree_mat, adj_mat);
[lambda2_t, eig_matrix] = fiedler_eig(laplac_mat);
% plot(lambda2_t);
% xlabel("Time (s)");
% ylabel("\lambda_2");
% title("\lambda_2 vs time for Subject chb01\_15");
P_t_all = sum(P_t, 1);
P_t_all_norm = P_t_all/max(P_t_all(:));
% plot Total power
% plot(P_t_all_norm, 'b');
% xlabel("Time (s)");
% ylabel("Normalized values");
% hold on;
% % plot seizure limits
% % plot lambda_2
% plot(lambda2_t/max(lambda2_t(:)),'r');
% line([1732,1732],[0,1],'Color','black','LineStyle','--' );
% line([1772,1772],[0,1],'Color','black','LineStyle','--' );
% title("Comparison Between Total Power and \lambda_2");
% legend("Sum of Power","\lambda_2","Seizure limits")