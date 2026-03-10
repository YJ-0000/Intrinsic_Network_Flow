clear; clc;

%% Load data

target_dim = 27;
g_ica_result_files = dir(['results/disc_INF_G_lev_',num2str(target_dim,'%03d'),'_MIGP_results_*.mat']);
results_discovery = load(fullfile(g_ica_result_files(end).folder,g_ica_result_files(end).name));

g_ica_result_files = dir(['results/repl_INF_G_lev_',num2str(target_dim,'%03d'),'_MIGP_results_*.mat']);
results_replication = load(fullfile(g_ica_result_files(end).folder,g_ica_result_files(end).name));

is_only_cortex = false;

if is_only_cortex
    target_idx = 1:59412;
else
    target_idx = 1:size(results_discovery.Phi_orig_DR,1);
end

Phi_orig_disc = results_discovery.Phi_orig_DR(target_idx,:);
Phi_orig_repl = results_replication.Phi_orig_DR(target_idx,:);
lambda_disc = results_discovery.lambda;
lambda_repl = results_replication.lambda;

Phi_orig_disc(:,abs(imag(lambda_disc))<1e-6) = [];
Phi_orig_repl(:,abs(imag(lambda_repl))<1e-6) = [];

g_ica_result_files = dir(['results/INF_G_lev_REST1_ALL_',num2str(target_dim,'%03d'),'_MIGP_results_*.mat']);
results_group_rest1 = load(fullfile(g_ica_result_files(end).folder,g_ica_result_files(end).name));

g_ica_result_files = dir(['results/INF_G_lev_REST2_ALL_',num2str(target_dim,'%03d'),'_MIGP_results_*.mat']);
results_group_rest2 = load(fullfile(g_ica_result_files(end).folder,g_ica_result_files(end).name));

Phi_orig_group_rest1 = results_group_rest1.Phi_orig_DL(target_idx,:);
Phi_orig_group_rest2 = results_group_rest2.Phi_orig_DL(target_idx,:);
lambda_group_rest1 = results_group_rest1.lambda;
lambda_group_rest2 = results_group_rest2.lambda;

Phi_orig_group_rest1(:,abs(imag(lambda_group_rest1))<1e-6) = [];
Phi_orig_group_rest2(:,abs(imag(lambda_group_rest2))<1e-6) = [];

%% Discovery set vs. Replication
disp(':=== Mode Reproducibility - Discovery vs. Replication ===:');

sac_disc_repl = SAC(Phi_orig_disc, Phi_orig_repl);

fprintf('SAC value (Disc vs. Repl) = %.6f \n',sac_disc_repl);

max_DMs = size(Phi_orig_disc,2);

mac_array_disc_repl = nan(max_DMs/2,max_DMs/2);
for i = 1:max_DMs/2
    for j = 1:max_DMs/2
        mac = MAC(Phi_orig_disc(:,2*i-1),Phi_orig_repl(:,2*j-1));
        mac_array_disc_repl(i,j) = mac;
    end
end
figure; heatmap(mac_array_disc_repl);
xlabel('discovery');
ylabel('replicate');
title('Modal Assurance Criterion (disc vs. repl)');

sac_5_6 = SAC(Phi_orig_disc(:,2*[5,6]-1), Phi_orig_repl(:,2*[5,6]-1));

fprintf('SAC value with 5th, 6th modes (Disc vs. Repl) = %.6f \n',sac_5_6);

%% REST1 vs. REST2

disp(':=== Mode Reproducibility - REST1 vs. REST2 ===:');

sac_rest1_rest2 = SpaceMAC(Phi_orig_group_rest1(:,1:2:end), Phi_orig_group_rest2(:,1:2:end));

fprintf('SpaceMAC value (REST1 vs. REST2) = %.6f \n',sac_rest1_rest2);

max_DMs = size(Phi_orig_disc,2);

mac_array_rest1_rest2 = nan(max_DMs/2,max_DMs/2);
for i = 1:max_DMs/2
    idx_i = 2*i-1;
    for j = 1:max_DMs/2
        if j < max_DMs/2
            idx_j = 2*j-1;
        else
            idx_j = 2*j;
        end
        mac = MAC(Phi_orig_group_rest1(:,idx_i),Phi_orig_group_rest2(:,idx_j));
        mac_array_rest1_rest2(i,j) = mac;
    end
end
fig = figure; heatmap(mac_array_rest1_rest2);
xlabel('REST1');
ylabel('REST2');
title('Modal Assurance Criterion (REST1 vs. REST2)');
print(fig, 'temp_images/MAC_REST1_REST2.jpg', '-djpeg', '-r300');

s2mac_9 = S2MAC(Phi_orig_group_rest2(:,2*9-1),...
        Phi_orig_group_rest1(:,2*9-1),Phi_orig_group_rest1(:,2*10-1));
s2mac_10 = S2MAC(Phi_orig_group_rest2(:,2*10-1),...
        Phi_orig_group_rest1(:,2*9-1),Phi_orig_group_rest1(:,2*10-1));

fprintf('S2MAC value with REST2 - 9th vs. REST1 - 9th, 10th modes (REST1 vs. REST2) = %.6f \n',s2mac_9);
fprintf('S2MAC value with REST2 - 10th vs. REST1 - 9th, 10th modes (REST1 vs. REST2) = %.6f \n',s2mac_10);

s2mac_11 = S2MAC(Phi_orig_group_rest2(:,2*11-1),...
        Phi_orig_group_rest1(:,2*11-1),Phi_orig_group_rest1(:,2*12-1));
s2mac_12 = S2MAC(Phi_orig_group_rest2(:,2*12-1),...
        Phi_orig_group_rest1(:,2*11-1),Phi_orig_group_rest1(:,2*12-1));

fprintf('S2MAC value with REST2 - 11th vs. REST1 - 11th, 12th modes (REST1 vs. REST2) = %.6f \n',s2mac_11);
fprintf('S2MAC value with REST2 - 12th vs. REST1 - 11th, 12th modes (REST1 vs. REST2) = %.6f \n',s2mac_12);
