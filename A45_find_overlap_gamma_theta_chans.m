
cd('/home/gattas/Desktop/e_loc')
% import text file data
%%
theta_specificty      = OFCFROTEMPCINGINSEC{:,1:3};
gamma_site_specificty =HCNCGammaDiscChans{:,1:3};
save('NC_theta_NC_gamma_discrim_elec_rois', 'theta_specificty', 'gamma_site_specificty')
%%

load('NC_theta_NC_gamma_discrim_elec_rois.mat')
dual_elecs = []
theta_elecs_to_delet = []
for iElec_gamma = 1:size(gamma_site_specificty,1)
    temp = (gamma_site_specificty(iElec_gamma,:));
    if ~isempty(find(ismember(theta_specificty, temp, 'rows')))
        dual_elecs = [dual_elecs iElec_gamma]
        theta_elecs_to_delet = [theta_elecs_to_delet find(ismember(theta_specificty, temp, 'rows'))];
    end
end
theta_idx = find(~ismember(1:size(theta_specificty,1),theta_elecs_to_delet))
theta_specificty = theta_specificty(theta_idx,:)

%%
Gamma_Theta_Eelecs = cat(1,gamma_site_specificty, theta_specificty)
for i = 1:size(Gamma_Theta_Eelecs,1)
column_1{i} = ['NC_Gamma_Theta' num2str(i)]
end
Dat_Table = table(column_1',Gamma_Theta_Eelecs(:,1),...
Gamma_Theta_Eelecs(:,2),Gamma_Theta_Eelecs(:,3));
%% make region color
colors = zeros(size(theta_specificty,1)+size(gamma_site_specificty,1),3);

% gamma in dark red
colors(1:size(gamma_site_specificty,1),1) = 255;

% theta in cyan
colors(size(gamma_site_specificty,1)+1:end,1) = 0;
colors(size(gamma_site_specificty,1)+1:end,2) = 255;
colors(size(gamma_site_specificty,1)+1:end,3) = 255;

% dual elecs in green
colors(dual_elecs,1) = 0;
colors(dual_elecs,2) = 255;
colors(dual_elecs,3) = 0;
