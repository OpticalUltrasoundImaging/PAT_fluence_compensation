clear;
close all;

% Load coordinates [cm]
load('x.mat'); load('z.mat');
% voxel size
dx = x(2) - x(1); dz = z(2) - z(1);
% Load convex cone and compensated linear unmixing fluence spectra
F = load('F.mat').F;
fc = load('fc.mat').fc;
% Load absorption spectra
meCS = load('meCS.mat').meCS;
meNS = load('meNS.mat').meNS;

% Load theoretical absorptin spectra
mua_sNi = load('mua_sNi.mat').mua_sNi;
% Processed experimental data
load('regionInfo.mat');

% ratio between NiSO4 and CuSO4
c = 14.28;
% p0 spectra for convex cone method
p0_sNi = zeros([size(mua_sNi,2) size(F)]);
for i = 1:size(mua_sNi,2)
    p0_sNi(i,:,:) = F .* repmat(mua_sNi(:,i),1,size(F,2));
    for j = 1:size(F,2)
        p0_sNi(i,:,j) = p0_sNi(i,:,j) / norm(squeeze(p0_sNi(i,:,j)),2);
    end
end

% sNi prediction
% convex cone
sNi_cc = zeros(numel(z), numel(x));
% Linear Unmixing
sNi_lu = zeros(numel(z), numel(x));
% corrected linear unmixing
sNi_clu = zeros(numel(z), numel(x));

for i = 1:size(regionInfo,2)
    % p0 spectrum
    p = regionInfo(i).p0_mean;
    % Linear unmixing
    lu = linear_unmixing(p, [meNS meCS/c]);
    % Convex cone
    cc = convex_cone(p,p0_sNi);
    % Corrected linear unmixing
    clu =linear_unmixing(p./fc, [meNS meCS/c]);
    for ic = 1:size(regionInfo(i).Coordinates,1)
        sNi_cc(regionInfo(i).Coordinates(ic,1),regionInfo(i).Coordinates(ic,2)) = cc;
        sNi_lu(regionInfo(i).Coordinates(ic,1),regionInfo(i).Coordinates(ic,2)) = lu;
        sNi_clu(regionInfo(i).Coordinates(ic,1),regionInfo(i).Coordinates(ic,2)) = clu;
    end
end

figure;
imagesc(x,z,sNi_lu);
clim([0 1]); colorbar;
set(gca,'fontsize',14);
xlabel('x [cm]', 'fontsize',18);
ylabel('z [cm]', 'fontsize',18);
title('Linear Unmixing');

figure;
imagesc(x,z,sNi_cc);
clim([0 1]); colorbar;
set(gca,'fontsize',14);
xlabel('x [cm]', 'fontsize',18);
ylabel('z [cm]', 'fontsize',18);
title('Convex Cone');

figure;
imagesc(x,z,sNi_clu);
clim([0 1]); colorbar;
set(gca,'fontsize',14);
xlabel('x [cm]', 'fontsize',18);
ylabel('z [cm]', 'fontsize',18);
title('Corrected Linear Unmixing');

