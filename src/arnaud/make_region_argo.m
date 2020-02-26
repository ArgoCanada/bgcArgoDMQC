clearvars

load /misc/1/home/laurenta/NEMO/workspace/matlab/mat/mask_nemo3oceans mask h lon lat

mask2 = mask;
mask2(:,1:83) = nan;
mask2(218:end,:) = nan;
mask2(209:end,189:end) = nan;
mask2(~isnan(mask2)) = 0;
mask2(isnan(mask2)) = 1;

[gridout,transect_xy] = transect_gridpoints_nemo_v2(mask2);