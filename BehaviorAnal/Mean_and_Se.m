function [mean_x, se_x, qL_x, qU_x] = Mean_and_Se(x, flg, pLevel)
% write by rxj @ 07/26/2016
% calculate the mean and standard error (it also called standard error of the mean, S.E.M.)
% whereas, the mean squared error is the square of S.E.M.
if nargin == 1
    flg = 1; % flg = 1, mean of lengthway; flg = 2, mean of crosswise;
end
mean_x = nanmean(x, flg);
% pos = find(isnan(x) == 1);
% x(pos) = [];
se_x = nanstd(x, 0, flg)./sqrt(size(x, flg));
%se_x = std(x, 0, "omitnan")./sqrt(size(x, flg));
      
if nargin == 3
    qL_x = quantile(x, pLevel/2, flg);
    qU_x = quantile(x, 1-pLevel/2, flg);
end
