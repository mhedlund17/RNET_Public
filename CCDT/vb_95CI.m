function [CI95] = vb_95CI(var)
% Calculate 95% CI for any variable distribution

sem = std(var)/sqrt(length(var));
CI95 = tinv([0.025 0.975], length(var)-1);                    
CI95 = bsxfun(@times, sem, CI95(:)); 


% example figure
% errorbar(1,mean(x),ci95(1),ci95(2))
% hold on
% bar(1,mean(x))