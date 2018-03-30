%Auxiliary script for predicting measurements for the sensor network problem

%Author: San Gultekin, Columbia University, san.gultekin@gmail.com

function [z_samp] = aux_predict_meas(samp_x,loc)

no_sns  = size(loc,2);
no_samp = size(samp_x,2);

z_samp = zeros(no_sns,no_samp);
for i=1:no_sns
    dif   = samp_x([1 3],:) - repmat(loc(:,i),1,no_samp);
    dif_r = sqrt( sum(dif.^2,1) );
    
    z_samp(i,:) = dif_r;
end