%Auxiliary script for computing Jacobian for the sensor network problem

%Author: San Gultekin, Columbia University, san.gultekin@gmail.com

function [H, z_predic] = aux_jac(x,loc)
num = size(loc,2);

dif   = repmat(x([1 3],:),1,num) - loc;
dif_r = sqrt( sum(dif.^2,1) );

z_predic = dif_r';

H      = zeros(num,4);
H(:,1) = dif(1,:)' ./ dif_r';
H(:,3) = dif(2,:)' ./ dif_r';
end