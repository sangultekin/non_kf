%Stochastic Search Kalman Filter
%Algorithm 1 from the paper:
%S. Gultekin and J. Paisley. Nonlinear Kalman filtering with divergence minimization, 
%IEEE Transactions on Signal Processing, vol. 65, no. 23, pp. 6319-6331.

%Author: San Gultekin, Columbia University, san.gultekin@gmail.com

function [OUT] = func_SKF(IN,MEAS)

x_filter = IN.x_filter;
P        = IN.P;
F        = IN.F;
Q_CV     = IN.Q_CV;

n   = size(x_filter,2);

x_bar = x_filter(:,1);
P_bar = P;
for t=2:n
    R = MEAS(t).R;
    z_true = MEAS(t).z_true;
    
    x_predic = F*x_bar;
    P_predic = F*P_bar*F' + Q_CV;  
    
    P_bar = P_predic; 
    x_bar = x_predic;
    for iter=1:20
        %reset gradients
        grd_x = zeros(size(x_bar));
        grd_P = zeros(size(P_bar));
        
        %use jacobian as approximation
        x_exp = x_bar;
        [H, z_exp] = aux_jac(x_predic, MEAS(t).loc);
        z_g = z_true - z_exp + H*x_exp;

        
        
        %stochastic search
        S = 500; %no samples to approximate integrals
        %note: with 20 iterations and 500 sampels the total size is 1e4,
        %which matches the other filters' sample size
        
        P_bar = .5*(P_bar+P_bar'); 
        samp_x = mvnrnd(x_bar,P_bar,S)';
        
        z_samp = aux_predict_meas(samp_x,MEAS(t).loc);
        
        tmp1 = repmat(z_true,1,S) - z_samp;
        f = sum(tmp1.*(inv(R)*tmp1),1);

        tmp2 = repmat(z_g,1,S) - H*samp_x;
        g = sum(tmp2.*(inv(R)*tmp2),1);
        
        tmp3 = inv(P_bar)*samp_x - repmat(P_bar\x_bar,1,S);
        grd_x = grd_x + (1/S) * sum( bsxfun(@times,f-g,tmp3) , 2);
        
        tmp4 = samp_x - repmat(x_bar,1,S);
        tmp_mat1 = bsxfun(@times, sign(f-g).*sqrt(abs(f-g)) , tmp4 );
        tmp_mat2 = bsxfun(@times, sqrt(abs(f-g)) , tmp4 );
        grd_P = grd_P + (1/S) * ( .5*inv(P_bar)*(tmp_mat1*tmp_mat2')*inv(P_bar) );
        grd_P = grd_P + (1/S) * sum(f-g) * (-.5*inv(P_bar));
        
        %these are the gradients of the objective
        grd_P = -.5*grd_P;
        grd_x = -.5*grd_x;
		
        %from the gradients calculate natural gradients
        %this is necessary for covariance stability
        inv_CP = P_bar;
        inv_Cx = P_bar; 
        
        grd_x = grd_x + P_predic\x_predic - P_predic\x_bar + ...
            (H'*inv(R)*z_g) - (H'*inv(R)*H)*x_bar; 
        grd_P = grd_P + .5*inv(P_bar) -.5*inv(P_predic) -.5*(H'*inv(R)*H);  
        
        %use standard O(1/t) step size
        rho = 1/iter;
        
        %posterior updates with covariance check for positive definiteness
        P_prop = P_bar + rho*inv_CP*grd_P*inv_CP;
        x_prop = x_bar + rho*inv_Cx*grd_x;
        while sum(eig(P_prop)<0) > 0
            rho = rho/2;           
            P_prop = P_bar + rho*inv_CP*grd_P*inv_CP;
            x_prop = x_bar + rho*inv_Cx*grd_x;
        end
        
        x_bar = x_prop;
        P_bar = P_prop;
    end
  
    %state estimation
    x_filter(:,t) = x_bar;
end

OUT = x_filter;
end