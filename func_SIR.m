%Particle Filter

%Author: San Gultekin, Columbia University, san.gultekin@gmail.com



function [OUT] = func_SIR(IN,MEAS)

x_filter = IN.x_filter;
P        = IN.P;
F        = IN.F;
Q_CV     = IN.Q_CV;

n = size(x_filter,2);

%number of particles and initial locations
ns = 1e4;
PART = mvnrnd(x_filter(:,1),P,ns)';

for t=2:n
    R = MEAS(t).R;
    data_measurement = MEAS(t).z_true;
    dim = size(x_filter,1);
    
    %propagate wrt prior
    PART = F*PART + mvnrnd(zeros(dim,1),Q_CV,ns)';
    PART_sub = PART([1 3],:);
    
    %predict measurements
    tmp_num = MEAS(t).num;
    h = zeros(tmp_num,ns);
    for i=1:tmp_num
        tmp_loc = MEAS(t).loc;
        tmp_dif = PART_sub - repmat(tmp_loc(:,i),1,ns);
        h(i,:) = sqrt( sum(tmp_dif.^2,1) );
    end

    %update importance weights
    tmp = repmat(data_measurement,1,ns) - h;
    weig = -.5 * sum(tmp.*(inv(R)*tmp) , 1);
    weig = weig - max(weig);
    weig = exp(weig);
    weig = weig/sum(weig);
        
    %state estimation
    x_filter(:,t) = sum( bsxfun(@times,weig,PART) , 2);
 
    %standard resampling using population distribution
    flag_resample=1;
    if flag_resample ==1
        PART_new = zeros(size(PART));
        weig_new = zeros(1,ns);
        
        weig_cdf = cumsum(weig);
        
        i = 1;
        u = zeros(1,ns);
        u(1) = rand(1)/ns;
        for j=1:ns
            u(j) = u(1) + (j-1)/ns;
            while u(j) > weig_cdf(i)
                i=i+1;
            end
            PART_new(:,j) = PART(:,i);
            weig_new(j) = 1/ns;
        end
        PART = PART_new;
        weig = weig_new; %not necessary to compute w_{t+1}
    end
end

OUT=x_filter;
end