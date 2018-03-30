%define state-space parameters of the synthetic data

%Author: San Gultekin, Columbia University, san.gultekin@gmail.com



T=1; n=300; MC_number=100; 

F2 = [1 T ; 0 1];
F = [F2 zeros(2,2) ; zeros(2,2) F2];

var_CV = 1e-2;
Q2 = [T^4/4 T^3/2 ; T^3/2 T^2];
Q_CV = var_CV * [Q2 zeros(2,2) ; zeros(2,2) Q2];  

mag2 = (20)^2;

no_sns = 200;
range_sns = 1e6;
min_sns = 3;
max_sns = 3;