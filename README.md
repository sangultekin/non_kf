1. PROJECT DESCRIPTION

This Github repository contains the Matlab scripts, and a sensor network demo for the paper:

S. Gultekin and J. Paisley. Nonlinear Kalman filtering with divergence minimization, IEEE Transactions on Signal Processing, vol. 65, no. 23, pp. 6319-6331.

The code is written by:
San Gultekin
Columbia University
san.gultekin@gmail.com

The scripts are all tested in Matlab R2016b and only depend on the core Matlab libraries, and should work without hassle. In particular the file paths are kept relative and usage of slash/backslash is avoided to have compatibility with both Windows and Unix-based operating systems. The processor and memory requirements are low and the code should run in reasonable amount of time in a mediocre computer. If, however, you run into any issues or bugs please send and email.

WARNING: The errors reported are mean errors (ME) but for these scripts we use root mean square error (RMSE), which corresponds to swapping the mean and sqrt functions in the outer part of the code block starting Line 117. For this reason the numbers here are slightly different than those in the paper, but all analysis and comparisons remain the same. (Note that both metrics are decomposing total sum of squared error, which is the same value independent of the decomposition aproach choosen, i.e. ME or RMSE.) We do so here, because RMSE is a more common metric and probably a better choice for the published code.



2. CONTENTS OF THIS REPOSITORY

- sample_data: The folder containing 100 sensor network measurements. Each file contains:
	- data_actual: correct target location
	- data_full: correct target location and velocity (full state)
	- loc_sns: location of sensors scattered across the field
	- target_state: initial state
	- MEAS: A struct with length equal to the total number of measurements, containing number of measurements (num), sensor identified (id), sensor location (loc), distance only measurements (z_true), and measurement covariance which is num-by-num (R).

- aux_jac: script to calculate Jacobian of measurement nonlinearity

- aux_predict_meas: script to calculate the corresponding measurement of a target state

- sample_data_constants: define state-space parameters for the provided dataset

- func_EKF: Extended Kalman Filter

- func_UKF: Unscented Kalman Filter

- func_SIR: Sequential Importance Resampling (SIR) Particle Filter

- func_NKF: Ensemble Kalman Filter

- func_SKF: Stochastic Search Kalman Filter (Algorithm 1) from the paper cited above.

- func_AKF: Alpha Divergence Kalman Filter (Algorithm 3) from the paper cited above. Also, setting alpha=1 here gives Moment Matching Kalman Filter (Algorithm 2).

- main_script: Run all the filters and calculate RMSE. This is the only script you need to run and evaluate the results.



3. RUNNING THE SCRIPTS

Open the file main_script in Matlab and execute it with the Run command. All filters will be run and evaluated simultaneously, and progress will be printed to the console.

WARNING: SKF is an iterative procedure and has a double loop, one for each time point, and one for optimization of posteriors. For longer simulations it can slow down the progress. In these cases, you may want to first comment it out during initial runs to determine optimal parameters for your setup, and add the SKF in the final run.

After running sample_data_constants at Line 5, the parameters that will be given to the filters are set. These are the optimal values. If you want to simulate parameter uncertainty as we did in our paper, you can overwrite the parameter value. In the default main_script note that at Line 9 we write 

Q_CV = 1e-1*eye(4)

This overwrites the true value 

var_CV = 1e-2;
Q2 = [T^4/4 T^3/2 ; T^3/2 T^2];
Q_CV = var_CV * [Q2 zeros(2,2) ; zeros(2,2) Q2]; 

in the sample_data_constants script.

Simulating with this parameter uncertainty and setting Q = 1e-1*eye(4) we get:

Root Mean Square Error Results
EKF RMSE: 16.9535
UKF RMSE: 11.9780
SIR RMSE: 10.9289
NKF RMSE: 13.2884
SKF RMSE: 10.9259
MKF RMSE: 10.5043
AKF RMSE: 9.1717

This is comparable to Table II of the paper. When parameter uncertainty is present, SIR is no longer the best state estimator. In particular, AKF does well in this case, as the alpha parameter has a "likelihood dampening" effect as we argue in the paper.

-----

Of course we can also simulate with matching parameters. Now simply comment out Line 9. We get the following:

Root Mean Square Error Results
EKF RMSE: 16.1015
UKF RMSE: 13.1994
SIR RMSE: 10.8891
NKF RMSE: 12.1924
SKF RMSE: 12.0647
MKF RMSE: 12.0247
AKF RMSE: 11.1030

This time the result agrees with Fig 4(c-d) of the paper. When no parameter uncertainty is present, particle filter does the best as expected. AKF is the only remaining filter that can achieve competitive results for this case, followed by SKF and MKF.