function result = Integration

Define_Constants ;

% ****** Step 1 : matrix for the GNNS positions ********
matrixGNNSResult = computeGNSSPosition ;
number_epoch = size(matrixGNNSResult,1) ;

matrixGNNSResult(1,3:4) =zeros(1,2) ;

% ****** Step 2 : matrix for the dead reckoning solutions ********
odometry_solutions = deadReckoningSolution ;

% ******* Cmpute the variance between both data *******
pos_var = (odometry_solutions(:,1:2)'*matrixGNNSResult(:,1:2))./number_epoch ;

vel_var = (odometry_solutions(:,3:4)'*matrixGNNSResult(:,3:4))./number_epoch ;

head_var=  (odometry_solutions(:,5)'*matrixGNNSResult(:,5))./number_epoch ;

measurement_noise_matrix = [pos_var, zeros(2), zeros(2,1) ; 
                            zeros(2),vel_var,  zeros(2,1)
                            zeros(1,2), zeros(1,2),head_var*deg_to_rad] ;

% ******* Step 3 : the kalman filter  *************
[x_est,P_matrix] = Initialise_Integration_KF ;

final_values = zeros(number_epoch, 6) ;
state_values = zeros(number_epoch, 5) ;
time = 0 ; 
for i=1:number_epoch
    % odometry result
    this_odometry_sol = odometry_solutions(i,:) ;
    % GNNS result
    GNNSCurrentValue = matrixGNNSResult(i,:) ;
        
    [newState, new_error_cov] = kalmanFilter(x_est, P_matrix,this_odometry_sol,GNNSCurrentValue, measurement_noise_matrix) ;
    
    state_values(i,:) = newState(1:5)' ;
    state_values(i,5) = state_values(i,5)*rad_to_deg ;
    final_values(i,1) = time ; 
    final_values(i,2:6) = this_odometry_sol - state_values(i,:) ;
    
    x_est = newState ;
    P_matrix = new_error_cov ;
    time = time+0.5 ;
end

result = final_values ;
dlmwrite('resultFinal.csv', result, 'delimiter', ',', 'precision', 9); 
end

function [newState, new_error_cov] = kalmanFilter(x_est, P_matrix,this_odometry_sol,GNNSCurrentValue, measurement_noise_matrix)

Define_Constants ;
prog_inter = 0.5 ;

% ***** step 1 : transition matrix *********
transMatrix = [eye(2), zeros(2), zeros(2,1)  , zeros(2,1), zeros(2,1) ; 
               zeros(2), eye(2),zeros(2,1), zeros(2,1), zeros(2,1)   ;    
               zeros(1,2), zeros(1,2), 1, prog_inter, 0 ;
               zeros(1,2), zeros(1,2), 0, 1, prog_inter
               zeros(1,2), zeros(1,2), 0, 0, 1] ;

% ***** step 2 : System noise covariance matrix *********
S_velocity = (0.07^2)*3 ;
S_heading = (2*deg_to_rad)^2 ;
S_gyro = (3*10^4)^2/prog_inter ;
S_bgd = ((1*deg_to_rad/prog_inter)^2)/prog_inter ;

systemNoiseCov = [zeros(2), zeros(2), zeros(2,1)  , zeros(2,1), zeros(2,1)  ;...
                  zeros(2), eye(2)*S_velocity, zeros(2,1), zeros(2,1), zeros(2,1)  ; ...
                  zeros(1,2), zeros(1,2), S_heading, 0,0 ; ...
                  zeros(1,2), zeros(1,2), 0, S_gyro*prog_inter, 0 ;
                  zeros(1,2), zeros(1,2),0 , 0, S_bgd*prog_inter] ;

% ***** step 3&4 : Propagate state and covariance *********
x_plus_est = transMatrix*x_est ;

error_cov =  transMatrix*P_matrix*transMatrix.'+systemNoiseCov ;

% ***** step 5 : Compute the measurement matrix *********
% measurement matrix
mesur_matrix = [-eye(2), zeros(2),  zeros(2,1), zeros(2,1), zeros(2,1) ; ...
                zeros(2), -eye(2), zeros(2,1)  , zeros(2,1), zeros(2,1); 
                zeros(1,2), zeros(1,2), -1, 0, 0] ;
            
% ***** step 6 : Measurement noise covariance matrix *********
%we computed it before with the variance between measurement 
noise_Cov = measurement_noise_matrix ;

% ***** step 7 :  kalman gain matrix *********
kalman_Gain = error_cov*mesur_matrix.'*(mesur_matrix*error_cov*mesur_matrix.'+noise_Cov)^-1 ;

% ***** step 8 :  Compute the measurement innovation *********
measInv = [GNNSCurrentValue(1:2)'-this_odometry_sol(1:2)' + x_plus_est(1:2) ;
           GNNSCurrentValue(3:4)'-this_odometry_sol(3:4)' + x_plus_est(3:4)
           GNNSCurrentValue(5)*deg_to_rad-this_odometry_sol(5)*deg_to_rad + x_plus_est(5)] ;

newState = x_plus_est+kalman_Gain*measInv ;
new_error_cov = (eye(7)- kalman_Gain*mesur_matrix)*error_cov ;
end

function [x_est,P_matrix] = Initialise_Integration_KF
deg_to_rad = 0.01745329252;

% Initialise state estimates
x_est = zeros(7,1);

% Initialise error covariance matrix
P_matrix =  zeros(7);
P_matrix(1:2,1:2) = eye(2) * 0.005^2; % initialize our position error to 0.005 rad
P_matrix(3:4,3:4) = eye(2) * 0.1^2; % initialize our position error to 0.1m.s-1
P_matrix(5,5) = (2*deg_to_rad)^2; % initialize the heading to 2 degree
P_matrix(6,6) = (1*deg_to_rad)^2 ; % 1deg/s
P_matrix(7,7) = (1*deg_to_rad)^2 ; %1 degree par seconde 
end
