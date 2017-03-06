function r = last_test_integration
matrixGNNSResult = computeGNSSPosition ;
number_epoch = size(matrixGNNSResult,1) ;
deltaT = 0.5 ;

first_pos = [0;0] ;
position = zeros(number_epoch,2) ;
position(1,:)= first_pos ;

heading = zeros(number_epoch,1) ;

for i=2:number_epoch
    heading(i) = atan(matrixGNNSResult(i,4)/ matrixGNNSResult(i,3)) ;
    position(i,:) = position(i-1,:) + matrixGNNSResult(i,4:5).*deltaT ; 
end

odometry_solutions = deadReckoningSol ;

x_est = [0;0;0] ;
P_matrix = zeros(3,3);

result = zeros(number_epoch,3) ;
for i=1:number_epoch
   this_odo = odometry_solutions(i,:) ;
   thisGNNS = [heading(i), position(i,:)] ;
    
   [newState, new_error_cov] = kalman_filter_odometry(x_est, P_matrix, this_odo,thisGNNS) ;
    
   
   result(i,:) = newState' ;
   
    x_est = newState ;
    P_matrix = new_error_cov ;
end


end


function [newState, new_error_cov] = kalman_filter_odometry(x_est, P_matrix, this_odo,thisGNNS)
Define_Constants ;

thisHead = this_odo(1) ;
thisOmega = this_odo(2) ;
thisVel =this_odo(3:4) ;
% ***** step 1 : transition matrix *********
% ts is the propagation interval
prop_time = 0.5 ;
transMatrix = [1 0 (-sind(thisHead)*thisVel(1)-cosd(thisHead)*thisVel(2))*prop_time ; ...
    0 1 (cosd(thisHead)*thisVel(1)+sind(thisHead)*thisVel(2))*prop_time ; ...
    0 0 1] ;

% ***** step 2 : System noise covariance matrix *********

systemNoiseCov = eye(3)*0.5 ;

% ***** step 3&4 : Propagate state and covariance *********
x_plus_est = transMatrix*x_est ;

error_cov =  transMatrix*P_matrix*transMatrix.'+systemNoiseCov ;

% ***** step 5 : Compute the measure matrix *********
thisGNNSPos = thisGNNS(2:3) ;
thisGNNSHead = thisGNNS(1) ;

% mesurement matrix
mesur_matrix = [-cosd(x_plus_est(3)), -sind(x_plus_est(3)) , -(x_plus_est(1)-thisGNNSPos(1))*sind(x_plus_est(3))+(x_plus_est(2)-thisGNNSPos(2))*cosd(x_plus_est(3)); ...
                 sind(x_plus_est(3)) -cosd(x_plus_est(3)) , -(x_plus_est(1)-thisGNNSPos(1))*cosd(x_plus_est(3))-(x_plus_est(2)-thisGNNSPos(2))*sind(x_plus_est(3)) ; ...
                 0 0 -1] ;

% ***** step 6 : Measurement noise covariance matrix *********
sigmaP = 2 ; % code tracking+multiplat error standard deviation 
sigmaR = 0.02 ; % range rate tracking error standard deviation 

% computation of measurement noise covariance matrix
noise_Cov = eye(3) ;


% ***** step 7 :  kalman gain matrix *********
kalman_Gain = error_cov*mesur_matrix.'*(mesur_matrix*error_cov*mesur_matrix.'+noise_Cov)^-1 ;

% ***** step 8 :  Compute the measurement innovation *********
hx = [(x_plus_est(1)-thisGNNSPos(1))*cosd(x_plus_est(3))+(x_plus_est(2)-thisGNNSPos(2))*sind(x_plus_est(3)) ; ...
       -(x_plus_est(1)-thisGNNSPos(1))*sind(x_plus_est(3))+(x_plus_est(2)-thisGNNSPos(2))*cosd(x_plus_est(3)) ; ...
       thisGNNSHead-x_plus_est(3)] ;
measurement_innov = thisGNNS' - hx ;

% ***** step 9 :  Compute the new state *********
newState = x_plus_est+kalman_Gain*measurement_innov ;


% ***** step 10 :  Compute the new covariance matrix *********
new_error_cov = (eye(3)- kalman_Gain*mesur_matrix)*error_cov ;

end

function matrixGNNSResult = computeGNSSPosition
% load the result for the GNSS position

Define_Constants ;

% load the files for the GNSS
pseudo_range_rate = load('Pseudo_range_rates.csv') ;
pseudo_range = load('Pseudo_ranges.csv') ;

% how many satellite they are and there number
satellites = pseudo_range(1,2:end) ;
number_Sat_Total = size(satellites,2) ;

% will store the position  of the satellitesfor epoch 0
satellite_Positions_Vel = zeros(3, number_Sat_Total,2) ;
for i=1:number_Sat_Total
    this_satellite = satellites(1,i) ;
    [position_Sat,velocity_Sat] = Satellite_position_and_velocity(0,this_satellite);
    satellite_Positions_Vel(:,i,1) = position_Sat ;
    satellite_Positions_Vel(:,i,2) = velocity_Sat ;

end

% initialization of the positions and matrix
[x_est,P_matrix] = Initialise_First_Position(pseudo_range(2,2:end),satellite_Positions_Vel, pseudo_range_rate(2,2:end)) ;

% how many epochs they are
number_epoch = size(pseudo_range_rate,1) ;

% will store the results of the GNSS (excluding the line with the number of
% satellites
matrixGNNSResult = zeros(number_epoch-1,6) ;

% store the first result
[L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(x_est(1:3),x_est(4:6)) ;
latPredic = L_b*rad_to_deg ;
longPredic = lambda_b*rad_to_deg ;
matrixGNNSResult(1,:) = [latPredic,longPredic,h_b, v_eb_n'] ;

for i=2:number_epoch-1
    % measured pseudo ranges and pseudo range rates
    measured_ranges = pseudo_range(1+i,2:end) ;
    measured_ranges_rates = pseudo_range_rate(1+i,2:end) ;
    
    % give the time for this epoch
    clockCurrent = pseudo_range(i+1,1) ;
    
    % compute the new state and cov matrix using a kalman filter
    [newState, new_error_cov, outliers] = kalman_filter(x_est, P_matrix, clockCurrent, measured_ranges,measured_ranges_rates, satellites) ;
    
    % if the index is true, it means outlier have been detected, we take
    % them out the measurement and compute again the position
    if ~isempty(outliers)
        [measured_range_new,measured_ranges_rates_new, satellites_new] = noOutliers(measured_ranges,measured_ranges_rates, satellites, outliers) ;
        [newState, new_error_cov, ~] = kalman_filter(x_est, P_matrix, clockCurrent, measured_range_new,measured_ranges_rates_new, satellites_new) ;
    end
    
    [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(newState(1:3),newState(4:6)) ;
    latPredic = L_b*rad_to_deg ;
    longPredic = lambda_b*rad_to_deg ;
    
    % store the result for position and velocity for north and east
    matrixGNNSResult(i,:) = [latPredic,longPredic,h_b, v_eb_n'] ;
    
    % update the filter to continue computing the position
    x_est = newState ;
    P_matrix = new_error_cov ;
end
end


function [newState, new_error_cov, outliers] = kalman_filter(x_est, P_matrix, clockCurrent, measured_ranges,measured_ranges_rates, numS)

Define_Constants ;

% ***** step 1 : transition matrix *********
% ts is the propagation interval
prop_time = 0.5 ;
transMatrix = [eye(3), eye(3)*prop_time, zeros(3,1), zeros(3,1) ; ...
    zeros(3), eye(3), zeros(3,1), zeros(3,1); ...
    zeros(1,3),zeros(1,3), 1, prop_time;...
    zeros(1,3),zeros(1,3), 0, 1] ;

% ***** step 2 : System noise covariance matrix *********
% acceleration power spectral density
Sae = 1 ; % value of the acceleration for a pedestrian 
% clock phase
Sacalpha = (1)^2/prop_time ; %(signalInSpace + residualionosphere+residualTrop)^2/prog_time
% clock frequency
Sacf = (200)^2/prop_time ; %(receiverClockSD)^2/prog_time

systemNoiseCov = [1/3*Sae*prop_time^3*eye(3), 1/2*Sae*prop_time^2*eye(3), zeros(3,1), zeros(3,1) ;...
    1/2*Sae*prop_time^2*eye(3), Sae*prop_time*eye(3), zeros(3,1), zeros(3, 1); ...
    zeros(1,3), zeros(1,3), Sacalpha*prop_time+1/3*prop_time^3*Sacf, 1/2*prop_time^2*Sacf ; ...
    zeros(1,3), zeros(1,3),  1/2*prop_time^2*Sacf, Sacf*prop_time] ;

% ***** step 3&4 : Propagate state and covariance *********
x_plus_est = transMatrix*x_est ;

error_cov =  transMatrix*P_matrix*transMatrix.'+systemNoiseCov ;

% ***** step 5 : Compute the measure matrix *********

% before we need to compute the predictited range and the line of sight

% number of satellites
number_sat = size(numS,2) ;
% compute and store the position and the velocity for each satallite at this epoch
satellite_positions_vellocity = zeros(3, number_sat, 2) ;
for i=1:number_sat
    this_sat = numS(1,i) ;
    [positionSat, sat_v_es_e] = Satellite_position_and_velocity(clockCurrent,this_sat);
    satellite_positions_vellocity(:,i,1) = positionSat ;
    satellite_positions_vellocity(:,i,2) = sat_v_es_e ;
end

% compute the predictided ranges for this epoch
predictRange = zeros(1,number_sat) ;
for i=1:number_sat
    % compute a first time the range without a sagnac correction
    satPos = satellite_positions_vellocity(:,i,1) ;
    x = satPos-x_plus_est(1:3,:) ;
    range1 = sqrt(x.'*x) ;
    
    % compute the sagnac compensation matrix with this previous predicted range
    C = [1 omega_ie*range1/c 0 ; -omega_ie*range1/c 1 0 ; 0 0 1] ;
    
    % compute and store the new range
    x2 = C*satPos-x_plus_est(1:3,:) ;
    predictRange(:,i) = sqrt(x2.'*x2);
end

% compute the line of sight for each satellite
lineOfSight = (satellite_positions_vellocity(:,:,1)-repmat(x_plus_est(1:3,:), [1, number_sat]))./repmat(predictRange, [3,1]) ;

% mesurement matrix
mesur_matrix = [-lineOfSight.' zeros(number_sat, 3) ones(number_sat,1) zeros(number_sat,1) ; ...
    zeros(number_sat, 3) -lineOfSight.' zeros(number_sat,1) ones(number_sat,1)] ;

% ***** step 6 : Measurement noise covariance matrix *********
sigmaP = 2 ; % code tracking+multiplat error standard deviation 
sigmaR = 0.02 ; % range rate tracking error standard deviation 

% computation of measurement noise covariance matrix
noise_Cov = [eye((number_sat),(number_sat)).*sigmaP^2 zeros(number_sat,number_sat) ; zeros(number_sat,number_sat) eye(number_sat,number_sat).*sigmaR^2] ;
noise_Cov = noise_Cov + [zeros(number_sat,2*(number_sat)) ; zeros(number_sat,2*(number_sat)-1) ones(number_sat,1)] ;
noise_Cov(2*(number_sat),2*(number_sat)) = noise_Cov(2*(number_sat),2*(number_sat))-1 ;


% ***** step 7 :  kalman gain matrix *********
kalman_Gain = error_cov*mesur_matrix.'*(mesur_matrix*error_cov*mesur_matrix.'+noise_Cov)^-1 ;

% compute the predicted range rates
rangeRatePredicted = zeros(1,number_sat) ;
for i=1:number_sat
    % the position and velocity for the satellite i
    satPos = satellite_positions_vellocity(:,i,1) ;
    satVel = satellite_positions_vellocity(:,i,2) ;
    
    % position and velocity of our propagated state
    pos = x_plus_est(1:3,:) ;
    vel = x_plus_est(4:6,:) ;
    
    % predicted range value for this satellite
    rangeValue = predictRange(1,i) ;
    
    % predict the range using the sagnac correction matrix
    xL = satVel+Omega_ie*satPos ;
    xR = vel+Omega_ie*pos ;
    u = lineOfSight(:,i).' ;
    C = [1, omega_ie*rangeValue/c, 0 ; -omega_ie*rangeValue/c 1 0 ; 0 0 1] ;
    
    rangeRatePredicted(1,i) = u*(C*xL-xR);
end

% ***** step 8 :  Compute the measurement innovation *********
measurement_innov = zeros(1, (number_sat)*2) ;
measurement_innov(1,1:number_sat) = measured_ranges - predictRange - ones(1,number_sat)*x_plus_est(end-1) ;
measurement_innov(1,number_sat+1:end) = measured_ranges_rates - rangeRatePredicted -ones(1,number_sat)*x_plus_est(end) ;

% ***** step 9 :  Compute the new state *********
newState = x_plus_est+kalman_Gain*measurement_innov.' ;


% ***** step 10 :  Compute the new covariance matrix *********
new_error_cov = (eye(8)- kalman_Gain*mesur_matrix)*error_cov ;


% ***** Extra step :  Look for outliers using only the position *********
% compute the measurement matrix and the measurement innovation vector with only the positions
mesur_matrix_positions = [mesur_matrix(1:number_sat , 1:3), ones(number_sat, 1)] ;
measurement_innov_position = measurement_innov(1,1:number_sat) ;
% compute the residual vector
residual = (mesur_matrix_positions*((mesur_matrix_positions.'*mesur_matrix_positions)^-1)*mesur_matrix_positions.'-eye(number_sat))*measurement_innov_position.' ;

% compute the measurement error standart deviation for each satellite in
% function of it's line of sight
measErrorSDZenith = repmat(2.2, [1,number_sat ]) ; % residual troposphere + ionosphere error
%compute the standard deviation and square it to obtain the variance
measErrorVariance = (measErrorSDZenith./lineOfSight(3,:)).^2 ;

% compute the residual covariance matrix (without the measurement error SD)
residual_cov_M = eye(number_sat)-mesur_matrix_positions*(mesur_matrix_positions.'*mesur_matrix_positions)^-1*mesur_matrix_positions.' ;

% multiply the diag of the residual covariance matrix with the measurement
% error variance for each satellite
diag_residual_cov_M = diag(residual_cov_M)'.*measErrorVariance ;

% compute the normalized residual
normalized_residual = residual'./sqrt(diag_residual_cov_M) ;

% compute a test with a threshold value of 6 to detect outliers
test = abs(normalized_residual) > (diag_residual_cov_M).*6 ;

if sum(test) > 0
    outliers = test ;
else
    outliers = [] ;
end

end


function [predicted_ranges,predicted_ranges_rates, satellites] = noOutliers(predicted_ranges_old,predicted_ranges_rates_old, satellites_old, outliers)
sumber_Sat = size(satellites_old,2) ;
predicted_ranges = predicted_ranges_old ;
predicted_ranges_rates = predicted_ranges_rates_old ;
satellites = satellites_old ;

% if a satellites is an outliers, we suppress its values
index = 0 ;
for i=1:sumber_Sat
    if outliers(i)==1
        index = index+1 ;
        predicted_ranges(i-index) = [] ;
        predicted_ranges_rates(i-index) = [] ;
        satellites(i-index) = [] ;
    end
end

end

function [x_est_final,P_matrix] = Initialise_First_Position(ranges,satPosVel, ranges_rates)

% initialize the first position of the user using least square approximation
% and use the first satellite position to linearize the problem
satPos = satPosVel(:,:,1) ;
numbSat  = size(satPos, 2) ;
Define_Constants ;
x_est = zeros(8,1) ;
% compute corrected position from sagnac effect
correctedPos = zeros(size(satPos)) ;

for i=1:numbSat
    currentRange = ranges(i) ;
    currentPos = satPos(:,i) ;
    sagnac_comp = [1 omega_ie*currentRange/c 0 ; -omega_ie*currentRange/c 1 0 ; 0 0 1] ;
    correctedPos(:,i) = sagnac_comp*currentPos ;
end

% compute the positions
correctedPos1 = correctedPos(:,1) ;
correctedPos(:,1) = [] ;
correctedPos1 = repmat(correctedPos1, [1, numbSat-1]) ;

% compute the ranges
range1 = ranges(1) ;
allRanges = ranges ;
ranges(1) = [] ;
range1 = repmat(range1, [1, numbSat-1]) ;

% compute the A matrix for least square approximation
A = (correctedPos-correctedPos1)' ;

% compute the distance between the first one and the other satellite
dist = sum((correctedPos(:,:)-correctedPos1(:,:)).^2, 1) ;

b = 1/2.*(range1.^2-ranges.^2+dist).' ;

x_result = A\b ;
x_est(1:3) = x_result + correctedPos1(:,1) ;

% I assumed the first velocite of the user to be zeros
x_est(4:6) = zeros(3,1) ;

% I now compute the clock offset and drift using single epoch positionning
clockOffsetStart = 1 ;
clockDriftCurrentStart = 0 ;

x_est(7:8) = [clockOffsetStart; clockDriftCurrentStart] ;

x_est_final = findClockValues(x_est, allRanges,satPosVel, ranges_rates) ;

% Initialise error covariance matrix
% as I initialized my value using a single epoch, I'm quite confidant with
% their values 
P_matrix =  zeros(8);
P_matrix(1,1) = x_est_final(7)^2 ;
P_matrix(2,2) = x_est_final(7)^2 ;
P_matrix(3,3) = x_est_final(7)^2 ;
P_matrix(4,4) = x_est_final(8)^2 ;
P_matrix(5,5) = x_est_final(8)^2 ;
P_matrix(6,6) = x_est_final(8)^2 ;
P_matrix(7,7) = x_est_final(7)^2 ;
P_matrix(8,8) = x_est_final(8)^2 ;
end


function new_state_final = findClockValues(x_est, ranges,satPosVel, ranges_rates)

Define_Constants ;

number_sat = size(ranges,2) ;
satellitePositions = satPosVel(:,:,1) ;
satelliteVelocity = satPosVel(:,:,2) ;

% use 4 iterations to assure a convergence 
while true
    % ********* Step 1 : compute the clock offset value ************
    % compute the predicted ranges for each satellite
    predictec_ranges= zeros(1,number_sat) ;
    for i=1:number_sat
        thisSatPos = satellitePositions(:,i) ;
        x = thisSatPos-x_est(1:3) ;
        range1 = sqrt(x.'*x) ;
        C = [1 omega_ie*range1/c 0 ; -omega_ie*range1/c 1 0 ; 0 0 1] ;
        x2 = C*thisSatPos-x_est(1:3) ;
        
        predictec_ranges(:,i) = sqrt(x2.'*x2);
    end
    
    % computation of the line of sight
    lineOfSight = (satellitePositions-repmat(x_est(1:3), [1, number_sat]))./repmat(predictec_ranges, [3,1]) ;
    
    % computation of the new position
    measurement_innov = ranges - predictec_ranges -ones(1,number_sat).*x_est(4) ;
    measurement_M = [-lineOfSight.' ones(number_sat,1)] ;
    new_state_pos = x_est(1:4) +(measurement_M.'*measurement_M)^-1*measurement_M.'*measurement_innov.' ;
    
    % ********* Step 2 : Compute the clock drift value ************
    % compute the predicted range rates values 
    predicted_range_rates= zeros(1,number_sat) ;
    for i=1:number_sat
        thisSatPos = satellitePositions(:,i) ;
        thisSatVel = satelliteVelocity(:,i) ;
        
        thisRange = ranges(1,i) ;
        
        xL = thisSatVel+Omega_ie*thisSatPos ;
        xR = x_est(5:7)+Omega_ie*x_est(1:3) ;
        u = lineOfSight(:,i).' ;
        C = [1, omega_ie*thisRange/c, 0 ; -omega_ie*thisRange/c 1 0 ; 0 0 1] ;
        
        predicted_range_rates(1,i) = u*(C*xL-xR);
    end
    
    % computation of the new velocity
    measurement_innov_vel = ranges_rates - predicted_range_rates - ones(1,number_sat)*x_est(8) ;
    measurement_M_vel = [-lineOfSight.' ones(number_sat,1)] ;
    new_state_vel = x_est(5:8) +(measurement_M_vel.'*measurement_M_vel)^-1*measurement_M_vel.'*measurement_innov_vel.' ;
    
    x_est_new = [new_state_pos ; new_state_vel];
    
    test = sqrt(sum((x_est-x_est_new).^2,1)) ;
    if test < 0.001
       break ;
    else
        x_est = x_est_new ;
    end
end
new_state_final = [new_state_pos;new_state_vel] ;
end

function odometry_solutions = deadReckoningSol

Define_Constants ;

% *********** step 1 : Load files and initialization ***********

% load the file
dead_reckoning_info = load('Dead_reckoning.csv') ;
deltaT = 0.5 ; % propagation time

heading =  dead_reckoning_info(:,7) ;
number_epoch = size(dead_reckoning_info,1) ;
matrixGNNSResult = computeGNSSPosition ;
odometry_solutions = zeros(number_epoch, 4) ;
odometry_solutions(:,1) = heading ;

% *********** step 2 : Compute the rear wheel-frame  speed ***********
wheel_speed_3 = dead_reckoning_info(:,4) ;
wheel_speed_4 = dead_reckoning_info(:,5) ;
vel_wheels = 1/2.*(wheel_speed_3+wheel_speed_4) ;

% *********** step 3 : Correct w_ib_b and vel_wheels from the scale factors and cross coupling errors ***********
w_ib_b_value = dead_reckoning_info(:,6) ;
w_ib_b = zeros(number_epoch,3) ;

% correction of the wheel speed sensors
vel_wheels = vel_wheels.*3 ; %scale factor sd of 3%

for i=1:number_epoch
    % correction of the gyroscope sensor
    Mg = [1, 0.1, 0.1 ;
          0.1, 1, 0.1 ;
          0.1, 0.1, 1] ;
    rot_vec = [0; 0; w_ib_b_value(i)] ;
    w_ib_b(i,:) = Mg*rot_vec ;
end

odometry_solutions(:,2) = w_ib_b_value ;

% *********** step 4 : Computing the v_eb_b ***********
center_lawnmower_array = repmat([-0.2,0,0], [number_epoch,1]) ;
comp1 = [vel_wheels, zeros(number_epoch,1),zeros(number_epoch,1)] ;
comp2 = cross(w_ib_b,center_lawnmower_array) ;
v_eb_b = comp1-comp2 ;

% *********** step 4 : Computing the v_eb_n  ***********
v_eb_n = zeros(number_epoch,2) ;
delta_r_eb_n = zeros(number_epoch,2) ;

for i = 1:number_epoch

     mat = [ cosd(heading(i)), -sind(heading(i)) ;
             sind(heading(i)), cosd(heading(i))] ;
     
     v_eb_n(i,:) = mat*v_eb_b(i,1:2)' ;
     
     delta_r_eb_n(i,:) = v_eb_n(i,:)*deltaT ;
end

odometry_solutions(:,3:4) = v_eb_n ;
end
