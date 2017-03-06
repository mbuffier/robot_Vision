function [odometry_solutions] = deadReckoningSolution

Define_Constants ;

% *********** step 1 : Load files and initialization ***********

% load the file
dead_reckoning_info = load('Dead_reckoning.csv') ;
deltaT = 0.5 ; % propagation time

heading =  dead_reckoning_info(:,7) ;
number_epoch = size(dead_reckoning_info,1) ;
matrixGNNSResult = computeGNSSPosition ;
odometry_solutions = zeros(number_epoch, 5) ;
odometry_solutions(:,5) = heading ; % heading en degree 

% *********** step 2 : Compute the rear wheel-frame  speed ***********
wheel_speed_3 = dead_reckoning_info(:,4) ;
wheel_speed_4 = dead_reckoning_info(:,5) ;
vel_wheels = 1/2.*(wheel_speed_3+wheel_speed_4) ;

% *********** step 3 : Correct w_ib_b from the scale factors and cross coupling errors ***********
w_ib_b_value = dead_reckoning_info(:,6) ;
w_ib_b = zeros(number_epoch,3) ;

for i=1:number_epoch
    % correction of the gyroscope sensor
    Mg = [1, 0.1, 0.1 ;
          0.1, 1, 0.1 ;
          0.1, 0.1, 1] ;
    rot_vec = [0; 0; w_ib_b_value(i)] ;
    w_ib_b(i,:) = Mg*rot_vec ;
end

% *********** step 4 : Computing the v_eb_b ***********
center_lawnmower_array = repmat([-0.2,0,0], [number_epoch,1]) ;
comp1 = [vel_wheels, zeros(number_epoch,1),zeros(number_epoch,1)] ;
comp2 = cross(w_ib_b,center_lawnmower_array) ;
v_eb_b = comp1-comp2 ;

% *********** step 5 : Computing the v_eb_n  ***********
v_eb_n = zeros(number_epoch,2) ;
delta_r_eb_n = zeros(number_epoch,2) ;

for i = 1:number_epoch

     mat = [ cosd(heading(i)), -sind(heading(i)) ;
             sind(heading(i)), cosd(heading(i))] ;
     
     v_eb_n(i,:) = mat*v_eb_b(i,1:2)' ;
     
     delta_r_eb_n(i,:) = v_eb_n(i,:)*deltaT ;
end

odometry_solutions(:,3:4) = v_eb_n ;

% *********** step 6 : Computing of new lattitude and longitude from the results   ***********
odometry_solutions(1,1:2) = matrixGNNSResult(1,1:2) ;
height = 0 ;
for i = 2:number_epoch
    [R_N,R_E]= Radii_of_curvature(odometry_solutions(i-1,1)) ;
    odometry_solutions(i,1) = odometry_solutions(i-1,1) + (delta_r_eb_n(i,1) / R_N+height)*rad_to_deg ;
    odometry_solutions(i,2) = odometry_solutions(i-1,2) + (delta_r_eb_n(i,2) / ((R_E+height)*cosd(odometry_solutions(i-1,1))))*rad_to_deg ;
end

end