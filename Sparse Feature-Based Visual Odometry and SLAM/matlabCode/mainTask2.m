%% *************** Task 2.1 : Generate a 3D plot of your results **************
% perform the algoritm for task 21 and 22
[groundTruth_position,initial_position,output_position, ~, error_location, ATE_result, accuracyArray] = bundle_adjustment(1, 1,1, 15) ;

fig1 = figure ;
pcshow(groundTruth_position, 'g', 'MarkerSize', 15) ;
hold on, 
pcshow(initial_position, 'r', 'MarkerSize', 6) ;
hold on 
pcshow(output_position, 'b', 'MarkerSize', 8) ;
legend('Ground truth', 'Initial position', 'Output position'), title('Evolution of the positions after the optimisation') ;

saveas(fig1,'task21.jpg') ;

%% *************** Task 2.2 : Evaluate the Performance of the Bundle Adjustment Algorithm **************

% evolution of the covariance along the axis 
fig2 = figure ;
set(fig2, 'Position', [0 0 1000 800]),
Y = [1 2:15] ;

% If I had access to the covariances 
% subplot(4,1,1), 
% plot(Y, spacial_uncert, 'r');
% legend('Mean of the magnitude of the standard deviation vector'), title('Spacial uncertainty at each iteration');

subplot(3,1,1), 
error_location = (error_location-min(error_location(:)))./(max(error_location(:)-min(error_location(:)))) ;
plot(Y, error_location, 'r');
legend('Mean of the distance between ground-truth and estimated position of the landmarks'), title('Error in the location of the landmarks at each iteration');

subplot(3,1,2), 
plot(Y, ATE_result, 'r');
legend('ATE for the camera positions'), title('ATE metrics for the camera position at each iteration');

subplot(3,1,3), 
plot(Y, accuracyArray, 'r');
legend('error (with respect to the initial error)'), title('Evolution of the error ');

saveas(fig2,'task22_covAxis.jpg') ;

%%  *************** Task 2.3 : Explore Methods to Reduce the Complexity of the Graph **************
% perform the test for several number of measurement 
numberIter = 30 ;
numberPoint = 250 ;
error_location100 = bundle_adjustment_test(100, numberIter, numberPoint) ;
error_location90 = bundle_adjustment_test(90, numberIter, numberPoint) ;
error_location80 = bundle_adjustment_test(80, numberIter, numberPoint) ;
error_location70 = bundle_adjustment_test(70, numberIter, numberPoint) ;
error_location60 = bundle_adjustment_test(60, numberIter, numberPoint) ;
error_location50 = bundle_adjustment_test(50, numberIter, numberPoint) ;
error_location40 = bundle_adjustment_test(40, numberIter, numberPoint) ;
error_location30 = bundle_adjustment_test(30, numberIter, numberPoint) ;
error_location20 = bundle_adjustment_test(20, numberIter, numberPoint) ;
error_location10 = bundle_adjustment_test(10, numberIter, numberPoint) ;
error_location5 = bundle_adjustment_test(5, numberIter, numberPoint) ;

% create the vector to plot and resize it between 0 and 1
errors = [error_location5 error_location10,error_location20, error_location30 ...
        error_location40,error_location50,error_location60, error_location70 ...
        error_location80,error_location90,error_location100] ;
    
lastError2 = (errors-min(errors(:)))./(max(errors(:))-min(errors(:))) ;

% Plot the result 
X = [5 10:10:100]' ;
fig3 = figure ; 
plot(X, lastError2), title('Evolution of the mean error in location for 250 features in function of the % of mesurements'), 
xlabel('Percentage of measurement created for each landmark'), ylabel('Error resized between 0 and 1');

saveas(fig3,'task23.jpg') ;

% ************** test the result above ********** 

[~,~,~, ~, ~, ATE_result40, ~] = bundle_adjustment_lessMeas(1, 1, 1, 20, 40, 10) ;
[~,~,~, ~, ~, ATE_result60, ~] = bundle_adjustment_lessMeas(1, 1, 1, 20, 60, 10) ;
[~,~,~, ~, ~, ATE_result100, ~] = bundle_adjustment(1, 1,1, 20) ;
[~,~,~, ~, ~, ATE_result10, ~] = bundle_adjustment_lessMeas(1, 1, 1, 20, 10, 10) ;

% curve for ATE results 
ATE_test = [ ATE_result10(20,1) ATE_result40(20,1), ATE_result60(20,1), ATE_result100(20,1)] ;
ATE_test2 = (ATE_test-min(ATE_test(:)))./(max(ATE_test(:))-min(ATE_test(:))) ;

% plot the result 
X2 = [10 40 60 100] ;
fig4=figure ; 
plot(X2, errors_test2), 
legend('ATE results between 0 and 1'), 
xlabel('Percentage of measurement accepted'), title('Comparaison off the ATE results with a different number of measurement') ;



