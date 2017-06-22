%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Estimate pendulum state with UKF and URTS as in Examples 5.3
% and 9.3 of the book
%
% Simo Sarkka (2013), Bayesian Filtering and Smoothing,
% Cambridge University Press. 
%
% Last updated: $Date: 2013/08/26 12:58:41 $.
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Data Generation
%
pendulum_sim;
%
% Filtering
%
%% AdUKF 1
Xprev = m0;
Pxx_prev = P0;

MM = zeros(size(Xprev,1),length(Y));
PP = zeros(size(Pxx_prev,1),size(Pxx_prev,2),length(Y));
f = @StateFunction;
h = @MeasurementFunction;
for k=1:length(Y)
    y = Y(k);
    [m, P] = AdUKF(Xprev, Pxx_prev, f, param_F, h, [], Q, R, y, 1, 0.1);
    Xprev = m;
    Pxx_prev = P;
    MM(:,k) = m;
    PP(:,:,k) = P;
end  

h = plot(T,Y,'k.',T,X(1,:),'r-',T,MM(1,:),'b--');
set(h,'Linewidth',5);
title('UKF estimate');
legend('Measurements','True','Estimate');

% rmse_ukf = sqrt(mean((X(1,:)-MM(1,:)).^2))
    
    
%%
% Plot the filtering result
%

    clf;
    h=plot(T,X(1,:),'k',T,Y,'bo',T,MM(1,:),'r');
    legend('True Angle','Measurements','UKF Estimate');
    xlabel('Time{\it t}');
    ylabel('Pendulum angle {\it{x}}_{1,{\it{k}}}') 
          