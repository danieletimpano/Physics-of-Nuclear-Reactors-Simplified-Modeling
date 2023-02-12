close all
clear all
clc

%% Exercise 5 - Fuel evolution
% ### Algorithm description ###

% For the solution of the problem the following procedure was employed:
% SECTION 1 - Data input from TP intro slides
% SECTION 2 - Define starting number of nuclides using information on
%             effective density of the fuel
% SECTION 3 - Build up the equations for each nuclide such to define mat. A
%         - Forward euler and Mat exponential are implemented in for cycles
% SECTION 4 - Analytical solutions are written as function handles
%           - Plots for analytical and numerical solutions are carried out
%           - Relative difference for each nuclide is also computed and
%             plotted
% SECTION 5 - The convergence properties are studied using a for cycle
%             changing the iteration step from 1 min to six hours
%           - The error between the analytical solution and the numerical
%             is stored at time=1 day for different time steps
% SECTION 6 - At every timestep the flux is recomputed accordingly to the
%             energy fission rate obtained after the depletion (i.e. with the
%             new material composition). The matrix A is recomputed at each
%             timestep accordingly to the new flux while the latter is
%             saved at each timestep.
% Please mind that section 6 may take a long time to run due to the small
% timestep (yet necessary to converge) and to the saving procedure
% implemented. The solution, while taking a bit of run time, shall be 
% consistent as discussed in class and was discussed in the report. A
% waitbar is shown in order to highlight the evolution of the simulation.


%% Section 1 - DATA

%U-235
sigma_c_U5= 12*10^(-24);
sigma_s_U5= 10*10^(-24);
nu_U5= 2.44;
sigma_f_U5= 56*10^(-24);
kappa_U5= 3.24E-11;
Thalf_U5= 0;
FY5= 0;
N5_start= [];

%U-238
sigma_c_U8= 4*10^(-24);
sigma_s_U8= 10*10^(-24);
nu_U8= 2.79;
sigma_f_U8= 1*10^(-24);
kappa_U8= 3.32E-11;
Thalf_U8= 0;
FY8= 0;
N8_start= [] ;

%Pu-239
sigma_c_P9= 81*10^(-24);
sigma_s_P9= 10*10^(-24);
nu_P9= 2.87;
sigma_f_P9= 144*10^(-24);
kappa_P9= 3.33E-11;
Thalf_P9= 0;
FY9= 0; 
N9_start=0 ;

%X
sigma_c_X= 5E6*10^(-24);
sigma_s_X= 10*10^(-24);
nu_X= 0;
sigma_f_X= 0*10^(-24);
kappa_X= 0;
Thalf_X= 9*3600; %in seconds
FYX= 0.06;
NX_start=0 ;

%Y 
sigma_c_Y= 50*10^(-24);
sigma_s_Y= 10*10^(-24);
nu_Y= 0;
sigma_f_Y= 0*10^(-24);
kappa_Y= 0;
Thalf_Y= 0;
FYY= 1.94;
NY_start= 0;

Flux= 10^13;
kappa= [kappa_U8, kappa_U5, kappa_P9, kappa_X, kappa_Y];
sigmaf= [sigma_f_U8, sigma_f_U5, sigma_f_P9, sigma_f_X, sigma_f_Y];

%geometry of the assembly
ass_pitch= 20;
height= 400;
r= 0.4;
n_ass=264;
N_A= 6.022*10^23;
rho_UO2= 10.4; %gcm^3
V_fuel= ass_pitch^2*height;
m_fuel= pi*r^2*height*n_ass*rho_UO2;
enrich= 0.03;

%% SECTION 2 - define the starting number of nuclides

m_UO2= m_fuel; 
rho_fuel= m_fuel/V_fuel; %effective density
M_UO2= 235*enrich + (1-enrich)*238 + 16*2;
N_UO2= rho_fuel*N_A/M_UO2;
N5_start= N_UO2*enrich; %use information on enrichment
N8_start= N_UO2*(1-enrich); %use information on enrichment

%% SECTION 3 - Building up the equations and solver

eq238= [-(sigma_c_U8+sigma_f_U8)*Flux, 0, 0, 0, 0];
eq235= [0,-(sigma_c_U5 + sigma_f_U5)*Flux, 0, 0, 0];
eqPu9= [Flux*sigma_c_U8, 0, -Flux*(sigma_c_P9+sigma_f_P9) , 0, 0];
eqX= [0, Flux*FYX*sigma_f_U5, 0, (-sigma_c_X*Flux - log(2)/Thalf_X), 0]; 
eqY= [0, Flux*FYY*sigma_f_U5, 0,0,(-sigma_c_Y*Flux)];

A= [eq238; eq235; eqPu9; eqX; eqY];

%Timestep
year= 365*24*3600;
timestep= 1*3600;
i= 1;
N_start= [N8_start; N5_start; N9_start; NX_start; NY_start];
N(:,1)= N_start;

%Forward Euler Method
for t= [timestep:timestep:(year-timestep)]
N(:, i+1)= (A*timestep*N(:,i) + N(:,i));
i= i+1; 
end

i= 1;
N_mat(:,1)= N_start;
%Matrix Exponential Method
for t= [timestep:timestep:(year-timestep)]
N_mat(:,i+1) = eye(size(A))*N_mat(:, i) + timestep.*A*N_mat(:, i) + (1/2).*(timestep.*A)^2*N_mat(:, i) + 1/6.*(timestep.*A)^3*N_mat(:, i);
i=i+1;
end

%% SECTION 4 - Analytical solutions
t= linspace(0, year, year/timestep);
N8_exact= @(t) N8_start*exp(-(sigma_c_U8 + sigma_f_U8)*Flux*t);
N5_exact= @(t) N5_start*exp(-(sigma_c_U5 + sigma_f_U5)*Flux*t);
K= (N8_start*sigma_c_U8)/((sigma_c_P9+sigma_f_P9)-(sigma_c_U8+sigma_f_U8));
NP9_exact= @(t) K*(exp(-(sigma_c_U8 + sigma_f_U8)*Flux*t) - exp(-(sigma_c_P9 + sigma_f_P9)*Flux*t));
K2= (Flux*FYX*sigma_f_U5*N5_start)/(Flux*sigma_c_X+log(2)/Thalf_X - (sigma_c_U5 +sigma_f_U5)*Flux);
NX_exact= @(t) K2*(exp(-(sigma_c_U5 + sigma_f_U5)*Flux*t) - exp(-(sigma_c_X*Flux + log(2)/Thalf_X)*t));
K3= (FYY*sigma_f_U5*N5_start)/(sigma_c_Y- (sigma_c_U5 +sigma_f_U5));
NY_exact= @(t) K3*(exp(-(sigma_c_U5 + sigma_f_U5)*Flux*t) - exp(-(sigma_c_Y*Flux)*t));
% prova= NX_exact(t);

%plot euler solution
figure(1)
plot(t, N(1,[1:8760]), 'linewidth', 2)
hold on
plot(t, N(2,[1:8760]), 'linewidth', 2)
plot(t, N(3,[1:8760]), 'linewidth', 2)
plot(t, N(4,[1:8760]), 'linewidth', 2)
plot(t, N(5,[1:8760]), 'linewidth', 2)
set(gca,'YScale', 'log')
grid on
title('Euler solution - Figure included in the report')
legend('Euler U-8', 'Euler U-5', 'Euler PU-239','Euler X', 'Euler Y', ...
    'Location','SouthEast', 'FontSize', 6)
%saveas(gcf,'euler_1year.png')
%saveas(gcf,'euler_1year','epsc')

%Plot of the analytical and numerical
figure(2)
plot(t, N8_exact(t), 'linewidth', 2)
hold on
plot(t, N5_exact(t), 'linewidth', 2)
plot(t, NP9_exact(t), 'linewidth', 2)
plot(t, NX_exact(t), 'linewidth', 2)
plot(t, NY_exact(t), 'linewidth', 2)
plot(t, N(1,[1:8760]), 'linewidth', 2)
plot(t, N(2,[1:8760]), 'linewidth', 2)
plot(t, N(3,[1:8760]), 'linewidth', 2)
plot(t, N(4,[1:8760]), 'linewidth', 2)
plot(t, N(5,[1:8760]), 'linewidth', 2)
plot(t, N_mat(1,[1:8760]), 'linewidth', 2)
plot(t, N_mat(2,[1:8760]), 'linewidth', 2)
plot(t, N_mat(3,[1:8760]), 'linewidth', 2)
plot(t, N_mat(4,[1:8760]), 'linewidth', 2)
plot(t, N_mat(5,[1:8760]), 'linewidth', 2)
set(gca,'YScale', 'log')
grid on
title('All the solutions - analytical and numericals')
legend('Exact U-8', 'Exact U-5', 'Exact Pu-9', 'Exact X','Exact Y', ....
    'Numerical U-8', 'Numerical U-5', 'Numerical PU-239','Numerical X', 'Numerical Y', ...
    'Mat Exp  U-8', 'Mat Exp U-5', 'Mat Exp Pu-239', 'Mat Exp X', 'Mat Exp Y', ...
    'Location','SouthEast', 'FontSize', 6)
%saveas(gcf,'Mat_and_num.png')
%saveas(gcf,'Mat_and_num','epsc')

%plot of the 1day focus
figure(3)
plot([0:(3600*24)], N8_exact([0:(3600*24)]), 'linewidth', 2)
hold on
plot([0:(3600*24)], N5_exact([0:(3600*24)]), 'linewidth', 2)
plot([0:(3600*24)], NP9_exact([0:(3600*24)]), 'linewidth', 2)
plot([0:(3600*24)], NX_exact([0:(3600*24)]), 'linewidth', 2)
plot([0:(3600*24)], NY_exact([0:(3600*24)]), 'linewidth', 2)
title('1-day  focus')
set(gca,'YScale', 'log')
grid on
legend('Exact U-8', 'Exact U-5', 'Exact Pu-9', 'Exact X','Exact Y', ....
    'Location','SouthEast', 'FontSize', 6)
%saveas(gcf,'1day.png')
%saveas(gcf,'1day','epsc')

%plot of the relative difference
t= [timestep:timestep:(year-timestep)];
t= [0,t];
err8=[];
err5=[];
err9=[];
errX=[];
errY=[];


err8= abs(N8_exact(t)-N(1,:))./N8_exact(t);
err5= abs(N5_exact(t)-N(2,:))./N5_exact(t);
err9= abs(NP9_exact(t)-N(3,:))./NP9_exact(t);
errX= abs(NX_exact(t)-N(4,:))./NX_exact(t);
errY= abs(NY_exact(t)-N(5,:))./NY_exact(t);


figure(4)
plot(t, err8, 'linewidth', 2)
hold on
plot(t, err5, 'linewidth', 2)
plot(t, err9, 'linewidth', 2)
plot(t, errX, 'linewidth', 2)
plot(t, errY, 'linewidth', 2)
set(gca,'YScale', 'log')
grid on
legend('Error U-8', 'Error U-5', 'Error Pu-9', 'Error X','Error Y', ...
    'Location','SouthEast', 'FontSize', 6)
title('Relative difference with respect to analytical solution')
%saveas(gcf,'relative_difference.png')
%saveas(gcf,'relative_difference','epsc')

%% SECTION 5 - Convergence properties of the two methods

err_euler= [];
err_mat= [];

for timestep=[60,(60*30), (60*60), (2*60*60), (4*60*60), (6*60*60)]
i=1;
Neuler_conv= [];
Nmat_conv= [];
Neuler_conv(:,1)= N_start;
Nmat_conv(:, 1)= N_start;

%Forward Euler Method for different time bins
for t= [timestep:timestep:(year-timestep)]
Neuler_conv(:, i+1)= (A*timestep*Neuler_conv(:,i) + Neuler_conv(:,i));
i= i+1; 
end

i=1;
%Matrix Exponential for different time bins
for t= [timestep:timestep:(year-timestep)]
Nmat_conv(:,i+1) = eye(size(A))*Nmat_conv(:, i) + timestep.*A*Nmat_conv(:, i) + (1/2).*(timestep.*A)^2*Nmat_conv(:, i) + 1/6.*(timestep.*A)^3*Nmat_conv(:, i);
i=i+1;
end

time=1*24*60*60;
index= 1+ time/timestep ; %not sure wheter it should be 1+ratio or just ratio

err_euler= [err_euler, abs(Neuler_conv(4,index)-NX_exact(time))];
err_mat= [err_mat, abs(Nmat_conv(4,index)-NX_exact(time))];

waitbar(timestep/21600)
end

time_bins= [60,(60*30), (60*60), (2*60*60), (4*60*60), (6*60*60)];

%plot of the error of the Euler Method - LOG SCALE to show order 1
figure(5)
loglog(time_bins,err_euler, 'linewidth', 2);
hold on
loglog(time_bins,err_euler, 'o', 'markersize', 5, 'color', 'b');
loglog(time_bins, time_bins, 'linewidth', 2, 'LineStyle', '--');
title('Error Euler Method')
legend('error', 'error', 'order 1', 'Location', 'SouthEast')
grid on
%saveas(gcf,'err_euler.png')
%saveas(gcf,'err_euler','epsc')

% plot of the error of the matrix exponential method - LOG SCALE
figure(6)
loglog(time_bins,err_mat, 'linewidth', 2);
hold on
loglog(time_bins,err_mat, 'o', 'markersize', 5, 'color', 'b');
loglog(time_bins, time_bins.^3, 'linewidth', 2, 'LineStyle', '--');
title('Error Matrix Exponential Method')
legend('error', 'error', 'order 3', 'Location', 'SouthEast')
grid on
%saveas(gcf,'err_mat_exp.png')
%saveas(gcf,'err_mat_exp','epsc')

%% SECTION 6 - CONSTANT POWER PART

%initial situation for the Flux
En_Fiss_Rate= (kappa.*sigmaf)*N_start;
Power= 937.5;
Flux= Power/En_Fiss_Rate;
Flux_vec= [Flux];
timestep= 60;

% BUILDING UP THE EQUATIONS
eq238_new= [-(sigma_c_U8+sigma_f_U8)*Flux, 0, 0, 0, 0];
eq235_new= [0,-(sigma_c_U5 + sigma_f_U5)*Flux, 0, 0, 0];
eqPu9_new= [Flux*sigma_c_U8, 0, -Flux*(sigma_c_P9+sigma_f_P9) , 0, 0];
eqX_new= [Flux*FYX*sigma_f_U8, Flux*FYX*sigma_f_U5, Flux*FYX*sigma_f_P9, (-sigma_c_X*Flux - log(2)/Thalf_X), 0]; 
eqY_new= [Flux*FYY*sigma_f_U8, Flux*FYY*sigma_f_U5, Flux*FYY*sigma_f_P9,0,(-sigma_c_Y*Flux)];
A_new= [eq238_new; eq235_new; eqPu9_new; eqX_new; eqY_new];


% Neuler_new(:, 1)= N_start;
Nmat_new(:, 1)= N_start;
i=1;

for t= [timestep:timestep:(year-timestep)]

    %select EULER
% Neuler_new(:, i+1)= (A_new*timestep*Neuler_new(:,i) + Neuler_new(:,i));

    %select MATRIX EXPONENTIAL
Nmat_new(:,i+1) = eye(size(A_new))*Nmat_new(:, i) + timestep.*A_new*Nmat_new(:, i) + (1/2).*(timestep.*A_new)^2*Nmat_new(:, i) + 1/6.*(timestep.*A_new)^3*Nmat_new(:, i);

% Update the Flux

    %select EULER
% En_Fiss_Rate= (kappa.*sigmaf)*Neuler_new(:,i+1);

    %select MATRIX EXPONENTIAL
En_Fiss_Rate= (kappa.*sigmaf)*Nmat_new(:,i+1);

Flux= Power/En_Fiss_Rate;
Flux_vec= [Flux_vec, Flux];

% Recompute the Matrix (shall we consider contribution of Pu and U?)
eq238_new= [-(sigma_c_U8+sigma_f_U8)*Flux, 0, 0, 0, 0];
eq235_new= [0,-(sigma_c_U5 + sigma_f_U5)*Flux, 0, 0, 0];
eqPu9_new= [Flux*sigma_c_U8, 0, -Flux*(sigma_c_P9+sigma_f_P9) , 0, 0];
eqX_new= [Flux*FYX*sigma_f_U8, Flux*FYX*sigma_f_U5, Flux*FYX*sigma_f_P9, (-sigma_c_X*Flux - log(2)/Thalf_X), 0]; 
eqY_new= [Flux*FYY*sigma_f_U8, Flux*FYY*sigma_f_U5, Flux*FYY*sigma_f_P9,0,(-sigma_c_Y*Flux)];
A_new= [eq238_new; eq235_new; eqPu9_new; eqX_new; eqY_new];

%move to next timestep
i= i+1; 

waitbar(i*timestep/(year-timestep))
end

% i=1;
% %Matrix Exponential for different time bins
% for t= [timestep:timestep:(year-timestep)]
% Nmat_conv(:,i+1) = (eye(size(A)) + timestep.*A + (1/2).*(timestep.*A)^2 + 1/6.*(timestep.*A)^3)*Nmat_conv(:, i);
% i=i+1;
% end


%% plotting of nuclide inventory in this case
t= linspace(0, year, year/timestep);
figure(7)

% PLOT EULER METHOD
% plot(t, Neuler_new(1,:), 'linewidth', 2)
% hold on
% plot(t, Neuler_new(2,:), 'linewidth', 2)
% plot(t, Neuler_new(3,:), 'linewidth', 2)
% plot(t, Neuler_new(4,:), 'linewidth', 2)
% plot(t, Neuler_new(5,:), 'linewidth', 2)

%PLOT MATRIX EXPONENTIAL
plot(t, Nmat_new(1,:), 'linewidth', 2)
hold on
plot(t, Nmat_new(2,:), 'linewidth', 2)
plot(t, Nmat_new(3,:), 'linewidth', 2)
plot(t, Nmat_new(4,:), 'linewidth', 2)
plot(t, Nmat_new(5,:), 'linewidth', 2)

title('Nuclide inventory when constant power')
set(gca,'YScale', 'log')
grid on
legend( 'Numerical U-8', 'Numerical U-5', 'Numerical PU-239','Numerical X', 'Numerical Y', ...
    'Location','SouthEast', 'FontSize', 6)
%saveas(gcf,'nuclide_const_pow.png')
%saveas(gcf,'nuclide_const_pow','epsc')

%% flux evolution
figure(8)
plot(t, Flux_vec, 'linewidth', 2)
title('Flux over time for constant power scenario')
grid on
%saveas(gcf,'flux_const_pow.png')
%saveas(gcf,'flux_const_pow','epsc')

