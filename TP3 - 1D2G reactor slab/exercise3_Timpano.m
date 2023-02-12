close all 
clear all 
clc

%% EXERCISE 3 - Slab Reactor 1D-2G

% ### Algorithm description ###

% For the solution of the problem the following procedure was employed:

% ### PREPROCESSING - ANALYTICAL SOLUTION ###

% SECTION 1 - data from the TP scriptum were reported;
% SECTION 2 - analytical solution was obtained:
% - the set of diffusion equations was written in matricial form
% - the solution for flux and keff was obtained using Matlab eigs 
% - the solution was normalised at its maximum value and written as
% function handle

% ### NUMERICAL SOLVER FOR BARE SYSTEM ###

% SECTION 3 - building of the numerical solver:
% - diffusion properties were defined in vectors to link them to mesh
%   position, the same approach is adopeted in section 5 considering the
%   reflector
% - diffusion properties were used to explicitly compute the beta
%   coefficients using for cycles on the given core mesh points 
% - diagonals of the M1 and M2 matrices are built using for cycles
% - submatrices for group 1 and group 2 were built using diag properties
% - boundary conditions were set accordingly to information on the symmetry
%   of the system and on zero incoming current at the boundary for both
%   groups


% ADDITIONAL NOTE ON MATRIX BUILDING
% definition of block matrices referring to group properties, 
% minding to consider the scattering source inside the M matrix (negative)
% so to leave a purely fissional term in the F matrix. Please note that this 
% was not the only possible approach and that alternatively we could have 
% included \textit{in-group scattering} on both sides, though lacking 
% further Matrix simplification and maybe affecting the number of iterations.

% - the two matrices were concatenated building matrix M
% - the downscattering term was added to M matrix, summing a diagonal matrix
%   whose diagonal index is nc [ i.e. diag(down_scatter, -nc)]
% SECTION 4 - direct solver (eigs) was used to solve the built problem
% - the solution was normalised according to the maximum value of the flux
%   in the fast group

% ### NUMERICAL SOLVER FOR REFLECTED SYSTEM ###
% SECTION 5 - repeat the previous steps for reflective system 
% SECTION 6 - direct solver (eigs) was used to solve the built problem

% ### ADDITIONAL REQUESTED OUTPUTS ###
% SECTION 7 - Requested outputs
% - function lsqcurvefit was employed to find the buckling 
% - for this sake, the numerical solution was fitted with a cosine function
% - absolute and relative error are obtained and plotted
% - finite difference approximation is used to get info on the current at
%   the boundary of the reflected system

%% Info on output
% All the output won't be shown in Command Window apart from those required
% in TP3 Slides + graphs 

%% SECTION 1 - Data

sigma_aC= [0.002, 0.060]; %cm-1
vsigma_fC= [0.001, 0.069]; %cm-1
sigma_s12C= 0.038;
sigma_aR= [0.00, 0.012]; %cm-1
sigma_s12R= 0.040;
sigma_tC= [sigma_aC(1)+ sigma_s12C, sigma_aC(2)+sigma_s12C];
sigma_tR= [sigma_aR(1)+ sigma_s12R, sigma_aR(2)+sigma_s12R];

a= 50; %cm
b= 10; %cm

Dc= [1.130, 0.160];
Dr= [1.130, 0.160];
dc= 2*Dc(1);
dr= 2*Dr(1);

epsilon=10^(-7);

%% SECTION 2 - Analytical solution
% dc= 2*Dc;
% dr= 2*Dr;

%bare reactor
B= pi/(2*a); %geometrical buckling excluding extrapolation length
% Lc1= (Dc/sigma_aC(1))^0.5;
x= linspace(0,a,1000);
A= [Dc(1)*B^2 + sigma_aC(1) + sigma_s12C, 0; -sigma_s12C, Dc(2)*B^2 + sigma_aC(2)]; %matricial form of the set of diffusion equations
F= [vsigma_fC(1), vsigma_fC(2); 0, 0]; %source term for the two equations
input0= A\F;
[PHI_a, k_a]= eigs(input0); %solving the problem using built-in Matlab properties
PHI1_K1= abs(PHI_a(1,1)); %extract maximum value of fast flux
PHI2_K1= abs(PHI_a(2,1)); %extract maximum value of thermal flux

%plot of the analytical solution
sol_gruppo1= @(x)(PHI1_K1/PHI1_K1)*cos(B*x); %fast flux as function handle
sol_gruppo2= @(x) (PHI2_K1/PHI1_K1)*cos(B*x); %thermal flux as function handle

figure(1)
plot(x, sol_gruppo1(x), x, sol_gruppo2(x), 'linewidth', 2);
title('Analytical solution for bare reactor');
legend('analytical solution group 1', 'analytical solution group 2', 'Location', 'NorthEast') 
grid on
% saveas(gcf,'an_groups.png')
% saveas(gcf,'an_groups','epsc')

%% SECTION 3 - Numerical solver for bare system
%% Dimension of the matrix
h= 0.1; %mesh size
nc= a/h; %dimension of the core
nr= b/h; %dimension of the reflector
nt= (a+b)/h; %total dimension

%% BARE SYSTEM
%% Definition of the beta coefficients
D_1=Dc(1)*ones(nc,1); %vectors bearing diffusion properties for mesh points in the core group 1
D_2=Dc(2)*ones(nc,1); %vectors bearing diffusion properties for mesh points in the core group 2

%explicit construction of the beta coefficients according to diffusion properties
%explicit construction is repeated for group 1 and group 2
for i= 1:(nc-1)
    beta_plus1(i)= (2*D_1(i)*D_1(i+1))/(D_1(i+1)+D_1(i));
    beta_plus2(i)= (2*D_2(i)*D_2(i+1))/(D_2(i+1)+D_2(i));
end

for i=2:nc
    beta_minus1(i)= (2*D_1(i)*D_1(i-1))/(D_1(i-1)+D_1(i));
    beta_minus2(i)= (2*D_2(i)*D_2(i-1))/(D_2(i-1)+D_2(i));
end 
beta_minus1(1)=[];
beta_minus2(1)=[];

%construction of the M matrix by construction of group 1 and 2 submatrices
%GROUP 1
main_diag1= zeros(nc,1);
up_diag1= zeros(nc-1,1);
low_diag1= zeros(nc-1,1);

%GROUP 2
main_diag2= zeros(nc,1);
up_diag2= zeros(nc-1,1);
low_diag2= zeros(nc-1,1);


%MATRIX GROUP 1
for i=1:(nc-1)
    main_diag1(i)= (beta_plus1(i)+ beta_minus1(i))/h^2 + sigma_tC(1);
    up_diag1(i)= -beta_plus1(i)/h^2;
    low_diag1(i)= -beta_minus1(i)/h^2;
end

%setting the boundary conditions
BC1= 0.5*(1/(h/(4*Dc(1)) + 1));

M1= diag(main_diag1,0) + diag(up_diag1,1) + diag(low_diag1, -1);
M1(1,1)= beta_plus1(1)/h^2 + sigma_tC(1);
M1(nc,nc)= (BC1*h + beta_minus1(end))/h^2 + sigma_tC(1);

%MATRIX GROUP 2
for i=1:(nc-1)
    main_diag2(i)= (beta_plus2(i)+ beta_minus2(i))/h^2 + sigma_aC(2);
    up_diag2(i)= -beta_plus2(i)/h^2;
    low_diag2(i)= -beta_minus2(i)/h^2;
end

%setting the boundary conditions
BC2= 0.5*(1/(h/(4*Dc(2)) + 1));

M2= diag(main_diag2,0) + diag(up_diag2,1) + diag(low_diag2, -1);
M2(1,1)= beta_plus2(1)/h^2 + sigma_aC(2);
M2(nc,nc)= (BC2*h + beta_minus2(end))/h^2 + sigma_aC(2);

%There is a need to concatenate the two matrices
%since the flux vector will now accomodate the two groups
M= zeros(2*nc,2*nc);
M= [M1, zeros(nc,nc); zeros(nc,nc), M2];

%we need to account for the downscattering term from group1 to group2
down_scatter= -sigma_s12C*ones(nc,1);
M= M + diag(down_scatter, -nc); 

%% Costruction of the F matrix 
ff= [(vsigma_fC(1))*ones(nc,1); zeros(nc,1)];
tf= (vsigma_fC(2))*ones(nc,1);
F= zeros(2*nc,2*nc);
F= F + diag(ff,0) + diag(tf,nc);

%% SECTION 4 - Solution (NO POWER ITERATION TECHNIQUE)

%solution 
input= M\F;
[PHI, KEFF_anal]= eigs(input);
keff_bare= KEFF_anal(1,1);

%we normalise accordingly to the maximal flux in the fast group
pos= [0.05:0.1:49.95]; %build position vector for plotting purposes
normalsol_fund1= PHI((1:500),1)/max(PHI((1:500),1));
normalsol_fund2= PHI((501:1000),1)/max(PHI((1:500),1));

%plot the flux as numerical solution FUNDAMENTAL MODE
figure(2)
plot(pos,normalsol_fund1, pos, normalsol_fund2, 'linewidth', 2)
title('Numerical solution for bare system')
legend('Flux Numerical Group 1', 'Flux Numerical Group 2','Location','SouthWest')
grid on
% saveas(gcf,'num_groups.png')
%saveas(gcf,'num_groups','epsc')

%% SECTION 5 - Numerical solver for REFLECTIVE SYSTEM 

%% Definition of the beta coefficients
% build the diffusion properties concatenating core and reflector
% properties, hence linking diffusion property to mesh position
D_1t=[Dc(1)*ones(nc,1); Dr(1)*ones(nr,1)];
D_2t=[Dc(2)*ones(nc,1); Dr(2)*ones(nr,1)];

%build the total and absorption cross-section as a function of mesh
%position
sigma_t1=[sigma_tC(1)*ones(nc,1); sigma_tR(1)*ones(nr,1)];
sigma_t2=[sigma_tC(2)*ones(nc,1); sigma_tR(2)*ones(nr,1)];
sigma_a1=[sigma_aC(1)*ones(nc,1); sigma_aR(1)*ones(nr,1)];
sigma_a2=[sigma_aC(2)*ones(nc,1); sigma_aR(2)*ones(nr,1)];

%explicit construction of the beta coefficients
for i= 1:(nt-1)
    beta_plus1t(i)= (2*D_1t(i)*D_1t(i+1))/(D_1t(i+1)+D_1t(i));
    beta_plus2t(i)= (2*D_2t(i)*D_2t(i+1))/(D_2t(i+1)+D_2t(i));
end

for i=2:nt
    beta_minus1t(i)= (2*D_1t(i)*D_1t(i-1))/(D_1t(i-1)+D_1t(i));
    beta_minus2t(i)= (2*D_2t(i)*D_2t(i-1))/(D_2t(i-1)+D_2t(i));
end 
beta_minus1t(1)=[];
beta_minus2t(1)=[];

%construction of the M matrix by construction of group 1 and 2 submatrices
%GROUP 1
main_diag1t= zeros(nt,1);
up_diag1t= zeros(nt-1,1);
low_diag1t= zeros(nt-1,1);

%GROUP 2
main_diag2t= zeros(nt,1);
up_diag2t= zeros(nt-1,1);
low_diag2t= zeros(nt-1,1);


%MATRIX GROUP 1
for i=1:(nt-1)
    main_diag1t(i)= (beta_plus1t(i)+ beta_minus1t(i))/h^2 + sigma_t1(i); 
    up_diag1t(i)= -beta_plus1t(i)/h^2;
    low_diag1t(i)= -beta_minus1t(i)/h^2;
end

%setting the boundary conditions
BC1t= 0.5*(1/(h/(4*Dc(1)) + 1));

M1t= diag(main_diag1t,0) + diag(up_diag1t,1) + diag(low_diag1t, -1);
M1t(1,1)= beta_plus1t(1)/h^2 + sigma_t1(1);
M1t(nt,nt)= (BC1t*h + beta_minus1t(end))/h^2 + sigma_t1(end);

%MATRIX GROUP 2
for i=1:(nt-1)
    main_diag2t(i)= (beta_plus2t(i)+ beta_minus2t(i))/h^2 + sigma_a2(i);
    up_diag2t(i)= -beta_plus2t(i)/h^2;
    low_diag2t(i)= -beta_minus2t(i)/h^2;
end

%setting the boundary conditions
BC2t= 0.5*(1/(h/(4*Dr(2)) + 1));

M2t= diag(main_diag2t,0) + diag(up_diag2t,1) + diag(low_diag2t, -1);
M2t(1,1)= beta_plus2t(1)/h^2 + sigma_a2(1);
M2t(nt,nt)= (BC2t*h + beta_minus2t(end))/h^2 + sigma_a2(end);

%There is a need to concatenate the two matrices
%since the flux vector will now accomodate the two groups
Mt= zeros(2*nt,2*nt);
Mt= [M1t, zeros(nt,nt); zeros(nt,nt), M2t];

%we need to account for the downscattering term from group1 to group2
down_scatter_t= [-sigma_s12C*ones(nc,1); -sigma_s12R*ones(nr,1)] ; %we need to add the part referring to the reflector
Mt= Mt + diag(down_scatter_t, -nt); 

%% Costruction of the F matrix 
ff2= [(vsigma_fC(1))*ones(nc,1); zeros(nr,1); zeros(nc,1); zeros(nr,1)];
tf2= [(vsigma_fC(2))*ones(nc,1); zeros(nr,1)];
Ft= zeros(2*nt,2*nt);
Ft= Ft + diag(ff2,0) + diag(tf2,nt);

%% SECTION 6 - Solution for reflective system (NO POWER ITERATION)

%analytical solution for comparison keff
input2= Mt\Ft;
[PHIt, KEFF_anal_t]= eigs(input2);
keff_t= KEFF_anal_t(1,1);

%we normalise accordingly to the maximal flux in the fast group
pos_t= [0.05:0.1:59.95]; %build position vector for plotting purposes
normalsol_fund1t= PHIt((1:600),1)/max(PHIt((1:600),1));
normalsol_fund2t= PHIt((601:1200),1)/max(PHIt((1:600),1));
normalsol_second1t= PHIt((1:600),2)/max((abs(PHIt((1:600),1))));
normalsol_second2t= PHIt((601:1200),2)/(max(abs(PHIt((1:600),1))));

%plot the flux as numerical solution FUNDAMENTAL MODE (1st armonics)
figure(3)   
plot(pos_t,normalsol_fund1t, pos_t, normalsol_fund2t, 'linewidth', 2)
title('Numerical solution for reflected system - First Armonics')
legend('Flux Numerical Group 1', 'Flux Numerical Group 2','Location','NorthEast')
grid on
%saveas(gcf,'num_groups_refl.png')
%saveas(gcf,'num_groups_refl','epsc')

%plot the flux as numerical solution (2nd armonic)
figure(4)
plot(pos_t,normalsol_second1t, pos_t, normalsol_second2t, 'linewidth', 2)
title('Numerical solution for reflected system - Second Armonics')
legend('Flux Numerical Group 1', 'Flux Numerical Group 2','Location','NorthEast')
grid on
% saveas(gcf,'num_groups_refl_2.png')
% saveas(gcf,'num_groups_refl_2','epsc')

%plot the flux (TOGETHER 1st and 2nd armonic)
figure(5)
plot(pos_t,normalsol_fund1t, pos_t, normalsol_fund2t, pos_t,normalsol_second1t, pos_t, normalsol_second2t, 'linewidth', 2)
title('Numerical solution for reflected system - Both armonics both groups')
legend('Flux Numerical Group 1 - Fundamental', 'Flux Numerical Group 2 - Fundamental', 'Flux Numerical Group 1 - Second', 'Flux Numerical Group 2 - Second','Location','SouthEast')
grid on
%saveas(gcf,'armonics_def.png')
%saveas(gcf,'armonics_def','epsc')

%% SECTION 7 - Requested output 
%% Cosine fitting to bind the buckling
Bfit_fast= 2*pi/100; %first estimate
Bfit_thermal= 2*pi/100; %first estimate
A1= 1;
A2= 1;

p0=[Bfit_fast, A1]; %vector of coefficients input to lsqcurvefit
p1=[Bfit_thermal, A2]; %vector of coefficients input to lsqcurvefit

FUN_fast= @(p0,pos) p0(2)*cos(p0(1)*pos);
Buck_fast = lsqcurvefit(FUN_fast,p0,pos, normalsol_fund1'); %cosine function is fitted into the numerical solution

FUN_th= @(p0,pos) p0(2)*cos(p0(1)*pos);
Buck_th = lsqcurvefit(FUN_th,p1,pos, normalsol_fund2'); %cosine function is fitted into the numerical solution

%Fitted function
fit_fast= @(x) Buck_fast(2)*cos(Buck_fast(1)*x);
fit_thermal= @(x) Buck_th(2)*cos(Buck_th(1)*x);

x= linspace(0,50, 1000);
%Plot of the fit
figure(6)
plot(pos, normalsol_fund1,'o', x, fit_fast(x),'linewidth', 2);
hold on
plot(pos, normalsol_fund2,'o', x, fit_thermal(x),'linewidth', 2);
legend('numerical points fast', 'fit fast', 'numerical points th', 'fit thermal');
title('Cosine fitting of the numerical solution')
grid on
% saveas(gcf,'cosine_fitting.png')
% saveas(gcf,'cosine_fitting','epsc')

%% error on the keff
err_keff= abs(k_a(1,1)-keff_bare)*10^5;

%% error in percentage between analytical and numerical solution FLUX

sol_gruppo1= @(x)(PHI1_K1/PHI1_K1)*cos(B*x);
sol_gruppo2= @(x) (PHI2_K1/PHI1_K1)*cos(B*x);
err_fast= abs((normalsol_fund1 - sol_gruppo1(pos)')./(sol_gruppo1(pos)'))*100;
err_th= abs((normalsol_fund2 - sol_gruppo2(pos)')./(sol_gruppo2(pos)'))*100;

figure(7)
plot(pos([1:460]), err_fast([1:460]), pos([1:460]), err_th([1:460]),'linewidth', 2);
legend('error fast', 'error thermal', 'Location', 'SouthEast');
title('Error in percentage')
grid on
% saveas(gcf,'error_closeup.png')
% saveas(gcf,'error_closeup','epsc')

%% absolute error between analytical and numerical solution FLUX

sol_gruppo1= @(x)(PHI1_K1/PHI1_K1)*cos(B*x);
sol_gruppo2= @(x) (PHI2_K1/PHI1_K1)*cos(B*x);
err_fast= abs(normalsol_fund1 - sol_gruppo1(pos)');
err_th= abs(normalsol_fund2 - sol_gruppo2(pos)');

figure(8)
plot(pos, err_fast, pos, err_th,'linewidth', 2);
legend('error fast', 'error thermal', 'Location', 'SouthEast');
title('Absolute error')
grid on
%saveas(gcf,'abs_error.png')
%saveas(gcf,'abs_error','epsc')

%% current at the boundary using centered finite difference

%compute the value of the current at the boundary using finite difference
J1_refl= -(Dr(1)*normalsol_fund1t(501) - Dc(1)*normalsol_fund1t(500))/(h);
J2_refl= -(Dr(2)*normalsol_fund2t(501) - Dc(2)*normalsol_fund2t(500))/(h);

%% PRINT THE REQUESTED OUTPUT
format long

disp('REQUESTED RESULTS FOR BARE REACTOR')
disp('Keff analytical of the bare system:')
k_a(1,1)
disp('Keff numerical of the bare system:')
keff_bare

disp('Buckling of the fast neutron flux:')
Buck_fast(1)
disp('Buckling of the thermal neutron flux:')
Buck_th(1)
disp('Error on keff:')
err_keff

disp('REQUESTED RESULTS FOR REFLECTED SYSTEM')
disp('Keff numerical of the reflected system fundamental:')
keff_t
disp('Keff numerical of the reflected system second harmonic:')
KEFF_anal_t(2,2)
disp('Current at the boundary of reflected system FAST:')
J1_refl
disp('Current at the boundary of reflected system THERMAL:')
J2_refl
