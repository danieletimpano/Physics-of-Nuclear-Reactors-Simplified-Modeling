close all 
clear all 
clc

%% EXERCISE 2 - Slab Reactor 1D-1G

% ### Algorithm description ###

% For the solution of the problem the following procedure was employed:

% ### PREPROCESSING - ANALYTICAL SOLUTION ###

% SECTION 1 - data from the TP scriptum were reported;
% SECTION 2 - analytical solution was written as a function handle and analytical solution was plotted;

% ### NUMERICAL SOLVER FOR BARE SYSTEM ###

% SECTION 3 - building of the numerical solver:
% - use of diag command in order to build F matrix
% - construction of the M matrix according to what performed in TP1 (first version),
% please note that the beta coefficients were not explicited for the sake of
% brevity as the uniformity of the materials employed did not required
% that. Nonetheless their explicit form is included in the reports and used
% when needed in TP3.
% SECTION 4 - power iteration technique
% - set the first guess for flux and keff
% - set the source term according to the guess
% - solve the problem for the flux
% - use the flux to update keff accordingly to TP2 info
% - update the source term
% - compute the relative difference between keff at each iteration to set
%   the stop criteria
% SECTION 5 - OUPUT AND COMPARISON
% - compute the analytical solution for keff
% - use of current continuity to find flux and current at the boundary
% - use of forward finite difference method to compute the current
% - compute the error for the flux and for the current along the mesh

% ### NUMERICAL SOLVER FOR REFLECTED SYSTEM ###

% SECTION 6 - Building of the matrix of the complete system using block
% matrix construction. Block matrix approach was employed since there is no
% disomogenity in the properties of the two materials a part from fission
% properties. Beta coefficients were used in TP3. 
% SECTION 7 - Repeat Power iteration for the reflected system
% SECTION 8 - Print the desired output as requested by

%% Info on Output
% All the output won't be shown in Command Window apart from those required
% in TP2 Slides + graphs 
%% SECTION 1 - Data

sigma_aR= 0.01; %cm-1
sigma_sR= 0.3; %cm-1
sigma_aC= 0.01; %cm-1
sigma_sC=0.3; %cm-1
sigma_tC= sigma_sC+ sigma_aC;
sigma_tR= sigma_sR+ sigma_aR;
vsigma_fC=0.015; %cm-1
a= 20; %cm
b= 10; %cm

Dc= 1/(3*(sigma_sC + sigma_aC)); %cm
Dr= 1/(3*(sigma_sR + sigma_aR)); %cm
dc= 2*Dc; %cm
dr= 2*Dr; %cm

epsilon=10^(-7); % convergence criteria

%% SECTION 2 - UNREFLECTED SYSTEM ANALYTICAL SOLUTION
%analytical solution for the flux
%we normalise to the value at 0.05 in order to compare the analytical and
%the numerical solution

B= pi/(2*a+2*dc); %definition of the buckling
Lc= (Dc/sigma_aC)^0.5; %definition of the diffusion length in the core
x= linspace(0,(a+dc),1000);
sol= @(x) cos(B*x)/cos(B*0.05); %analytical solution normalised
kinf_sol= Lc^2*B^2 + 1; %kinf according to criticality condition 

%analytical solution for the current
sol_curr= @(x) Dc*B*sin(B*x)/cos(B*0.05); %application of Fick's Law in order to have the current

%plot the analytical solution
figure(1)
plot(x, sol(x),x, sol_curr(x),'linewidth', 2)
legend('flux', 'current')
title('Analytical solution')
grid on
%saveas(gcf,'anal_both.png')
%saveas(gcf,'anal_both','epsc')

%% SECTION 3 - NUMERICAL SOLVER FOR BARE SYSTEM

%% choice for the mesh
%size of the mesh
n= 200;

%% Construction of the F matrix
f= (vsigma_fC)*ones(n,1);
F= -diag(f,0); 

%% Construction of the M matrix
h= a/n; %definition of the mesh size for unreflected system

%in the following, the same structure of "exercise1_Timpano" was repeated,
%with no explicit use of the beta coefficients. Note that the explicit form
%of the beta coefficients was still included in TP1 report and was employed
%when the code required it for TP3.
main_diag= (-2*Dc/h^2 - sigma_aC)*ones(n,1);
up_diag= Dc/h^2*ones(n-1,1);
low_diag= Dc/h^2*ones(n-1,1);

%costruction of the matrix
M= diag(main_diag,0) + diag(up_diag,1) + diag(low_diag, -1);

%boundary conditions at x=a
M(n,n)= (2*Dc)/h^2*(1/(h/(4*Dc) + 1))-(3*Dc/h^2) - sigma_aC;

%boundary condition at x=0
M(1,1)= -Dc/h^2 - sigma_aC;

%% SECTION 4 - Power iteration technique

%first guess - build empty vector solutions
phi=[];
k=[];
S=zeros(n,n);

%set the first guess
diff=1; %set the difference so to enter in the while cycle
it=1; %save the number of iterations updating it
k(1)=1; %first guess for k
phi(:,1)=ones(n,1); %take the vector of fluxes=1 for the first try
S(:,1)= 1/k(1)*F*phi(:,1); %set the source term accordingly to first try flux

while (diff > epsilon)
    phi(:,it+1)= M\S(:,it);
    k(it+1)= k(it)*sum(phi(:,it+1))/sum(phi(:,it)); %compute the new keff accordingly to the suggested update
    S(:,it+1)= 1/k(it+1)*F*phi(:,it+1); %update the source term accordingly to the new keff
    diff= abs((k(it+1)-k(it))/(k(it+1))); %convergence criteria shall be set accordingly to relative difference between keff
    it= it+1;
end 

%% SECTION 5 - Output and comparison 
%normalise the solution according to maximum flux value
pos= [0.05:0.1:19.95];
normalsol= (phi(:,end)/max(phi(:,end)));

%compute analytical solution for comparison keff
input= M\F;
[PHI, KEFF_anal]= eig(input);
keff_anal= KEFF_anal(1,1);

%compute the value of the flux at the boundary
phi_n_a= (1/(0.25 + Dc/h))*(Dc/h)*normalsol(end);
%compute the value of the current at the boundary
J_na_dx= -Dc*(phi_n_a - normalsol(end))/(0.5*h); 

%error on the current at the boundary
err_curr_b= abs(sol_curr(a)-J_na_dx);

%Compute the current everywhere with FORWARD FINITE DIFFERENCE
% %finite difference method for current
for i=[1:199]
curr(i)= -Dc*(normalsol(i+1) - normalsol(i))/(h);
% curr(i)= -Dc*(phi(i+1,end) - phi(i,end))/(h);
end 
poscurr= pos([1:(end-1)]);

%CENTERED FINITE DIFFERENCE
% %finite difference method for current
% for i=[1:198]
% curr(i)= -Dc*(phi(i+2,end) - phi(i,end))/(2*h);
% end 
% poscurr= pos([2:(end-1)]);

%plot of the numerical solution flux and current
figure(2)
plot(pos,normalsol,poscurr, curr, 'linewidth', 2)
title('Numerical solution for bare system')
legend('flux numerical', 'current numerical')
grid on
%saveas(gcf,'num_both.png')
%saveas(gcf,'num_both','epsc')

%plot of the analytical and numerical solutions together for the FLUX
figure(3)
plot(x, sol(x),pos,normalsol, 'linewidth', 2)
legend('analytical', 'numerical')
title('Numerical and Analytical solution for bare system FLUX')
grid on
%saveas(gcf,'num_an_flux.png')
%saveas(gcf,'num_an_flux','epsc')

%plot of the analytical and numerical solution together for the CURRENT
figure(4)
plot(x, sol_curr(x),poscurr,curr, 'linewidth', 2)
legend('analytical', 'numerical', 'Location', 'SouthEast')
title('Numerical and Analytical solution for bare system CURRENT')
grid on
%saveas(gcf,'num_an_current.png')
%saveas(gcf,'num_an_current','epsc')

%error flux
err= abs(sol(pos)'-normalsol);
figure(5)
plot(pos,err, 'linewidth', 2, 'color', 'g')
title('error along the mesh for flux')
grid on
%saveas(gcf,'err_flux.png')
%saveas(gcf,'err_flux','epsc')

%error current
errcurr= abs(sol_curr(poscurr)-curr);
figure(6)
plot(poscurr,errcurr, 'linewidth', 2, 'color', 'g')
title('error along the mesh for current')
grid on
%saveas(gcf,'err_current.png')
%saveas(gcf,'err_current','epsc')


%% SECTION 6 - System with the reflector
%construction of the M2 matrix
%% Construction of the M matrix
n1= 200; %dimension of the core
n2= 100; %dimension of the reflector
ntot= n1+n2; %total dimension of the system
h= (a+b)/ntot;

M2=zeros(ntot, ntot);

%block 1 matrix (core)
main_diag= (-2*Dc/h^2 - sigma_aC)*ones(n1,1);
up_diag= Dc/h^2*ones(n1-1,1);
low_diag= Dc/h^2*ones(n1-1,1);

%block 2 matrix (reflector)
main_diag2= (-2*Dr/h^2 - sigma_aR)*ones(n2,1);
up_diag2= Dr/h^2*ones(n2,1);
low_diag2= Dr/h^2*ones(n2,1);

%total matrix (core + reflector)
main_diagtot= [main_diag; main_diag2];
up_diagtot= [up_diag; up_diag2];
low_diagtot= [low_diag; low_diag2];

%costruction of the matrix
M2= diag(main_diagtot,0) + diag(up_diagtot,1) + diag(low_diagtot, -1);

%boundary conditions at x=a
M2(ntot,ntot)= (2*Dr)/h^2*(1/(h/(4*Dr) + 1))-(3*Dr/h^2) - sigma_aR;

%boundary condition at x=0
M2(1,1)= -Dc/h^2 - sigma_aC;

f2= [(vsigma_fC)*ones(n1,1); zeros(n2,1)] ;
F2= -diag(f2,0);

%% SECTION 7 - Power iteration technique for reflector system

%first guess 
phi2=[];
k2=[];
S2=zeros(ntot,ntot);

%set the first guess
diff2=1;
it2=1;
k2(1)=1;
phi2(:,1)=ones(ntot,1);
S2(:,1)= 1/k2(1)*F2*phi2(:,1);

while (diff2 > epsilon)
    phi2(:,it2+1)= M2\S2(:,it2);
    k2(it2+1)= k2(it2)*sum(phi2(:,it2+1))/sum(phi2(:,it2));
    S2(:,it2+1)= 1/k2(it2+1)*F2*phi2(:,it2+1);
    diff2= abs((k2(it2+1)-k2(it2))/(k2(it2+1)));
    it2= it2+1;
end 

pos2= [0.05:0.1:29.95];
x2= linspace(0,(a+b),1000);


%plot of the numerical solution
figure(7)
plot(pos2,(phi2(:,end)/max(phi2(:,end))), 'linewidth', 2)
xline(20, '--', 'linewidth', 2)
legend('flux with reflector', 'boundary core-reflector', 'Location', 'SouthWest')
title('Numerical solution for reflected system')
grid on
%saveas(gcf,'flux_refl.png')
%saveas(gcf,'flux_refl','epsc')

normalsol2= (phi2(:,end)/max(phi2(:,end)));


%compute the value of the current at the boundary using finite difference
J_na_dx_refl= -(Dr*normalsol2(201) - Dc*normalsol2(200))/(h);


%% SECTION 8 - Print all the desired output

disp('REQUESTED OUTPUT')
format long

disp('Keff for bare reactor:')
keff_bare= k(end)

disp('Keff bare analytical:')
k_analytical=keff_anal

disp('error on keff')
err_keff=abs(keff_bare-keff_anal)

disp('Current at the core boundary for bare core NUMERICAL')
current_bound= J_na_dx

disp('Current at the core boundary for bare core ANALYTICAL')
current_analytical_bound= sol_curr(a)

disp('Error on boundary current')
err_curr= err_curr_b

disp('Keff for reflected system:')
keff_refl= k2(end)

disp('Number of iterations for the reflected system')
iteration2=it2

disp('Current at the core boundary for reflected core NUMERICAL')
current_reflected_bound= J_na_dx_refl




