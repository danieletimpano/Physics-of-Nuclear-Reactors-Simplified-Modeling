close all
clear all
clc
%% Exercise 4 - Optimization of loading pattern 
% ### Algorithm description ###

% The general script solves the problem for any INPUT loading pattern
% The permut script selects the solutions that respects the constraints
% The LUP and LLP scripts solve the initial and optimised configurations
% For the solution of the problem the following procedure was employed:

% ### PREPROCESSING - SELECTION, DATA AND MESH PROPERTIES ###

% SECTION 0 - selection of the configuration;
% SECTION 1 - Data
% - the data in the thermal and fast groups are stored in vectors
% - the index in the properties' vectors indicates the assembly type

% SECTION 2 - Building properties along the mesh
% - the properties are linked to mesh position accordingly to the chosen
%   assembly pattern, the procedure is repeated for diffusion properties,
%   total cross-sections, up and down scattering cross-sections and source
%   terms

% ### NUMERICAL SOLVER ###

% SECTION 3 - Construction of the matrix M
% - the beta coefficients are computed accordingly to the diffusion
%   properties along the mesh built in SECTION 2.
% - diagonals of the M1 and M2 matrices are built using for cycles
% - submatrices for group 1 and group 2 were built using diag properties
% - boundary conditions were set accordingly to information on the symmetry
%   of the system and on zero incoming current at the boundary for both
%   groups
% - %upscattering and downscattering terms are added
% - the F matrix is built accordingly to diag properties, considering both
%   the fast and thermal fission contribution

% SECTION 4 - Power Iteration scheme
% - set the first guess for flux and keff
% - set the source term according to the guess
% - solve the problem for the flux
% - use the flux to update keff accordingly to TP2 info
% - update the source term
% - compute the relative difference between keff at each iteration to set
%   the stop criteria

% SECTION 5 - Output and normalisation
% - output plots are produced
% - normalisation of the flux is carried out according to power info
%   as specified in the report and in TP4 intro slides

%% Info on Output
% All the output won't be shown in Command Window apart from those required
% in TP4 Slides + graphs 

%% SOLUTION FOR THE CHOSEN CONFIGURATION
% choose a given pattern
%% SECTION 0 - STARTING POINT: CHOICE OF THE PATTERN
%choose pattern
ass= [3 3 2 2 1 1 4];

%% SECTION 1 - Data
%geometric parameters for the definition of the core
a=120; %core size
b=20; %reflector size
tot= a+b; %total size
ass_size=20;
h=1; %mesh size
n_ass=ass_size/h; %mesh points in each assembly
nt=(a+b)/h; %total mesh points
height= 400; %height of the assembly

%Fast data - index of the vector equals assembly type
D1= [];
D1= [1.48240, 1.48540E+00, 1.48500E+00, 1.20000E+00 ];
sigma_a1= [9.61590e-3, 1.05770E-02, 1.11090E-02, 1.00000E-03 ];
nusigma_f1=[];
nusigma_f1= [7.16950E-03, 6.00220E-03, 5.11280E-03, 0 ];
ksigma_f1=[];
ksigma_f1= [9.17230E-14, 7.44960E-14, 6.19580E-14, 0];
sigma_s12= [];
sigma_s12= [1.73690E-02, 1.63500E-02, 1.58150E-02, 2.50000E-02 ]; %from group 1 to group 2

%Thermal Data - index of the vector equals assembly type
D2= [];
D2= [3.81380e-1, 3.70450E-01, 3.67600E-01, 4.00000E-01];
sigma_a2= [8.21530e-02, 9.56160E-02, 9.30040E-02, 2.00000E-02 ];
nusigma_f2=[];
nusigma_f2= [1.40380E-01, 1.42670E-01, 1.27650E-01, 0];
ksigma_f2=[];
ksigma_f2= [1.86440E-12, 1.78730E-12, 1.54970E-12, 0];
sigma_s21= [];
sigma_s21= [1.62710E-03, 1.84670E-03, 1.81120E-03,0.00000E+00]; %from group 2 to group 1

%% SECTION 2 - Building properties along the mesh

%vectors containing diffussion coefficients along the mesh

D_1t=[D1(ass(1))*ones(n_ass,1); D1(ass(2))*ones(n_ass,1); ...
      D1(ass(3))*ones(n_ass,1); D1(ass(4))*ones(n_ass,1); ...
      D1(ass(5))*ones(n_ass,1); D1(ass(6))*ones(n_ass,1); ...
      D1(ass(7))*ones(n_ass,1)];
D_2t=[D2(ass(1))*ones(n_ass,1); D2(ass(2))*ones(n_ass,1); ...
      D2(ass(3))*ones(n_ass,1); D2(ass(4))*ones(n_ass,1); ...
      D2(ass(5))*ones(n_ass,1); D2(ass(6))*ones(n_ass,1); ...
      D2(ass(7))*ones(n_ass,1)];

%vectors containing removal cross-sections along the mesh

sigma_t1=[ (sigma_a1(ass(1))+ sigma_s12(ass(1)))*ones(n_ass,1); ...
           (sigma_a1(ass(2))+ sigma_s12(ass(2)))*ones(n_ass,1); ...
           (sigma_a1(ass(3))+ sigma_s12(ass(3)))*ones(n_ass,1); ...
           (sigma_a1(ass(4))+ sigma_s12(ass(4)))*ones(n_ass,1); ...
           (sigma_a1(ass(5))+ sigma_s12(ass(5)))*ones(n_ass,1); ...
           (sigma_a1(ass(6))+ sigma_s12(ass(6)))*ones(n_ass,1); ...
           (sigma_a1(ass(7))+ sigma_s12(ass(7)))*ones(n_ass,1)];
sigma_t2=[ (sigma_a2(ass(1))+ sigma_s21(ass(1)))*ones(n_ass,1); ...
           (sigma_a2(ass(2))+ sigma_s21(ass(2)))*ones(n_ass,1); ...
           (sigma_a2(ass(3))+ sigma_s21(ass(3)))*ones(n_ass,1); ...
           (sigma_a2(ass(4))+ sigma_s21(ass(4)))*ones(n_ass,1); ...
           (sigma_a2(ass(5))+ sigma_s21(ass(5)))*ones(n_ass,1); ...
           (sigma_a2(ass(6))+ sigma_s21(ass(6)))*ones(n_ass,1); ...
           (sigma_a2(ass(7))+ sigma_s21(ass(7)))*ones(n_ass,1)];

%we need to account for the downscattering term from group1 to group2
down_scatter_t= [-sigma_s12(ass(1))*ones(n_ass,1); -sigma_s12(ass(2))*ones(n_ass,1);...
                 -sigma_s12(ass(3))*ones(n_ass,1); -sigma_s12(ass(4))*ones(n_ass,1); ...
                 -sigma_s12(ass(5))*ones(n_ass,1); -sigma_s12(ass(6))*ones(n_ass,1); ...
                 -sigma_s12(ass(7))*ones(n_ass,1)] ;

%we need to account for the upscattering term from group2 to group1
up_scatter_t= [  -sigma_s21(ass(1))*ones(n_ass,1); -sigma_s21(ass(2))*ones(n_ass,1);...
                 -sigma_s21(ass(3))*ones(n_ass,1); -sigma_s21(ass(4))*ones(n_ass,1); ...
                 -sigma_s21(ass(5))*ones(n_ass,1); -sigma_s21(ass(6))*ones(n_ass,1); ...
                 -sigma_s21(ass(7))*ones(n_ass,1)] ;

%source term data
ff= [(nusigma_f1(ass(1)))*ones(n_ass,1); (nusigma_f1(ass(2)))*ones(n_ass,1); ...
     (nusigma_f1(ass(3)))*ones(n_ass,1); (nusigma_f1(ass(4)))*ones(n_ass,1); ...
     (nusigma_f1(ass(5)))*ones(n_ass,1); (nusigma_f1(ass(6)))*ones(n_ass,1); ...
     (nusigma_f1(ass(7)))*ones(n_ass,1)];
kf1= [(ksigma_f1(ass(1)))*ones(n_ass,1); (ksigma_f1(ass(2)))*ones(n_ass,1); ...
    (ksigma_f1(ass(3)))*ones(n_ass,1); (ksigma_f1(ass(4)))*ones(n_ass,1); ...
    (ksigma_f1(ass(5)))*ones(n_ass,1); (ksigma_f1(ass(6)))*ones(n_ass,1); ...
    (ksigma_f1(ass(7)))*ones(n_ass,1)];
ff= [ff; zeros(nt,1)];
tf= [(nusigma_f2(ass(1)))*ones(n_ass,1); (nusigma_f2(ass(2)))*ones(n_ass,1);...
     (nusigma_f2(ass(3)))*ones(n_ass,1); (nusigma_f2(ass(4)))*ones(n_ass,1); ...
     (nusigma_f2(ass(5)))*ones(n_ass,1); (nusigma_f2(ass(6)))*ones(n_ass,1); ...
     (nusigma_f2(ass(7)))*ones(n_ass,1)];
kf2= [(ksigma_f2(ass(1)))*ones(n_ass,1); (ksigma_f2(ass(2)))*ones(n_ass,1);...
      (ksigma_f2(ass(3)))*ones(n_ass,1); (ksigma_f2(ass(4)))*ones(n_ass,1); ...
      (ksigma_f2(ass(5)))*ones(n_ass,1); (ksigma_f2(ass(6)))*ones(n_ass,1); ...
      (ksigma_f2(ass(7)))*ones(n_ass,1)];
  
             
%% SECTION 3 - Costruction of the matrix M

%explicit construction of the beta coefficients
beta_plus1t=zeros(1,nt);
beta_plus2t=zeros(1,nt);

for i=1:(nt-1) %Definition of Beta i+1
        beta_plus1t(i+1)=1.*2*D_1t(i)*D_1t(i+1)/(D_1t(i+1)+D_1t(i));
        beta_plus2t(i+1)=1.*2*D_2t(i)*D_2t(i+1)/(D_2t(i+1)+D_2t(i));
end

beta_minus1t=zeros(1,nt);
beta_minus2t=zeros(1,nt);

for i=2:nt %Definition of Beta i-1
        beta_minus1t(i-1)=1.*2*D_1t(i)*D_1t(i-1)/(D_1t(i-1)+D_1t(i));
        beta_minus2t(i-1)=1.*2*D_2t(i)*D_2t(i-1)/(D_2t(i-1)+D_2t(i));
end

M1t=zeros(nt,nt);
M2t= zeros(nt,nt);

%MATRIX GROUP 1
for i=2:(nt-1)
    M1t(i,i)= (beta_plus1t(i+1)+ beta_minus1t(i-1))/h^2 + sigma_t1(i); %ragiona su questo sigma
    M1t(i,i+1)= -beta_plus1t(i+1)/h^2;
    M1t(i,i-1)= -beta_minus1t(i-1)/h^2;
end

%setting the boundary conditions
BC1t= 0.5*(1/(h/(4*D1(end)) + 1));

%implementing the boundary condition in the matrix for group1 
M1t(1,1)= beta_plus1t(2)/h^2 + sigma_t1(1);
M1t(1,2)= M1t(2,3);
M1t(nt,nt)= (BC1t*h + beta_minus1t(end-1))/h^2 + sigma_t1(end);
M1t(nt,nt-1)= M1t(nt-1,nt-2);

%MATRIX GROUP 2
for i=2:(nt-1)
    M2t(i,i)= (beta_plus2t(i+1)+ beta_minus2t(i-1))/h^2 + sigma_t2(i); %ragiona su questo sigma
    M2t(i,i+1)= -beta_plus2t(i+1)/h^2;
    M2t(i,i-1)= -beta_minus2t(i-1)/h^2;
end

%setting the boundary conditions
BC2t= 0.5*(1/(h/(4*D2(end)) + 1));

%boundary conditions 
M2t(1,1)= beta_plus2t(2)/h^2 + sigma_t2(1);
M2t(1,2)=M2t(2,3);
M2t(nt,nt)= (BC2t*h + beta_minus2t(end-1))/h^2 + sigma_t2(end);
M2t(nt,nt-1)= M2t(nt-1,nt-2);

%There is a need to concatenate the two matrices
%since the flux vector will now accomodate the two groups
Mt= zeros(2*nt,2*nt);
Mt= [M1t, zeros(nt,nt); zeros(nt,nt), M2t];
%upscattering and downscattering terms are added
Mt= Mt + diag(down_scatter_t, -nt) + diag(up_scatter_t, nt);

%definition of the F matrix
Ft= zeros(2*nt,2*nt);
Ft= Ft + diag(ff,0) + diag(tf,nt);

%% SECTION 4 - Power Iteration
% 
%first guess 
phi=[];
k=[];
S=zeros(2*nt,2*nt);
epsilon= 1e-7;
% 
% set the first guess
diff=1; %set diff=1 to enter the cycle
it=1;
k(1)=1;
phi(:,1)=ones(2*nt,1);
S(:,1)= 1/k(1)*Ft*phi(:,1);

%two possible convergence criteria were taken into account
while (diff > epsilon)
phi(:,it+1)= Mt\S(:,it);
%k(it+1)= k(it)*sum(Ft*phi(:,it+1))/sum(Ft *phi(:,it)); %option 2 for convergence
k(it+1)= k(it)*sum(phi(:,it+1))/sum(phi(:,it)); %option 1 for convergence
S(:,it+1)= 1/k(it+1)*Ft*phi(:,it+1); %update the source term
diff= abs((k(it+1)-k(it))/(k(it+1))); %stopping criteria
it= it+1;
end 

%% "Exact solution" for comparison
% % analytical solution for comparison keff
% input2= Mt\Ft;
% [PHIt, KEFF_anal_t]= eigs(input2);
% keff_t= KEFF_anal_t(1,1);

%% SECTION 5 - Output and normalisation
pos_t= [0.5:1:139.5];

% figure(2)
% plot(pos_t,PHIt((1:140),end),pos_t,PHIt((141:280),end), 'linewidth', 2)
% title('Analytical solution')
% legend('Flux Analytical Group 1', 'Flux Analytical Group 2','Location','NorthEast')
% grid on

figure(1)
plot(pos_t,phi((1:140),end),pos_t,phi((141:280),end), 'linewidth', 2)
xline(20, 'g--')
xline(40, 'b--')
xline(60, 'g--')
xline(80, 'b--')
xline(100, 'r--')
xline(120, 'r--')
xline(140, 'k--')
title('Numerical solution selected configuration')
grid on

%% NORMALISATION PROCEDURE
P=[];
for i= 1:140
    P(i)= height*1*(kf1(i)*phi(i,end) + kf2(i)*phi((i+140),end));
end

%compute total power
P_tot= sum(P);

%relative power computation
P_mean= P_tot/120; %compute mean power along the mesh
P_rel= P/P_mean; %compute relative power with respect to the mean
Peak= max(P_rel); %select peak power as maximum of the relative power

%normalization factor
F_norm= (45e6)/P_tot; %compute the normalisation factor

%normalized flux
phi_norm= F_norm*phi(:,end);

%% PLOTTING - NORMALIZED FLUX
figure(2)
plot(pos_t,phi_norm((1:140),end),pos_t,phi_norm((141:280),end), 'linewidth', 2)
xline(20, 'g--')
xline(40, 'b--')
xline(60, 'g--')
xline(80, 'b--')
xline(100, 'r--')
xline(120, 'r--')
xline(140, 'k--')
title('Normalised solution selected configuration')
grid on

%% PLOTTING - RELATIVE POWER
figure(3)
plot(pos_t,P_rel, 'linewidth', 2)
xline(20, 'g--')
xline(40, 'b--')
xline(60, 'g--')
xline(80, 'b--')
xline(100, 'r--')
xline(120, 'r--')
xline(140, 'k--')
title('Relative Power selected configuration')
grid on

%% REQUESTED OUTPUT

format long
disp('SOLUTION FOR THE LOW-LEAKAGE PATTERN:')
disp('The solution in 10.5 for normalised flux FAST:')
phi_norm(11,end)
disp('The solution in 10.5 for normalised flux THERMAL:')
phi_norm(151,end)
disp('The flux at the reflector boundary:')
phi_norm(140, end)
disp('Keff is:')
k(end)
disp('Peak power relative to the average:')
Peak


