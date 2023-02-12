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

% SECTION 3 - Construction of the matrix M and F
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
% in TP4 printouts + graphs 

%% SOLUTION FOR THE GIVEN LOW UNIFORMITY PATTERN
% code for the starting case in TP4 slides

%% SECTION 0 - STARTING POINT
%low uniformity pattern
% [3 3 2 2 1 1 4]

pattern= [3 3 2 2 1 1];
permut_pattern = perms(pattern);
permut_pattern = [permut_pattern, 4*ones(size(permut_pattern,1),1)];

%% SECTION 1 - Data

%geometric parameters for the definition of the core
a=120;
b=20;
tot= a+b;
ass=20;
h=1;
n_ass=ass/h;
nt=(a+b)/h; 
height= 400;

%ID 1
D1= [];
D1(1)= 1.48240;
D1(2)= 3.81380e-1;
sigma_a1= [];
sigma_a1(1)= 9.61590e-3;
sigma_a1(2)= 8.21530e-02;
nusigma_f1=[];
nusigma_f1(1)= 7.16950E-03;
nusigma_f1(2)= 1.40380E-01;
ksigma_f1=[];
ksigma_f1(1)= 9.17230E-14;
ksigma_f1(2)= 1.86440E-12;
sigma1_s1= [];
sigma1_s1(1)= 1.97880E-01;
sigma1_s1(2)= 1.73690E-02;
sigma1_s2= [];
sigma1_s2(1)= 1.62710E-03;
sigma1_s2(2)= 7.90230E-01;

%ID 2
D2= [];
D2(1)= 1.48540E+00;
D2(2)= 3.70450E-01;
sigma_a2= [];
sigma_a2(1)= 1.05770E-02;
sigma_a2(2)= 9.56160E-02;
nusigma_f2=[];
nusigma_f2(1)= 6.00220E-03;
nusigma_f2(2)= 1.42670E-01;
ksigma_f2=[];
ksigma_f2(1)= 7.44960E-14;
ksigma_f2(2)= 1.78730E-12;
sigma2_s1= [];
sigma2_s1(1)= 1.97470E-01;
sigma2_s1(2)= 1.63500E-02;
sigma2_s2= [];
sigma2_s2(1)= 1.84670E-03;
sigma2_s2(2)= 8.02350E-01;

%ID 3
D3= [];
D3(1)= 1.48500E+00;
D3(2)= 3.67600E-01;
sigma_a3= [];
sigma_a3(1)= 1.11090E-02;
sigma_a3(2)= 9.30040E-02;
nusigma_f3=[];
nusigma_f3(1)= 5.11280E-03;
nusigma_f3(2)= 1.27650E-01;
ksigma_f3=[];
ksigma_f3(1)= 6.19580E-14;
ksigma_f3(2)= 1.54970E-12;
sigma3_s1= [];
sigma3_s1(1)= 1.97540E-01;
sigma3_s1(2)= 1.58150E-02;
sigma3_s2= [];
sigma3_s2(1)= 1.81120E-03;
sigma3_s2(2)= 8.11960E-01;

%ID4
D4= [];
D4(1)= 1.20000E+00;
D4(2)= 4.00000E-01;
sigma_a4= [];
sigma_a4(1)= 1.00000E-03;
sigma_a4(2)= 2.00000E-02;
nusigma_f4=[];
nusigma_f4(1)= 0;
nusigma_f4(2)= 0;
ksigma_f4=[];
ksigma_f4(1)= 0;
ksigma_f4(2)= 0;
sigma4_s1= [];
sigma4_s1(1)= 2.51780E-01;
sigma4_s1(2)= 2.50000E-02;
sigma4_s2= [];
sigma4_s2(1)= 0.00000E+00;
sigma4_s2(2)= 8.13330E-01;

%% SECTION 2 - Building properties along the mesh

%vectors containing diffussion coefficients along the mesh

D_1t=[D3(1)*ones(n_ass,1); D3(1)*ones(n_ass,1); ...
    D2(1)*ones(n_ass,1); D2(1)*ones(n_ass,1); ...
    D1(1)*ones(n_ass,1); D1(1)*ones(n_ass,1); ...
    D4(1)*ones(n_ass,1)];
D_2t=[D3(2)*ones(n_ass,1); D3(2)*ones(n_ass,1); ...
    D2(2)*ones(n_ass,1); D2(2)*ones(n_ass,1); ...
    D1(2)*ones(n_ass,1); D1(2)*ones(n_ass,1); ...
    D4(2)*ones(n_ass,1)];

%vectors containing removal cross-sections along the mesh

sigma_t1=[(sigma_a3(1)+ sigma3_s1(2))*ones(n_ass,1); ...
          (sigma_a3(1)+ sigma3_s1(2))*ones(n_ass,1); ...
         (sigma_a2(1)+ sigma2_s1(2))*ones(n_ass,1); ...
         (sigma_a2(1)+ sigma2_s1(2))*ones(n_ass,1); ...
         (sigma_a1(1)+ sigma1_s1(2))*ones(n_ass,1); ...
         (sigma_a1(1)+ sigma1_s1(2))*ones(n_ass,1); ...
         (sigma_a4(1)+ sigma4_s1(2))*ones(n_ass,1)]; ...
sigma_t2=[(sigma_a3(2)+ sigma3_s2(1))*ones(n_ass,1); ...
          (sigma_a3(2)+ sigma3_s2(1))*ones(n_ass,1); ...
         (sigma_a2(2)+ sigma2_s2(1))*ones(n_ass,1); ...
         (sigma_a2(2)+ sigma2_s2(1))*ones(n_ass,1); ...
         (sigma_a1(2)+ sigma1_s2(1))*ones(n_ass,1); ...
         (sigma_a1(2)+ sigma1_s2(1))*ones(n_ass,1); ...
         (sigma_a4(2)+ sigma4_s2(1))*ones(n_ass,1)];


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

%% SECTION 3 - Building matrix M and F

M1t=zeros(nt,nt);
M2t= zeros(nt,nt);

%MATRIX GROUP 1
for i=2:(nt-1)
    M1t(i,i)= (beta_plus1t(i+1)+ beta_minus1t(i-1))/h^2 + sigma_t1(i); %ragiona su questo sigma
    M1t(i,i+1)= -beta_plus1t(i+1)/h^2;
    M1t(i,i-1)= -beta_minus1t(i-1)/h^2;
end


%setting the boundary conditions
BC1t= 0.5*(1/(h/(4*D4(1)) + 1));

%boundary condition + filling the matrix holes due to for cycle
%construction
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
BC2t= 0.5*(1/(h/(4*D4(2)) + 1));

%boundary condition + filling the matrix holes due to for cycle
%construction
M2t(1,1)= beta_plus2t(2)/h^2 + sigma_t2(1);
M2t(1,2)=M2t(2,3);
M2t(nt,nt)= (BC2t*h + beta_minus2t(end-1))/h^2 + sigma_t2(end);
M2t(nt,nt-1)= M2t(nt-1,nt-2);

%There is a need to concatenate the two matrices
%since the flux vector will now accomodate the two groups
Mt= zeros(2*nt,2*nt);
Mt= [M1t, zeros(nt,nt); zeros(nt,nt), M2t];

%we need to account for the downscattering term from group1 to group2
down_scatter_t= [-sigma3_s1(2)*ones(n_ass,1); -sigma3_s1(2)*ones(n_ass,1); ...
    -sigma2_s1(2)*ones(n_ass,1); -sigma2_s1(2)*ones(n_ass,1); ...
    -sigma1_s1(2)*ones(n_ass,1); -sigma1_s1(2)*ones(n_ass,1); ...
    -sigma4_s1(2)*ones(n_ass,1)] ;

%we need to account for the upscattering term from group2 to group1
up_scatter_t= [-sigma3_s2(1)*ones(n_ass,1); -sigma3_s2(1)*ones(n_ass,1); ...
    -sigma2_s2(1)*ones(n_ass,1); -sigma2_s2(1)*ones(n_ass,1); ...
    -sigma1_s2(1)*ones(n_ass,1); -sigma1_s2(1)*ones(n_ass,1); ...
    -sigma4_s2(1)*ones(n_ass,1)] ;

Mt= Mt + diag(down_scatter_t, -nt) + diag(up_scatter_t, nt);

%% Costruction of the F matrix 
ff= [(nusigma_f3(1))*ones(n_ass,1); (nusigma_f3(1))*ones(n_ass,1); ...
    (nusigma_f2(1))*ones(n_ass,1); (nusigma_f2(1))*ones(n_ass,1); ...
    (nusigma_f1(1))*ones(n_ass,1); (nusigma_f1(1))*ones(n_ass,1); ...
    (nusigma_f4(1))*ones(n_ass,1)];
kf1= [(ksigma_f3(1))*ones(n_ass,1); (ksigma_f3(1))*ones(n_ass,1); ...
    (ksigma_f2(1))*ones(n_ass,1); (ksigma_f2(1))*ones(n_ass,1); ...
    (ksigma_f1(1))*ones(n_ass,1); (ksigma_f1(1))*ones(n_ass,1); ...
    (ksigma_f4(1))*ones(n_ass,1)];
ff= [ff; zeros(nt,1)];
tf= [(nusigma_f3(2))*ones(n_ass,1); (nusigma_f3(2))*ones(n_ass,1); ...
    (nusigma_f2(2))*ones(n_ass,1); (nusigma_f2(2))*ones(n_ass,1); ...
    (nusigma_f1(2))*ones(n_ass,1); (nusigma_f1(2))*ones(n_ass,1); ...
    (nusigma_f4(2))*ones(n_ass,1)];
kf2= [(ksigma_f3(2))*ones(n_ass,1); (ksigma_f3(2))*ones(n_ass,1); ...
    (ksigma_f2(2))*ones(n_ass,1); (ksigma_f2(2))*ones(n_ass,1); ...
    (ksigma_f1(2))*ones(n_ass,1); (ksigma_f1(2))*ones(n_ass,1); ...
    (ksigma_f4(2))*ones(n_ass,1)];
Ft= zeros(2*nt,2*nt);
Ft= Ft + diag(ff,0) + diag(tf,nt);

%% Section 4 - Power Iteration method
% 
%first guess 
phi=[];
k=[];
S=zeros(2*nt,2*nt);
epsilon= 1e-7;
% 
% set the first guess
diff=1;
it=1;
k(1)=1;
phi(:,1)=ones(2*nt,1);
S(:,1)= 1/k(1)*Ft*phi(:,1);


while (diff > epsilon)
phi(:,it+1)= Mt\S(:,it);
%k(it+1)= k(it)*sum(Ft*phi(:,it+1))/sum(Ft *phi(:,it)); %option 2 for convergence
k(it+1)= k(it)*sum(phi(:,it+1))/sum(phi(:,it)); %option 1 for convergence
S(:,it+1)= 1/k(it+1)*Ft*phi(:,it+1);
diff= abs((k(it+1)-k(it))/(k(it+1)));
it= it+1;
end 

%% "Exact solution" for comparison
% % analytical solution for comparison keff
% input2= Mt\Ft;
% [PHIt, KEFF_anal_t]= eigs(input2);
% keff_t= KEFF_anal_t(1,1);



%% Section 5 - Output and normalisation
pos_t= [0.5:1:139.5];

% figure(2)
% plot(pos_t,PHIt((1:140),end),pos_t,PHIt((141:280),end), 'linewidth', 2)
% title('Analytical solution')
% legend('Flux Analytical Group 1', 'Flux Analytical Group 2','Location','NorthEast')
% grid on

figure(1)
plot(pos_t,phi((1:140),end),pos_t,phi((141:280),end), 'linewidth', 2)
xline(20, 'r--')
xline(40, 'r--')
xline(60, 'b--')
xline(80, 'b--')
xline(100, 'g--')
xline(120, 'g--')
xline(140, 'k--')
title('Numerical solution - Low-non uniformity pattern')
legend('Flux Numerical Group 1', 'Flux Numerical Group 2', 'Assembly 3', ...
    'Assembly 3', 'Assembly 2', 'Assembly 2', 'Assembly 1', ' Assembly 1', ...
    'Reflector', 'Location','NorthWest', 'FontSize', 6)
grid on


%% NORMALISATION PROCEDURE
P=[];
for i= 1:140
    P(i)= height*1*(kf1(i)*phi(i,end) + kf2(i)*phi((i+140),end));
end

P_tot= sum(P);

%relative power computation
P_mean= P_tot/120;
P_rel= P/P_mean;
Peak= max(P_rel);

%normalization factor
F_norm= (45e6)/P_tot;

%normalized flux
phi_norm= F_norm*phi(:,end);

%% PLOTTING - NORMALIZED FLUX
figure(2)
plot(pos_t,phi_norm((1:140),end),pos_t,phi_norm((141:280),end), 'linewidth', 2)
xline(20, 'r--', 'linewidth', 2);
xline(40, 'r--', 'linewidth', 2);
xline(60, 'b--', 'linewidth', 2);
xline(80, 'b--', 'linewidth', 2);
xline(100, 'g--', 'linewidth', 2);
xline(120, 'g--', 'linewidth', 2);
xline(140, 'k--', 'linewidth', 2);
title('Normalised solution - Low non-uniformity pattern')
legend('Flux Normalised Group 1', 'Flux Normalised Group 2', 'Assembly 3', ...
    'Assembly 3', 'Assembly 2', 'Assembly 2', 'Assembly 1', ' Assembly 1', ...
    'Reflector', 'Location','NorthWest', 'FontSize', 6)
grid on
% saveas(gcf,'LUP_flux.png')
% saveas(gcf,'LUP_flux','epsc')

%% PLOTTING - RELATIVE POWER
figure(3)
plot(pos_t,P_rel, 'linewidth', 3)
xline(20, 'r--', 'linewidth', 2);
xline(40, 'r--', 'linewidth', 2);
xline(60, 'b--', 'linewidth', 2);
xline(80, 'b--', 'linewidth', 2);
xline(100, 'g--', 'linewidth', 2);
xline(120, 'g--', 'linewidth', 2);
xline(140, 'k--', 'linewidth', 2);
title('Relative power - Low non-uniformity pattern')
legend('Relative Power', 'Assembly 3', ...
    'Assembly 3', 'Assembly 2', 'Assembly 2', 'Assembly 1', ' Assembly 1', ...
    'Reflector', 'Location','NorthWest', 'FontSize', 6)
grid on
% saveas(gcf,'LUP_relativepower.png')
% saveas(gcf,'LUP_relativepower','epsc')

%% HINT PARAMETERS

format long
disp('SOLUTION FOR THE LOW-NON UNIFORMITY PATTERN:')
disp('Number of iterations:')
it
disp('The solution in 10.5 for normalised flux FAST:')
phi_norm(11,end)
disp('The solution in 10.5 for normalised flux THERMAL:')
phi_norm(151,end)
disp('The fast flux at the reflector boundary:')
phi_norm(140, end)
disp('Keff is:')
k(end)
disp('Peak power relative to the average:')
Peak


