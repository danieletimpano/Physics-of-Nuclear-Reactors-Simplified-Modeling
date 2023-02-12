close all
clear all 
clc

%% EXERCISE 1 - Attenuation of neutrons from a 1D planar source
% ### Algorithm description ###

% For the solution of the problem the following procedure was employed:

% ### PREPROCESSING - ANALYTICAL SOLUTION ###
% SECTION 1 - data from the TP scriptum were reported;
% SECTION 2 - analytical solution was written as a function handle and analytical solution was plotted;
% SECTION 2 - OUTPUT analytical solution values for x=0 and x=x0 were asked as an output;

% ### SOLUTION OF THE PROBLEM ###
% SECTION 3 - MATRIX AND SYSTEM SOLUTION BUILDING 
% n was set to 100 (i.e. size of the matrix = x0 / mesh size | mesh size = 0.1 cm)
% matrix A was assembled using diag command;
% values were changed for A(1,1) and A(n,n) accordingly to BCs;
% known values vector was set to zero apart from the first value which equals ?C/?x
% SECTION 4 - Matlab solver for the linear system was used (i.e. x= A \ b);
% SECTION 4 - OUTPUT Flux values were plotted using a linspace defined as linspace
% (h/2, x0-h/2, n) were h is the size of the mesh;
% SECTION 5 - The value of the flux in x0 and 0 was obtained using the information on the BCs

% ### ERROR COMPUTATION ###
% SECTION 6 - The error along the mesh was computed subtracting vector
% solution to linspace analytical solution evaluated in mesh points. The
% output plots do not include the values at the boundary.
% SECTION 7 - Repeats the previous steps (mesh building and solution) using different mesh sizes
% The error is stored at each iteration step, considering both the infinite norm error (max error) and the error at the boundary x0, as requested.
% Difference between analytical and numerical solution was finally plotted and printed

%% Info on Output
% All the output won't be shown in Command Window apart from those required
% in TP1 Slides + graphs 

%% SECTION 1 - Data
sigma_a = 0.02; %cm-1
sigma_s= 4; %cm-1
x0= 10; %cm 
S= 1000; %n cm-2 s-1

D= 1/(3*(sigma_s+sigma_a)); %cm
L= (1/(3*(sigma_s+sigma_a)*sigma_a))^0.5; %cm
d= 2*D; %cm

%% SECTION 2 - exact analytical solution
x= linspace(0,10, 1000); %definition of the linspace for plotting purpose
phi= @(x)(S*L)/(2*D*cosh((x0+d)/L))*sinh((x0+d-x)/L); %use of the function handle to define the analytical solution
figure(1) %plot of the analytical solution
plot(x,phi(x), 'linewidth',2 )
grid on
title("Analytical solution")
xlabel('distance [cm]')
ylabel('Flux [n.s-1.cm-2]')
% saveas(gcf,'Analytical.png')
% saveas(gcf,'Analytical','epsc')

%% SECTION 2 - OUTPUT - exact values of the flux
format long
disp('Analytical solution at 0 and x0=')
phi_0= phi(0)
phi_x0= phi(x0) 

%% SECTION 3 - Construction of the matrix
C= S/2; %definition of the known source term
n=100; %size of the matrix
h= x0/n; %size of the mesh
main_diag= (-2*D/h^2 - sigma_a)*ones(n,1); %definition of the main diagonal of the matrix
up_diag= D/h^2*ones(n-1,1); %definition of the upper diagonal of the matrix
low_diag= D/h^2*ones(n-1,1); %definition of the lower diagonal of the matrix

%costruction of the matrix summing diagonal matrices
A= diag(main_diag,0) + diag(up_diag,1) + diag(low_diag, -1);

%setting boundary conditon at x=0
A(1,1)= -D/h^2 - sigma_a;
b= zeros(n,1);
b(1)= - C/h;
%setting boundary conditions at x=x0
A(n,n)= (2*D)/h^2*(1/(h/(4*D) + 1))-(3*D/h^2) - sigma_a;

%solving the system using standard MATLAB approach
sol=A\b;

%plot the numerical solution
pos= linspace(h/2,x0-h/2,n);

%% SECTION 4 - plot of the numerical solution
figure(2)
plot(pos,sol, 'linewidth', 2, 'color','r')
grid on
title("Numerical solution")
xlabel('distance [cm]')
ylabel('Flux [n.s-1.cm-2]')

%% SECTION 4 - OUTPUT - plot of the two together
figure(2)
plot(x, phi(x),pos,sol, 'linewidth', 2)
grid on
title("Comparison between Analytical and Numerical solution")
xlabel('distance [cm]')
ylabel('Flux [n.s-1.cm-2]')
legend('Analytical solution', 'Numerical solution')
% saveas(gcf,'Numerical_vs_Analytical.png')
% saveas(gcf,'Numerical_vs_Analytical','epsc')

%% SECTION 5 - flux values for 0 and x0 for the numerical solution

%application of the boundary conditions to find the flux
disp('Flux value in 0 and x0 for 0.1 mesh - NUMERICAL SOLUTION')
phi_n_0= (h*C/(2*D))+ sol(1)
phi_n_x0= (1/(0.25 + D/h))*(D/h)*sol(100)

%computation of the error with respect to the analytical solution
disp('Flux error in 0 and x0 for 0.1 mesh - NUMERICAL SOLUTION')
err_0= abs(phi(0)-phi_n_0)
err_x0= abs(phi(x0)-phi_n_x0)

%% SECTION 6 - plot the error along the mesh for size 0.1 (i.e. n=100)
%FOLLOWING UNCOMMENTED PLOT DON'T INCLUDE THE ERROR AT THE BOUNDARY
%if I don't include the source and phi(x0)
err= abs(sol' - phi(pos));

%Flux value and error for sample point as shown in TP1 Slides
disp('Flux value and error for sample point as shown in TP1 Slides:')
sol(11)
err(11)

%if I do include the source and phi(x0)
%err= [abs(phi(0)-phi_n_0),abs(sol' - phi(pos)),abs(phi(x0)-phi_n_x0)];

%if I don't include the source and phi(x0)
pos1=  pos;

%if I do include the source and phi(x0)
%pos1= [0, pos, 10];

figure(3)
plot(pos1, err, 'linewidth', 2, 'color', 'g')
grid on
title("Error as a function of mesh position")
xlabel('distance [cm]')
ylabel('error [n.s-1.cm-2]')
% saveas(gcf,'Error_mesh_position.png')
% saveas(gcf,'Error_mesh_position','epsc')

%% SECTION 7 - find the error as a function of mesh size

% %Repeat all the previous steps for several mesh refinements
err_h=[];
err_phix0=[];
meshsize=[];

% % Construction of the matrix
for i=[1 2 3 4 5 6 7 8 9 10]
C= S/2;
n=2^i;
h= x0/n;
meshsize=[meshsize, h];
main_diag= (-2*D/h^2 - sigma_a)*ones(n,1);
up_diag= D/h^2*ones(n-1,1);
low_diag= D/h^2*ones(n-1,1);

% %costruction of the matrix
A= diag(main_diag,0) + diag(up_diag,1) + diag(low_diag, -1);
% %boundary conditon at x=0
A(1,1)= -D/h^2 - sigma_a;
b= zeros(n,1);
b(1)= - C/h;
% %boundary conditions at x=x0
A(n,n)= (2*D)/h^2*(1/(h/(4*D) + 1))-(3*D/h^2) - sigma_a;

% %solving the system
sol=A\b;

phi_nit_0= (h*C/(2*D))+ sol(1); %compute the solution in 0 using BCs
phi_nit_x0= (1/(0.25 + D/h))*(D/h)*sol(end); %compute the solution in x0 using BCs

err_it_0= abs(phi(0)-phi_nit_0);
err_it_x0= abs(phi(x0)-phi_nit_x0);
 
%vector of points in the mesh
pos= linspace(h/2,x0-h/2,n);

% %find the error
%if I choose to include the source and x0
% err_max= max([err_it_0, abs(sol'-phi(pos)), err_it_x0]); 
%if I don't include the source and x0
err_max= max(abs(sol'-phi(pos)));
err_h=[err_h, err_max];
err_phix0=[err_phix0, err_it_x0];
end 
 
%plot the infinite norm error
figure(4)
loglog(meshsize, err_h, 'linewidth', 2)
hold on
loglog(err_h, err_h.^2, 'linewidth', 2, 'LineStyle', '--') % order 2 is plotted to check convergence order
xlim([10^-2 10^1])
grid on
title("Infinite norm Error as a function of mesh size (log scale)")
xlabel('mesh size [cm]')
ylabel('error [n.s-1.cm-2]')
legend('error', 'order 2', 'Location', 'SouthEast')
% % saveas(gcf,'Error_mesh_size.png')
% % saveas(gcf,'Error_mesh_size','epsc')

%plot the error at x0 as a function of mesh size
figure(5)
loglog(meshsize, err_phix0, 'linewidth', 2)
hold on
loglog(err_h, err_h.^2, 'linewidth', 2, 'LineStyle', '--') % order 2 is plotted to check convergence order
xlim([10^-2 10^1])
grid on
title("Phi(x0) Error as a function of mesh size (log scale)")
xlabel('mesh size [cm]')
ylabel('error [n.s-1.cm-2]')
legend('error', 'order 2', 'Location', 'SouthEast')
%saveas(gcf,'Error_phi_x0.png')
%saveas(gcf,'Error_phi_x0','epsc')
