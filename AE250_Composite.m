%% Composite Laminate Stiffness & Compliance Analysis

%Given play material properties & layup configuration (number of plies &
%ply orientation angle (theta)
%Calculates laminate stiffness & compliance matrices, Off Axis & On Axis
%Stresses

clear all; close all; clc;

%% Copmosite Laminate Properties
%Laminate Configurations
lam1 = [0 90 0 90]*pi/180; %Cross Ply
lam2 = [53 -53 53 -53]*pi/180; %Angle Ply
lam3 = [45 -45 -45 45]*pi/180; %Symmetric Angle Ply
%Ply Properties
mat = [25.8 2.48 0.295 0.783 0.00492];
     %[E1 E2 v12 G12 t_ply]
%Cylindrical Tube Loading
p=217e-6;%Msi
q=5;%in
N=[p*q/2;p*q;0] %[Nx Ny Nxy] forces
M=[0;0;0];      %[Mx My Mxy] M=moments
Load=[N;M]
       
s1 = 3.47;%Msi
s2 = 10.4e-3;%Msi
s12 = 7e-3;%Msi

%% A,B,D matrices wrt Laminate Midplane
[A,B,D,QB] = ABDmatrices(mat, lam3) 

%% alpha, beta, delta matrices
stiffness = [A B; B D];
compliance = inv(stiffness);
alpha = compliance(1:3,1:3);
beta = compliance(4:6,1:3);
delta = compliance(4:6,4:6);

%% Laminate Strains
lam_strain = compliance*Load
%strain = [epsilonx epsilony gammaxy kx ky kxy]

%% Stress & Strain of Plys
h=linspace(-mat(5)*size(lam1,2)/2,mat(5)*size(lam1,2)/2,size(lam1,2)+1); 
for a=1:1:size(lam1,2)      
    z=0;
    for q=1:1:size(lam1,2)
        for b=1:2
            c = cos(lam1(a));
            s = sin(lam1(a));
            %Stress Transformation Matrix
             T_sigma(:,:,q) = [c^2   s^2  2*c*s;
                               s^2   c^2  -2*c*s;
                               -c*s  c*s  c^2-s^2];
            v=b-1;
            z=z+1;
            position(z)=h(q+v);
            %Off Axis ply strain
            ply_strain_off(:,z)=lam_strain(1:3)+position(z)*lam_strain(4:6);
            %Off Axis ply stress
            ply_stress_off(:,z)=QB(:,:,q)*ply_strain_off(:,z);
            %On Axis ply stress
            ply_stress_on(:,z)=inv(T_sigma(:,:,q))*ply_stress_off(:,z);
        end
    end
    stress_off(:,:,a)=ply_stress_off*1e3;
    stress_on(:,:,a)=ply_stress_on*1e3;
end

%% Plots
figure
subplot(3,1,1)
plot(stress_off(1,:,1),position)
title('Cross ply laminate: Off axis normal stresses')
xlabel( 'stress (Ksi)')
ylabel( ' laminate thickness (inches)')
axis([0 150 h(1) h(end) ])


subplot(3,1,2)
plot(stress_on(2,:,1),position)
hold on
title('Cross ply laminate: On axis normal stresses')
xlabel( 'stress (Ksi)')
ylabel( ' laminate thickness (inches)')
axis([0 150 h(1) h(end) ])

function [A,B,D,QB] = ABDmatrices(mat,lam)
    %'mat' should be a matrix which includes material properties 
    %[E1 E2 v12 G12 t_ply]
    %'lam' should be a matrix which describes ply layout in degrees
    %ex: [-53 53 -53 53]
    %multiple laminate layouts may be processed with multiple row matrix
    
    %Define Laminate Property Values
    E1=mat(1,1); E2=mat(1,2); v12=mat(1,3); G12=mat(1,4); t_ply=mat(1,5);
    [m,n]=size(lam); %n=number of plys
    t_lam=n*t_ply; %laminate thickness
    h=linspace(-t_lam/2,t_lam/2,n+1) %h vector describes ply distace from center
    %Define On-Axis Compliance (S) and Stiffness (Q) Matrices
    S = [1/E1     -v12/E1  0;
         -v12/E1  1/E2     0;
         0        0        1/G12];
    Q = inv(S);
    %Define A,B,D Matrices
    A=zeros(3,3); B=A; D=A;
    for k=1:1:n
        c = cos(lam(k));
        s = sin(lam(k));
        %Define Off-Axis Stiffness Transition Matrix  
        A_bar = [c^4      s^4      2*c^2*s^2    4*c^2*s^2;
                 s^4      c^4      2*c^2*s^2    4*c^2*s^2;
                 c^2*s^2  c^2*s^2  c^4+s^4      -4*c^2*s^2;
                 c^3*s    -c*s^3   c*s^3-c^3*s  2*(c*s^3-c^3*s);
                 c*s^3    -c^3*s   c^3*s-c*s^3  2*(c^3*s-c*s^3);
                 c^2*s^2  c^2*s^2  -2*c^2*s^2   (c^2-s^2)^2];
        %Define Off-Axis Compliance Transition Matrix
        B_bar = [c^4        s^4        2*c^2*s^2        c^2*s^2;
                 s^4        c^4        2*c^2*s^2        c^2*s^2;
                 c^2*s^2    c^2*s^2    c^4+s^4          -c^2*s^2;
                 4*c^2*s^2  4*c^2*s^2  -8*c^2*s^2       (c^2-s^2)^2;
                 2*c^3*s    -2*c*s^3   2*(c*s^3-c^3*s)  (c*s^3-c^3*s);
                 2*c*s^3    -2*c^3*s   2*(c^3*s-c*s^3)  c^3*s-c*s^3];
        %Off Axis Stiffness Matrices
        Q_bar(:,:) = A_bar*[Q(1,1); Q(2,2); Q(1,2); Q(3,3)];
        %Off Axis Compliance Matrces
        S_bar(:,:) = B_bar*[S(1,1); S(2,2); S(1,2); S(3,3)];
        %Rearange Q_bar for A,B,D summation formula
        QB(:,:,k) = [Q_bar(1,1) Q_bar(3,1) Q_bar(4,1);
                     Q_bar(3,1) Q_bar(2,1) Q_bar(5,1);
                     Q_bar(4,1) Q_bar(5,1) Q_bar(6,1)];
        for i=1:1:3
            for j=1:1:3
                %A,B,D Matrix Summations
                A(i,j)=A(i,j)+QB(i,j,k)*(h(k+1)-h(k));
                B(i,j)=B(i,j)+(1/2)*QB(i,j,k)*(h(k+1)^2-h(k)^2);
                D(i,j)=D(i,j)+(1/3)*QB(i,j,k)*(h(k+1)^3-h(k)^3);
            end
        end
    end
end
