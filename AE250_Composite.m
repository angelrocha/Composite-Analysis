%% Composite Laminate Stiffness & Compliance Analysis

%Given play material properties & layup configuration (number of plies &
%ply orientation angle (theta)
%Calculates laminate stiffness & compliance matrices, Off Axis & On Axis
%Stresses

clear all; close all; clc;

%% Copmosite Laminate Properties
%Laminate Configurations
lam = [0 90 0 90]*pi/180; %Cross Ply
lam = [53 -53 53 -53]*pi/180; %Angle Ply
lam1 = [45 -45 -45 45]*pi/180; %Symmetric Angle Ply
lam = [60 -60]*pi/180; %example
%Ply Properties
mat = [25.8 2.48 0.295 0.783 0.00492];
     %[E1 E2 v12 G12 t_ply]
%Cylindrical Tube Loading
p=217e-6;%Msi
q=5;%in
alpha1 = -0.7*1e-6; %1/C thermal expansion coefficient along fibers
alpha2 = 25*1e-6;   %1/C thermal expansion coefficient of matrix
s1 = 3.47;    %Msi
s2 = 10.4e-3; %Msi
s12 = 7e-3;   %Msi
dT=65.5555;   %150F to C

NMmechanical=[p*q/2;  p*q; 0;  %Nx Ny Nxy forces
              0;      0;   0]; %Mx My Mxy M=moments

%% Main Code
[Ts,Te,A,B,D,Q,QB,SB,h] = ABDmatrices(mat, lam1);
[stiffness, compliance] = StiffComp(A,B,D);
[alphaOff,TstrainOn,TstressOn,TstrainOff,TstressOff]=ThermalEffects(lam1,Te,alpha1,alpha2,dT,Q,SB);
[NMthermal,alpha] = ThermalStress(lam1,dT,h,QB,alpha1,alpha2,alphaOff);


Load=NMmechanical-NMthermal;
% Laminate Strains
lam_strain = compliance*Load;
%strain = [epsilonx epsilony gammaxy kx ky kxy]

%% Stress & Strain of Plys
for a=1:1:size(lam1,2)
    z=0;
    for q=1:1:size(lam1,2)
        for b=1:2
            c = cos(lam1(q));
            s = sin(lam1(q));
            v=b-1;
            z=z+1;
            position(z)=h(q+v);
            %Mechanical Off Axis ply strain
            MstrainOff(:,z)=lam_strain(1:3)+position(z)*lam_strain(4:6);
            %Mechanical Off Axis ply stress
            MstressOff(:,z)=QB(:,:,q)*MstrainOff(:,z);
            %Mechanical On Axis ply stress
            MstressOn(:,z)=(Ts(:,:,q))*MstressOff(:,z);
            %On Axis Total Effect
            stressOff=MstressOff;
            %Off Axis Total Effect
            stressOn=MstressOn;
        end
    end
    PstressOff(:,:,a)=stressOff*1e3;
    PstressOn(:,:,a)=stressOn*1e3;
end

%% Thermal Effect


%% Plots
figure
subplot(4,1,1)
plot(PstressOff(1,:,1),position)
title('Cross ply laminate: Off axis normal stresses')
xlabel( 'stress (Ksi)')
ylabel( ' laminate thickness (inches)')
axis([0 150 h(1) h(end) ])


subplot(4,1,2)
plot(PstressOn(1,:,1),position)
hold on
title('Cross ply laminate: On axis normal stresses (x direction)')
xlabel( 'stress (Ksi)')
ylabel( ' laminate thickness (inches)')
axis([0 150 h(1) h(end) ])

subplot(4,1,3)
plot(PstressOn(2,:,1),position)
hold on
title('Cross ply laminate: On axis normal stresses (transverse direction)')
xlabel( 'stress (Ksi)')
ylabel( ' laminate thickness (inches)')
axis([0 20 h(1) h(end) ])

subplot(4,1,4)
plot(PstressOn(3,:,1),position)
hold on
title('Cross ply laminate: On axis normal stresses (shear)')
xlabel( 'stress (Ksi)')
ylabel( ' laminate thickness (inches)')
axis([-20 20 h(1) h(end) ])

function [T_sigma,T_epsilon,A,B,D,Q,QB,SB,h] = ABDmatrices(mat,lam)
%'mat' should be a matrix which includes material properties 
%[E1 E2 v12 G12 t_ply]
%'lam' should be a matrix which describes ply layout in degrees
%ex: [-53 53 -53 53]

%Define Laminate Property Values
E1=mat(1,1); E2=mat(1,2); v12=mat(1,3); G12=mat(1,4); t_ply=mat(1,5);
t_lam=size(lam,2)*t_ply; %laminate thickness
h=linspace(-t_lam/2,t_lam/2,size(lam,2)+1); %h vector describes ply distace from center
%Define On-Axis Compliance (S) and Stiffness (Q) Matrices
S = [1/E1     -v12/E1  0;
     -v12/E1  1/E2     0;
     0        0        1/G12];
Q = inv(S);
%Define A,B,D Matrices
A=zeros(3,3); B=A; D=A;
for a=1:1:size(lam,2)
    c = cos(lam(a));
    s = sin(lam(a));
    %Stress Transformation Matrix
    T_sigma(:,:,a) = [c^2   s^2  2*c*s;
                      s^2   c^2  -2*c*s;
                      -c*s  c*s  c^2-s^2];
    %Strain Transformation Matrix
    T_epsilon(:,:,a) = [c^2     s^2    c*s;
                        s^2     c^2    -c*s;
                        -2*c*s  2*c*s  c^2-s^2];
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
    %Rearange Q_bar for A,B,D summation formula
    QB(:,:,a) = [Q_bar(1,1) Q_bar(3,1) Q_bar(4,1);
                 Q_bar(3,1) Q_bar(2,1) Q_bar(5,1);
                 Q_bar(4,1) Q_bar(5,1) Q_bar(6,1)];
    SB(:,:,a) = inv(QB(:,:,a));
    for i=1:1:3
        for j=1:1:3
            %A,B,D Matrix Summations
            A(i,j)=A(i,j)+QB(i,j,a)*(h(a+1)-h(a));
            B(i,j)=B(i,j)+(1/2)*QB(i,j,a)*((h(a+1))^2-(h(a))^2);
            D(i,j)=D(i,j)+(1/3)*QB(i,j,a)*((h(a+1))^3-(h(a))^3);
        end
    end
end
end
function [stiffness, compliance] = StiffComp(A,B,D)
%stiffness matrix transforms strains to stresses
%compliance matrix transforms stresses to strains
stiffness = [A B; B D];
compliance = inv(stiffness);
alpha = compliance(1:3,1:3);
beta = compliance(4:6,1:3);
delta = compliance(4:6,4:6);
end
function [alphaOff,TstrainOn,TstressOn,TstrainOff,TstressOff]=ThermalEffects(lam,T_epsilon,alpha1,alpha2,dT,Q,SB)
%On Axis Thermal Effects
TstrainOn=[alpha1;alpha2;0]*dT;
TstressOn=Q*[alpha1;alpha2;0]*dT;
%Off Axis Thermal Effects
for a=1:1:size(lam,2)
    SBi(:,:,a)=inv(SB(:,:,a));
    alphaOff(:,1,a)=T_epsilon(:,:,a)*[alpha1;alpha2;0]; %alphaOff=[a_x;a_y;a_s]
    TstrainOff(:,1,a)=alphaOff(:,1,a)*dT;
    TstressOff(:,:,a)=SBi(:,:,a)*alphaOff(:,1,a)*dT;
end
end
function [NMthermal,alpha] = ThermalStress(lam,dT,h,QB,alpha1,alpha2,alphaOff)
for a=1:1:size(lam,2)
    Nthermal=zeros(3,1);
    Mthermal=zeros(3,1);
    c = cos(lam(a));
    s = sin(lam(a));
    alphaxy=-2*c*s*alpha1+2*c*s*alpha2
    alpha=[alphaOff(1,1,a);alphaOff(2,1,a);alphaxy];
    %Off Axis Forces & Moments
    z=h(a+1)-h(a)
    Nthermal(:,:)=Nthermal(:,:)+dT*(h(a+1)-h(a))*QB(a).*alpha;
    Mthermal(:,:)=Mthermal(:,:)+dT*(((h(a+1))^2-(h(a))^2)/2)*QB(a)*alpha;
    NMthermal=[Nthermal;Mthermal];
end
end