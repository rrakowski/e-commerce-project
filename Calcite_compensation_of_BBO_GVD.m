% Script calculates calcite's compensation of group velocity delay (GVM) using BBO (type_I negative uniaxial) crystal,
% acceptance angle, ... all for 100 um BBO
% R.R. 22.May.2012

clear all
close all

% Fundamental wavelength
lambda1=0.8;                         % [um]

% Second harmonic
lambda2=0.4;                         % [um]
theta=(29.2/180)*pi;

the1=0;
the2=pi/2;                           % [rad]
thestep=1000;
d_lambda=1/thestep;
d_theta=(the2-the1)/thestep;

% BBO - Sellmeier index equations from EKSMA (rounded from indexinfo site)

	no  = sqrt(2.7405 + 0.0184 / (lambda1^2-0.0179) - 0.0155*lambda1^2); 
    ne  = sqrt(2.3730 + 0.0128 / (lambda1^2-0.0156) - 0.0044*lambda1^2);
    n2o = sqrt(2.7405 + 0.0184 / (lambda2^2-0.0179) - 0.0155*lambda2^2);
    n2e = sqrt(2.3730 + 0.0128 / (lambda2^2-0.0156) - 0.0044*lambda2^2);    
   
% Angular acceptance angle

    n2e_theta=sqrt(  1/((sin(theta)/(n2e)).^2+(cos(theta)/(n2o)).^2) ); % n2e - for phasematching angle 29.2 deg for 800/400 for BBO
    
    dtheta=theta+d_theta;
    dn_2e_theta= n2e_theta-sqrt(1/((sin(dtheta)/(n2e)).^2+(cos(dtheta)/(n2o)).^2));
    
    double_theta=2*theta;
    delta_theta=(lambda1/(2*100))*(((n2e_theta)^3)*((n2e)^-2-(n2o)^-2)*(sin(double_theta)))^-1    % [rad per 100 um]
    
    % phasematching angle for type_I negative uniaxial crystals:
    content=(((no.^-2)-(n2o.^-2))/((n2e.^-2)-(n2o.^-2))).^0.5;
    theta_PM=asin(content)        % [rad] 
    theta_PM_deg=asind(content)   % or [deg] 
    
    % angle between propagation and optic axis 'c'  
    theta=linspace(the1,the2,thestep); 

    % angular walk-off
    ro_BBO=(1/n2e_theta)*(dn_2e_theta/d_theta)*1e3                              % [mrad]  
    
% Now compensate BBO's GVD using a calcite plate (CaCO3)   

% Sellmeier index taken from refrindexinfo site (Handbook of Optics) with VisualBasic syntax, theta=pi./2
% indices for w(lambda1=800nm) & 2w(lambda2=400nm)

n_o = sqrt( 1 +( 0.8559*lambda1^2/(lambda1^2-0.0588^2) + 0.8391*lambda1^2/(lambda1^2-0.141^2) + 0.0009*lambda1^2/(lambda1^2-0.197^2) + 0.6845*lambda1^2/(lambda1^2-7.005^2) )); 
n_2o= sqrt( 1 +( 0.8559*lambda2^2/(lambda2^2-0.0588^2) + 0.8391*lambda2^2/(lambda2^2-0.141^2) + 0.0009*lambda2^2/(lambda2^2-0.197^2) + 0.6845*lambda2^2/(lambda2^2-7.005^2) )); 
n_e = sqrt( 1 +( 1.0856*lambda1^2/(lambda1^2-0.07897^2) + 0.0988*lambda1^2/(lambda1^2-0.142^2) + 0.317*lambda1^2/(lambda1^2-11.468^2) ));
n_2e= sqrt( 1 +( 1.0856*lambda2^2/(lambda2^2-0.07897^2) + 0.0988*lambda2^2/(lambda2^2-0.142^2) + 0.317*lambda2^2/(lambda2^2-11.468^2) ));

% n calculated for n+dn(f(d_lambda))
n_o_dlambda=sqrt( 1 +( 0.8559*(lambda1+lambda1*d_lambda)^2/((lambda1+lambda1*d_lambda)^2-0.0588^2) + 0.8391*(lambda1+lambda1*d_lambda)^2/((lambda1+lambda1*d_lambda)^2-0.141^2) + 0.0009*(lambda1+lambda1*d_lambda)^2/((lambda1+lambda1*d_lambda)^2-0.197^2) + 0.6845*(lambda1+lambda1*d_lambda)^2/((lambda1+lambda1*d_lambda)^2-7.005^2) )) ;
n_2e_dlambda=sqrt( 1 +( 1.0856*(lambda2+lambda2*d_lambda)^2/((lambda2+lambda2*d_lambda)^2-0.07897^2) + 0.0988*(lambda2+lambda2*d_lambda)^2/((lambda2+lambda2*d_lambda)^2-0.142^2) + 0.317*(lambda2+lambda2*d_lambda)^2/((lambda2+lambda2*d_lambda)^2-11.468^2) ));
n_2o_dlambda=sqrt( 1 +( 0.8559*(lambda2+lambda2*d_lambda)^2/((lambda2+lambda2*d_lambda)^2-0.0588^2) + 0.8391*(lambda2+lambda2*d_lambda)^2/((lambda2+lambda2*d_lambda)^2-0.141^2) + 0.0009*(lambda2+lambda2*d_lambda)^2/((lambda2+lambda2*d_lambda)^2-0.197^2) + 0.6845*(lambda2+lambda2*d_lambda)^2/((lambda2+lambda2*d_lambda)^2-7.005^2) )); 


for k=1:length(theta)

% n = f(theta) for extraordinary (e) @ 800&400 nm
% n_e_theta=sqrt(1/((sin(theta)/(n_e))^2+(cos(theta)/(n_o))^2)); - not needed here
n_2e_theta(k)=sqrt(1/((sin(theta(k))/(n_2e)).^2+(cos(theta(k))/(n_2o)).^2)); 

% dn for 800 nm (o)
dn_o = n_o-n_o_dlambda;
dn_2e(k)= n_2e_theta(k)-sqrt(1/((sin(theta(k))/(n_2e_dlambda)).^2+(cos(theta(k))/(n_2o_dlambda)).^2));

% Calcite compensation 400nm (e), 800nm (o)
dn_dlambda_2e(k)=dn_2e(k)/(lambda2*d_lambda);                        
dn_dlambda_o=dn_o/(lambda1*d_lambda); % for ordinary 800 nm not theta dependent

% Calcite group velocity delay, [fs/100 microns]
GVD(k)=(1/3)*1.e-12*(  (n_2e_theta(k)+lambda2*dn_dlambda_2e(k))-(n_o+lambda2*dn_dlambda_o)  )*1.e15; 

% Walk-off for caclite as f(theta)
dtheta(k)=theta(k)+d_theta;
dn_2e_theta(k)=n_2e_theta(k)-sqrt(    1/((sin(dtheta(k))/(n_2e)).^2+(cos(dtheta(k))/(n_2o)).^2)   );
ro_calcite(k)=(1/n_2e_theta(k))*(dn_2e_theta(k)/d_theta);  %in [rad]
end

% Plots

% Bare calcite GVD, [fs/100um]
figure(1);      
plot(theta,GVD)

% 100um BBO + calcite GVD, [fs/100um]
figure(2);      
GVD_compensated(1:length(theta))= GVD(1:length(theta))+ 19.3; % GVD for 100 micron BBO included
plot(theta,GVD_compensated)

% Theta and GVD in phasematching condition
m=find(theta>=theta_PM,1);
ro_calcite=(1/n_2e_theta(m))*(dn_2e(m)/d_theta)
theta(m)
GVD_compensated(m)

% Angular walk-off in calcite for angle between fundamental propagation direction and optic axis
% c, [mrad]
figure(3);      
ro_calcite=ro_calcite*1e3;
plot(theta,ro_calcite)

