
% Script calculates 3D intensities as f(r,z,t) and ionization
% degrees(gives input data for the second and main part of the PMT code)
% R.R. 23.Sep.2010

% q harmonic order
% r harmonic field radius, m
% z gas cell length, m

clear all
close all

tic; 
input_data;

% Initialization of big arrays

I=zeros(rstep,zstep);  
Int=zeros(rstep,zstep,tstep); 
degree=zeros(rstep,zstep,tstep); 
dipole_moments_temp=zeros(rstep,zstep,tstep);

% Loop to calculate laser intensities  in 3D and ionization degrees for the time
% steps within a laser pulse
 
    for rr=1:length(r)
       for zz=1:length(z)
  
         I(rr,zz)=Io*((1/(sqrt(1+(z(zz)/Zr)^2)))^2)*exp((-2*(r(rr)^2))/((Zr*lambda/pi)*(1+(z(zz)/Zr)^2)));   % Gaussian laser Int. peaked in time (exp(t)=1)
            II=I(rr,zz);
            [p,tmin,tmax,ionization,tcentral,t1,t2,dipole_moments_temps]=ADK_degrees(II,tau); 
            order=length(p);                                                                                 % gives no of factors p in polynomial fit (in fact: poly.order+1)
            
            t=linspace(tmin,tmax,tstep);
            tre=linspace(-(tcentral-tmin),tmax-tcentral,tstep);
            
            Int(rr,zz,1:length(t))=II*exp((-2.78*(tre(1:length(tre)).^2))/((tau*e-15)^2));
            degrees=polyval(p,t);
            
            firstcut=find(t>=t1,1);       
            degrees(1:firstcut)=max(0,degrees(firstcut));           
            secondcut=find(t>=t2,1);          
            for i=firstcut+1:secondcut-1
                if degrees(i)<0
                degrees(i)=degrees(firstcut);
                end
            end
            
            maxi=max(degrees(1:secondcut));
            degrees(secondcut:end)=max(maxi,degrees(secondcut));
            
            dipole_moments_temp(rr,zz,1:length(t))=dipole_moments_temps(1:end);
            
            degree(rr,zz,1:length(t))=degrees(1:end);      
       end 
       rr
    end
    
    [B]=beta_MetteG(Io);
    degree=degree*B;
save('degree_ionization_en4e14_511G','degree');
save('intensities_en4e14_511G','Int');
save('ionization_probability_en4e14_511G','dipole_moments_temp');
save('tmin','tmin');
save('tmax','tmax');

disp(strcat('Execution time:',num2str(toc/3600),'h'));

