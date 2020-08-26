
function [p,tmin,tmax,ionization,tcentral,t1,t2,dipole_moments_temps]=ADK_degrees(I,tau)

%clear all   % needed when this script is tested separately and not used as a function            
%tau=30.;    %
%I=4E13;     %

clear dipole_moments_temp

input_data; 

tau=tau*1e15;                                % FWHM pulse durationin, fs
I=I*1e-4;                                    % Partial amplitude Intensity(r,z), W/cm2

%Eion = 0.58065;                             % Ionization Potential for Argon, a.u.
%Eion2=27.6/27.21;                           % 2nd Ionization Potential for Argon, a.u.
Eion = 21.5646/27.21;                        % Ionization Potential for Neon, a.u.
Eion2=40.1/27.21;                            % 2nd Ionization Potential for Neon, a.u.
QN = 1./sqrt(2*Eion);
QN2 = 2./sqrt(2*Eion2);
%Ain = 2.0E14;                               % Peak intensity, W/cm2
Ain = I;                                     % Local peak intensity I=f(R,z), W/cm2
Ain = sqrt(Ain/3.51E16);                     % Electric field amplitude
Tpul= tau;                                   % Pulse duration (FWHM), fs
Tpul = Tpul*1E-15*4.1322E16;
Tf = Tpul*0.84932;
wo = 0.05695;
TT = 2.*pi/wo;
Nx = 10000;                                  % Number of time steps
HT = 4.*Tf/Nx;

%### ADK kernel for atoms taken from the web###
for i=1:Nx
    T = i*HT;
    E(i)=Ain*exp(-(T-2.0*Tf)^2/Tf^2)*cos(wo*(T-2*Tf));
    Int(i)=(Ain*exp(-(T-2.0*Tf)^2/Tf^2))^2;
end
 
N(1)=1.;
N1(i)=0.;
N2(1)=0.;
T(1)=0.;
Rint=0.;
for i=1:Nx-1
           T(i+1) = i*HT;
           F = abs(E(i))+1.E-30;
           C2=2^(2*QN)/(QN*gamma(QN+1)*gamma(QN));
           C22=2^(2*QN2)/(QN2*gamma(QN2+1)*gamma(QN2));
           Rate1(i)=Eion*C2*(4.0*Eion*sqrt(2.0*Eion)/F)^(2.0*QN-1.0)*...
           exp(-4.0*Eion*sqrt(2.0*Eion)/(3.0*F));
          %Rate2(i)=Eion2*C22*(4.0*Eion2*sqrt(2.0*Eion2)/F)^(2.0*QN2-1.0)*exp(-4.0*Eion2*sqrt(2.0*Eion2)/(3.0*F));
           Rint=Rint+Rate1(i)*HT;
           N(i+1)=exp(-Rint);                % Number of Neutrals,*(1-N1(i));
           N1(i+1)=1.-N(i);                  % Number of Single Ions
          %N2(i+1)= N2(i)+Rate2(i)*N1(i)*HT; % Number of Double Ions
           Rate1(i)=Rate1(i)*N(i);
          %Rate2(i)=Rate2(i)*N1(i);
end
  
 Rate1(10000)=0.;     % for interpolation to equal vector lengths
%Neutral=1-N;
 Ion1=N1;
%Ion2=N2;
 degree=N1;

%T1=(T-2*Tf)/TT*2.67; % Convert to fs
 T1=T/TT*2.67;        % Convert to fs

%figure(1),plot(T1,degree);
%%%axis([0  4*tau  0  1]);
%figure(2),plot(1:10000,degree);
%figure(3),plot(1:10000,Int);
%figure(4),plot(1:10000,E,1:10000,Int);
%figure(5),plot(1:10000,Ion2);
%figure(6),plot(1:10000,N);

%########END of ADK ############

% Find smallest index within Nx where ionization starts

 j=find(degree>=0.001*max(degree),1);
 plus=find(degree>=0.,1);
  if plus==1          % let it count here from degree=0, later to be limited to tmin
     j=plus;
  end
 t1=T(j)/TT*2.67;
 t1=int8(t1);

% Find the smallest index within Nx where ionization gets flat

ionization=max(degree);
i=10000;
while (degree(end)-degree(i))<0.0001*max(degree) % or i=1%
      i=i-1;
end
t2=T(i)/TT*2.67;
t2=int8(t2);   

% Find polynomial fit of the order n (n=1 - linear interpolation with least squares)

n=10;                       % order of the polynomial fit

a=(T(j)/TT*2.67);
b=(((T(i)/TT*2.67)-(T(j)/TT*2.67))/(i-j));
c=(T(i)/TT*2.67);
X=a:b:c;

% figure(7);
p = polyfit(X, degree(j:i), n);
% plot(X,degree(j:i),X,polyval(p,X), 'r-');

%s=size(p);                 % gives two numbers s(1) as no. of rows and s(2) as no of columns
%N=s(2);                    % N gives no. of factors in polynomial (N=n+1)

tcentral=T(5000)/TT*2.67;   % absolute center of the input IR pulse in [fs]

% Calculate time vector for which all (ion. degree and final harmonic fild)
% would be considered(calculated)

switch tau   
       case 30
       tmin=tcentral-(tau/2)-10;    % in fs 
       tmax=tcentral+(tau/2)+10;    % in fs 
       case 50
       tmin=tcentral-(tau/2)-20;    % in fs 
       tmax=tcentral+(tau/2)+10;    % in fs 
       case 70  
       tmin=tcentral-(tau/2)-30;    % in fs 
       tmax=tcentral+(tau/2)+10;    % in fs            
end   

kmin=find(T>=tmin*TT/2.67,1);
kmax=find(T>=tmax*TT/2.67,1);
exline=linspace(T(kmin),T(kmax),tstep); % remember change it to (t1-t2)/tstep and zeros as dipole moments out of t1-t2 complementary to tmin-tmax
dipole_moments_temps=interp1(T,Rate1,exline);
%real_exline=linspace(T(kmin)/TT*2.67,T(kmax)/TT*2.67,tstep);
%figure(13),plot(real_exline,dipole_moments_temp);
%axis([tmin  tmax  0  max(dipole_moments_temp)]);

% New part below for the envelope in time calculations only!

if  sum(dipole_moments_temps)>0

step=(tmax-tmin)/tstep;
steps_between_ionization_peaks=pi/(w*step)*1e15;
interval=ceil(steps_between_ionization_peaks);   
     
Maximum=max(dipole_moments_temps);
J=find(dipole_moments_temps==Maximum,1);

AA(J)=Maximum;

k=J+fix(interval/2);
while k<(tstep-interval),
    
      k=k+interval;
      pe=(k-interval);
      value=max(dipole_moments_temps(pe:k));
      
      l=find(dipole_moments_temps(pe:k)==value,1);
      ll=pe+(l-1);
      AA(ll)=value;
      %BB(ll)=ll;
end
AA(tstep)=dipole_moments_temps(tstep);
%BB(tstep)=tstep

k=J-fix(interval/2);
while k>(1+interval),
    
      k=k-interval;
      py=k+interval;
      value=max(dipole_moments_temps(k:py));
      
      bla=dipole_moments_temps(k:py);
      l=find(bla==value,1);
      ll=k+(l-1);
      AA(ll)=value;
      %BB(ll)=ll;
end
AA(1)=dipole_moments_temps(1);
%BB(1)=1;
%figure(15),plot(AA);

No=1;
for m=1:tstep
    if AA(m)> 0
       data(No)=AA(m);
       ex(No)=exline(m);
       No=No+1;     
    else
        
    end
   
end

if (No > 2)     
    dipole_moments_temps=interp1(ex, data, exline);   
end
       
end

TF = isnan(dipole_moments_temps);
Index=find(TF);
dipole_moments_temps(Index)=0;

 

