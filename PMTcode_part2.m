
% PHASE MATCHING IN TIME (PMT)code 
% Script calculates the temporally resolved harmonic signal for neon 
% in 3D as a function of intensity Iq(r,z,t) assuming axial symmetry and phasematching terms
% R.R. 25.Oct.2010

clear all
close all

tic; 
input_data;

% q=31:2:89;                                % harmonic order range
 
% Initialization of big arrays
            
deltak=zeros(rstep,zstep,tstep);            % resultant wavevector mismatch   
dipole_moments=zeros(rstep,zstep,tstep);    % short trajectory 
atomicphaseS=zeros(rstep,zstep,tstep);      % short trajectory dipole phase
atomicphaseL=zeros(rstep,zstep,tstep);      % long trajectory dipole phase
deltakdisp=zeros(rstep,zstep,tstep);        % plasma dispersion phases
deltakgouy=zeros(zstep);                    % geometrical phase
deltakneutr=0;                              % neutral gas dispersion phase
wp=zeros(rstep,zstep,tstep);                % plasma frequency

load('degree_ionization_en4e14_511G');
load('intensities_en4e14_511G');
load('ionization_probability_en4e14_511G');
load('neutral_refractive_indices');     


            Na=pressure*Avo/R/T;                                           % Density of Ne atoms, m^-3
            density=Na*atommass/Avo*1.e-6;                                 % Neon density, g/cm^3
            Ne=degree.*Na;                                                 % Density of electrons, m^-3
         
            n1=1+Na*alpha/2/epsilon0;                                      % Refractive index for fundamental
            wp=sqrt(Ne.*e^2/m/epsilon0);                                   % Plasma frequency 

            load absorption_coefficients.m                                 % takes and calculates the absorption K coeff.
            mi_crosssections=absorption_coefficients;
            mis=mi_crosssections(:,3);                                     % crosssections, cm^2/g, for q=31:2:89
            
            % Absorption factors
            absorption=mis.*density/100/2;                                 % Kappaq coefficients, (1/m)
            sigmas=mis.*atommass/Avo;                                      % cross-sections, cm^2/atom
            
            i=sqrt(-1);
            
% Nested loop for vector mismatches calculation

for qq=1:length(q)          
  for tt=1:tstep
    for zz=1:zstep
       for rr=1:rstep
                                                      
          % Single atom dipole-moments
          
            dipole_moments(rr,zz,tt)=dipole_moments_temp(rr,zz,tt).*sqrt(q(qq)).*sigmas(qq); 

          % Single atom phase as f(Io,alfa(q),z,r,t) for short & long traj.
          
            ne=load('ne.txt');
            alphaS_q=ne(q(qq),5);
            alphaL_q=ne(q(qq),1);                                         
            
            atomicphaseS(rr,zz,tt)=alphaS_q.*Int(rr,zz,tt);
            atomicphaseL(rr,zz,tt)=alphaL_q.*Int(rr,zz,tt);

          % Plasma dispersion as f(q,wp(Io,z,r),t)
          
            np1=sqrt(1-(wp(rr,zz,tt)./w).^2);  
            npq=sqrt(1-(wp(rr,zz,tt)./q(qq)/w).^2);
            deltakdisp(rr,zz,tt)=wp(rr,zz,tt)^2*((1/(2*c*q(qq).*w))-(q(qq)./(2*w*c)));
            kqdisp(rr,zz,tt)=-wp(rr,zz,tt)^2*(1/(2*c*q(qq).*w));
            
          % Neutral dispersion as f(q)       
   
            deltakneutr=(n1-real(nq(qq))).*q(qq).*w/c;
            kqneutr=real((nq(qq)-1)).*q(qq).*w/c;
            
          % Gouy phase = Gaussian deltak as f(q,z)
            deltakgouy(zz)=(1-q(qq)).*(1/Zr*(1/(1+(z(zz)/Zr)^2)));
            kqgouy(zz)= (-1/Zr*(1/(1+(z(zz)/Zr)^2)))  +  ( ((q(qq).*w*(r(rr)^2))/(2*c))*((-z(zz)+Zr)/((z(zz)+Zr)^3)) ) ;
            
          % Resultant deltaks
            
            deltak(rr,zz,tt)=deltakneutr+deltakdisp(rr,zz,tt)+deltakgouy(zz);       % 4Apr.2012 extracted deltakgouy(zz)
            kq(rr,zz,tt)=kqgouy(zz)+kqneutr+kqdisp(rr,zz,tt);                       % 9Apr. 2012  extracted kqgouy(zz)+ inserted kggoy(zz) again on 10Apr
            
            totaldeltak(rr,zz,tt)=deltak(rr,zz,tt)-i*absorption(qq); 
            
        end
    end
 end

 save('deltak_en4e14_511G','deltak');
 save('kqgouy_en4e14_511G','kqgouy');
 save('deltakgouy_en4e14_511G','deltakgouy');
 save('kqneutr_en4e14_511G','kqneutr');
 save('deltakneutr_en4e14_511G','deltakneutr');
 save('kqdisp_en4e14_511G','kqdisp'); 
 save('deltakdisp_en4e14_511G','deltakdisp');
 
 % Now main loop
                                                      
 for ttt=1:tstep
    for rrr=1:rstep
        for zzz=1:zstep
            path=zre(1:zzz);
                                                                          
            if  (length(path)>1)
                phase(rrr,zzz)=trapz(path,i*totaldeltak(rrr,1:length(path),ttt))*1.; 
            else
                phase(rrr,zzz)=0.;
            end
            
                phase_S(rrr,zzz)=phase(rrr,zzz)+i*atomicphaseS(rrr,zzz,ttt);                       
                phase_L(rrr,zzz)=phase(rrr,zzz)+i*atomicphaseL(rrr,zzz,ttt);                      
                
                Eqtemp_S(rrr,zzz)=exp(phase_S(rrr,zzz)).*dipole_moments(rrr,zzz,ttt)*(1-degree(rrr,zzz,ttt))*Na*pi*dr^2*dz;
                Eqtemp_L(rrr,zzz)=exp(phase_L(rrr,zzz)).*dipole_moments(rrr,zzz,ttt)*(1-degree(rrr,zzz,ttt))*Na*pi*dr^2*dz;
            
            %%%
            ...place for partial sums if there is a future need
            %%%
        end
            
    end
            % In separate script calculate Hankel-Furier transforms for
            % short & long trajectory of electrons
            
            [EqtempT_S]=HT_function(Eqtemp_S,Nfft);
            [EqtempT_L]=HT_function(Eqtemp_L,Nfft);
    
        for rrrr=1:56 
            Eq_S(rrrr,ttt)=trapz(zre,EqtempT_S(rrrr+56,1:length(zre)));
            Eq_L(rrrr,ttt)=trapz(zre,EqtempT_L(rrrr+56,1:length(zre)));   
                                       
            abscorrection=trapz(zre,imag(totaldeltak(rrrr,1:length(zre),ttt)));
            phasecorrection=trapz(zre,-real(kq(rrrr,1:length(zre),ttt)));
            
            Eqtot_S(qq,ttt,rrrr)=exp(abscorrection).*exp(i*phasecorrection).*Eq_S(rrrr,ttt);
            Eqtot_L(qq,ttt,rrrr)=exp(abscorrection).*exp(i*phasecorrection).*Eq_L(rrrr,ttt);
        end
 end
    
    % Output harmonic field for further coherent combining of both short % long
    % trajectory contributions
    save('harmonic_field_S_en4e14_511G','Eqtot_S');
    save('harmonic_field_L_en4e14_511G','Eqtot_L');
    qq
end

% Save last iteration's fields and phases to review for bugs

save('dipole_moments_en4e14_511G','dipole_moments');
save('atomicphaseS_en4e14_511G','atomicphaseS');
save('atomicphaseL_en4e14_511G','atomicphaseL');
save('kq_en4e14_511G','kq');
save('totaldeltak_en4e14_511G','totaldeltak');
save('phase_S_en4e14_511G','phase_S');
save('phase_L_en4e14_511G','phase_L');
save('Eqtemp_S_en4e14_511G','Eqtemp_S');
save('Eqtemp_L_en4e14_511G','Eqtemp_L');
save('EqtempT_S_en4e14_511G','EqtempT_S');
save('EqtempT_L_en4e14_511G','EqtempT_L');
save('Eq_S_en4e14_511G','Eq_S');
save('Eq_L_en4e14_511G','Eq_L');
save('abscorrection_en4e14_511G','abscorrection');
save('phasecorrection_en4e14_511G','phasecorrection');

disp(strcat('Execution time:',num2str(toc/3600),'h'));

