
% Script does Gnomomic projection of Laue hard X-ray ultrafast reflections and
% in the next step convert the pattern into sinusoids by Hough transform - to
% filter out noisy (mainly due to electron bremsstrahlung and cosmic rays) Lauegrams at 100 keV and above
% R.R. 20.Sep.2014

clear all
close all
count=1;
% projection_values=zeros(1e4,1e4);
dim=2000;
projection_values=zeros(dim,dim);
load('Lauegram');  % import here sum_frames from Process...Upd_MF_RR_processes_Si.m
%load('Lauegram_filtered'); % import filtered data

% Pre-filtering prior to sending to Gnomonic processing to get rid of hits that
% give messy background 

D=55;                                        % distance crystal-detector, [mm]
a0= 5.431;                                   % lattice constant, [Angstrom]
spot_jitter=6;                               % around Laue spot (+-spot_jitter/2) along x or y, [pixels]
Energy_jitter=1;                             % +-Energy_jitter/2, [keV]
sum_frames_filtered = zeros(80,80);          % initialization of fileterd below sum_frames for later processing in Gnomonic projection subrutine

for h = 1:10
     for k = 1:10
        for l = 1:10
             if (2*fix(h/2))-h==0
                 if (2*fix(k/2))-k==0
                     if (2*fix(l/2))-l==0
                         if h+k+l ~= 4
                             if h~=k
                                 if k~=l
                                    x=abs(((2*h*l)/(l^2-(h^2+k^2))))*D;   % [mm]                       
                                    y=abs(((2*k*l)/(l^2-(h^2+k^2))))*D;     
                                    S = sqrt(x^2+y^2);                    %  S is a distance from zero_order-Laue_spot(reflection, [mm]
                                    if x< 55*0.25                         % [mm]
                                       if y< 61*0.25                      % [mm]
                                        if S>Smin
                                        x_pixel =fix(x/0.25);             % in pixels
                                        y_pixel =fix(y/0.25); 
                                        theta = atan(S/D)/2; 
                                        d=a0/(l^2+h^2+k^2);
                                        lambda=2*d*sin(theta);            % lambda for hkl spot in Angstrom
                                        Energy = 12.4/lambda              % photon energy for hkl spot in keV
                                    
                                        Ymin=16+y_pixel-(spot_jitter/2);  % row pixel_range (for y axis)
                                        Ymax=16+y_pixel+(spot_jitter/2);
                                    
                                        Xmin=58-x_pixel-(spot_jitter/2);  % col pixel_range (for x axis)
                                        Xmax=58-x_pixel+(spot_jitter/2);
                                    
                                        E1=Energy-(Energy_jitter/2);
                                        E2=Energy+(Energy_jitter/2);
                                        for row=Ymin:Ymax 
                                            for col=Xmin:Xmax
 
                                                sum_frames_filtered(row,col)=sum_frames(row,col);
                                            end
                                        end           
                                    end
                                 end
                              end                                 
                          end
                     end
                 end
             else
              if (2*fix(k/2))-k~=0
                   if (2*fix(l/2))-l~=0
                       if h~=k
                          if k~=l
                             x=abs(((2*h*l)/(l^2-(h^2+k^2))))*D;          % [mm]                       
                             y=abs(((2*k*l)/(l^2-(h^2+k^2))))*D;     
                             S = sqrt(x^2+y^2);  
                             if x< 55*0.25
                                 if y< 61*0.25                        
                                        
                                    x_pixel =fix(x/0.25);                 % [pixels]
                                    y_pixel =fix(y/0.25); 
                                    theta = atan(S/D)/2; 
                                    d=a0/(l^2+h^2+k^2);
                                    lambda=2*d*sin(theta);                % lambda for hkl spot, [Angstrom]
                                    Energy = 12.4/lambda;                 % photon energy for hkl spot, [keV]
                                    
                                    % row pixel_range (for y axis)
                                    Ymin=16+y_pixel-(spot_jitter/2);      
                                    Ymax=16+y_pixel+(spot_jitter/2);
                                    
                                    % col pixel_range (for x axis)
                                    Xmin=58-x_pixel-(spot_jitter/2);      
                                    Xmax=58-x_pixel+(spot_jitter/2);
                                    
                                    E1=Energy-(Energy_jitter/2);
                                    E2=Energy+(Energy_jitter/2);
                                    
                                    for row=Ymin:Ymax 
                                        for col=Xmin:Xmax
                                             
                                            sum_frames_filtered(row,col)=sum_frames(row,col);
                                        end
                                    end 
                                 end                                  
                              end
                           end
                       end                         
                   end 
                end
             end
          end
       end
    end
end

% Gnomonic projection

for r=1:80                                                                % r for row, y
    for c=1:48                                                            % c for column, x
                                                                          % Laue_frame (c-58,r-16)       
    c_pixel=(0.25/2)+(c-1)*0.25;                                              % x, 0.25mm is a pixel size of Hexitec CdTe RAL camera
    r_pixel=(0.25/2)+(r-1)*0.25;                                              % y
        
    c_center_pixel=(0.25/2)+(58-1)*0.25;   
    r_center_pixel=(0.25/2)+(16-1)*0.25; 

    c_distance =abs(c_pixel-c_center_pixel);
    r_distance =abs(r_pixel-r_center_pixel);

    S= sqrt(c_distance^2 + r_distance^2);

    theta = atan(S/55)/2;                                                     % D = 55mm distance crystal-detector

    fi = atan(r_distance/c_distance);

    xG = cot(theta)* cos(fi);                                                 % = colG
    yG = cot(theta)* sin(fi);                                                 % = rowG

    xG_integer=round(xG*1e1);
    yG_integer=round(yG*1e1);

    if xG_integer == 0
       xG_integer = 1;
    end
    if yG_integer == 0
       yG_integer = 1;
    end

    IND = sub2ind(size(projection_values), yG_integer,xG_integer);
    projection_values(IND) = sum_frames_filtered(c,r);  % sum_frames(c,r); - this can be put here to look at not filtered data after Gnomonic projection
    end
end

 figure (1)
imagesc(sum_frames)
caxis([100 600])

 figure (2)
imagesc(sum_frames_filtered)
caxis([100 600])

 figure (11)
%imagesc(projection_values(1:95,1:80))    % for sum_frames
imagesc(projection_values(1:200,20:300))  % for sum_frames_filtered
caxis([100 600])

figure (111)
imagesc(projection_values)
caxis([100 600])
%caxis([E1 E2])

figure (1111)
imagesc(projection_values(66:70,22:34))    % for sum_frames
caxis([100 600])

% Start Hough transformation

Fi = linspace(0, pi, 180);
figure(1111); hold on
for i = 1:200                              % 66:70 dim
    for j = 20:300                         % 22:34 col, XG
        if projection_values(i,j)>0
           for k = 1:length(Fi)
               Hough_angle = Fi(k);
               Hough_transform(k) = j*cos(Hough_angle)+ i*sin(Hough_angle);
           end
           plot(Fi,Hough_transform);
        end
    end
end
hold off
 