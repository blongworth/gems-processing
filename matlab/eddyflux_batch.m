%Eddy flux Script
%Date created: 6/5/2013; updated 9/14/18
%Modified for batch input 11/21/24
%User inputs: 
%infile: file with time (h) ,x,y,z (m s-1),o2 (umol L-1),ph (units),pressure (dbar) 
%outfile2: output filename for fluxes (mmol m-2 h-1)and burst means
%outfile3: output filename with cumulative fluxes (mmol m-2 H-1)and o2 and ph means
%flagrotate: yes = 1, no = 0
%flagstorage: yes = 1, no = 0
%mheight: measuring height (m)
%heading: vector heading from .sen file (degrees)
%hz: sampling frequency
%timestamps: input timestamps for beginning and end of each burst
%flagwrite: yes = 1 (slow), no = 0
%Hs = waveheight (m), TKE = total kinetic energy (m s-1), 
%Ustar (m s-1), Zzero (m), EddyDiff (m2 s-1), wave velocity (m s-1), Cd =
%drag coefficient (dimensionless)

function eddyflux_batch(varargin)
  try
    if nargin < 2
      error('Need data and timestamp filenames');
    end
    infile = varargin{1}; %input filename
    timestampfile = varargin{2};
    outfile = ['outfileCPSD_', infile];
    outfile2 = ['outfile2_', infile]; %output file with fluxes
    outfile3 = ['outfile3_', infile]; %output file with cumulative fluxes
    outfile4 = ['outfile2Rot_', infile]; %output file with fluxes for rotated data
    outfile5 = ['outfile3Rot_', infile]; %output file with cumulative fluxes for rotated data
    flagwrite = 0; %write cumulative fluxes? 0=no, yes =1 (slow)
    if nargin == 3
      disp("Using planar rotation")
      flagrotate = 2;
      rotations=[-0.507, -0.477, -0.040, -0.579, -0.570, 0.305];
      planarzrot = rotations(varargin{3}); % planar rotation angle for 2nd rotation; planar z rotation
    else
      flagrotate = 1;
      planarzrot = -5.5725; % planar rotation angle for 2nd rotation; planar z rotation
    end
    mheight = 0.6; %measuring height for storage calc (m)
    heading = 0; % leave at 0  if using ENU; or heading of vector from .sen file (degrees) from column 11 if fixed frame and XYZ corrdinates
    hz = 8; %frequency of measured data
    windowsize = hz*60*5+1; %averaging window for running mean
    window = hz*60*14+1; %Hamming window that defines window size in pwelch and cpsd function
    Td = 5; % wave frequency for accumulativing at frequencies below the wave band using CPSD
    B = 4; % number of burst per hour (e.g. 15 minute bursts = 4)
    dist = 0.025; %distance between sensors (m)
    flowlag = 0.5; % (sec) time it takes water to get from tip of sensor to the sensing surface based on pump flow rate (mL/sec) and flow path volume (mL)
    %flagshift routine causes index error with test data
    flagshift = 1; %flag to timeshift data based on mean velocity
    flagwv = 0.045; %minimum velocity for estimating Cd and phi2 (for determining planarzrot)

    %input timestamps for flux calculation;
    timestamps = readmatrix(timestampfile);
    [n, m] = size(timestamps);

    disp('Importing data');
    data = importdata(infile,' ',1);

    %data columns
    time = data.data(:,1); vx = data.data(:,2); vy = data.data(:,3); vz = data.data(:,4);
    o2 = data.data(:,5); ph = data.data(:,6); pres = data.data(:,7); SNR = data.data(:,8); correlation = data.data(:,9);
    clear data;

    disp('Calculating fluxes');
    %preallocate arrays
    j1 = zeros(n,1); j2 = zeros(n,1); timebegin = zeros(n,1); timeend = zeros(n,1); timemean = zeros(n,1);
    vxmean = zeros(n,1); vymean = zeros(n,1); vzmean = zeros(n,1); vmean = zeros(n,1); Flux_cpsd_xy = zeros(n,1); Flux_cpsd_xy_low = zeros(n,1);
    o2mean = zeros(n,1); phmean = zeros(n,1); Hmean = zeros(n,1); presmean = zeros(n,1); Hs = zeros(n,1); SNRmean = zeros(n,1); correlationmean = zeros(n,1);
    flux2o2 = zeros(n,1); flux2ph = zeros(n,1); Cd = zeros(n,1); phi2 = zeros(n,1); Ustar = zeros(n,1); EddyDiff = zeros(n,1); Zzero = zeros(n,1);
    flux3o2 = zeros(n,1); flux3ph = zeros(n,1); shift = zeros(n,1); flux3o2stor = zeros(n,1); flux3phstor = zeros(n,1);
    TKE2 = zeros(n,1); TKE3 = zeros(n,1); Flux3cpsd = zeros(n,1); Flux3cpsdLow = zeros(n,1); Flux3cpsdph = zeros(n,1); Flux3cpsdphLow = zeros(n,1);
    wv1 = zeros(n,1); wv2 = zeros(n,1); wv3 = zeros(n,1); flux3vx = zeros(n,1); Flux3cpsdLowStor = zeros(n,1); Flux3cpsdphLowStor = zeros(n,1);


    %find index of burst begin and end (using timestamps)
    for i = 1:n
        jall1 = find(time > timestamps(i,1));
        j1(i) = jall1(1);
        jall2 = find(time < timestamps(i,2));
        j2(i) = jall2(length(jall2));
    end

    %convert pH
    H = (10.^(-ph))*10.^6;  %H = hydrogen ion concentration in uMol L-1
    o2storage = zeros(length(o2),1); pHstorage = zeros(length(o2),1);
    o2storage(j1(1):j2(end)) = movmean(o2(j1(1):j2(end)), hz*3600*6); 
    phstorage(j1(1):j2(end)) = movmean(H(j1(1):j2(end)), hz*3600*6);

    %calculate means for each burst
    for i = 1:n
        timemean(i) = mean(time(j1(i):j2(i)));
        vxmean(i) = mean(vx(j1(i):j2(i)));
        vymean(i) = mean(vy(j1(i):j2(i)));
        vzmean(i) = mean(vz(j1(i):j2(i)));
        o2mean(i) = mean(o2(j1(i):j2(i)));
        phmean(i) = mean(ph(j1(i):j2(i)));
        Hmean(i) = mean(H(j1(i):j2(i)));
        vmean(i) = (vxmean(i).^2+vymean(i).^2+vzmean(i).^2).^(1/2);
        presmean(i) = mean(pres(j1(i):j2(i)));
        Hs(i) = 4.*std(pres(j1(i):j2(i)));
        SNRmean(i) = mean(SNR(j1(i):j2(i)));
        correlationmean(i) = mean(correlation(j1(i):j2(i)));
    end

    %caluclate rotation angle and rotated velocities
    theta = zeros(n,1); phi = zeros(n,1); vrot1 = zeros(length(time),3); vrot2 = zeros(length(time),3);
    vxrot = zeros(n,1); vyrot = zeros(n,1); vzrot = zeros(n,1); thetad = zeros(n,1); phid = zeros(n,1);
    if flagrotate == 1
        %zrotate
        for i = 1:n
            theta(i) = atan2(vymean(i),vxmean(i));
            rot1  = [cos(theta(i)) -sin(theta(i)) 0; sin(theta(i)) cos(theta(i)) 0; 0 0 1];
            for j = j1(i):j2(i);
                vrot1(j,:) = [vx(j) vy(j) vz(j)]*rot1;
            end
            phi(i) = atan2(mean(vrot1(j1(i):j2(i),3)),mean(vrot1(j1(i):j2(i),1)));
            rot2  = [cos(phi(i)) 0 -sin(phi(i)); 0 1 0; sin(phi(i)) 0 cos(phi(i))];
            for j = j1(i):j2(i);
                vrot2(j,:) = vrot1(j,:)*rot2;
            end
            vxrot(i) = mean(vrot2(j1(i):j2(i),1));
            vyrot(i) = mean(vrot2(j1(i):j2(i),2));
            vzrot(i) = mean(vrot2(j1(i):j2(i),3));
        end
        %convert radian angle to degrees
        thetad = rad2deg(theta);
        phid = rad2deg(phi);    
    end

    if flagrotate == 2
        %zrotate
        planarzrot = deg2rad(planarzrot);
        for i = 1:n
            theta(i) = atan2(vymean(i),vxmean(i));
            rot1  = [cos(theta(i)) -sin(theta(i)) 0; sin(theta(i)) cos(theta(i)) 0; 0 0 1];
            for j = j1(i):j2(i);
                vrot1(j,:) = [vx(j) vy(j) vz(j)]*rot1;
            end
            phi(i) = planarzrot;
            rot2  = [cos(phi(i)) 0 -sin(phi(i)); 0 1 0; sin(phi(i)) 0 cos(phi(i))];
            for j = j1(i):j2(i);
                vrot2(j,:) = vrot1(j,:)*rot2;
            end
            vxrot(i) = mean(vrot2(j1(i):j2(i),1));
            vyrot(i) = mean(vrot2(j1(i):j2(i),2));
            vzrot(i) = mean(vrot2(j1(i):j2(i),3));
        end
        %convert radian angle to degrees
        thetad = rad2deg(theta);
        phid = rad2deg(phi);    
    end

    %calculate horizontal angle - this is incorrect (originally designed for
    %xyz ADV sampling and a fixed instrument heading)
    % Need to incorporate this:
    % CurrDir = 90-atan2d(V2(:,3:end),V1(:,3:end)); 
    % CurrDir(CurrDir < 0) = CurrDir(CurrDir<0)+360;
    horz = zeros(n,1);
    thetad2 = zeros(n,1);
    if flagrotate == 0;
        for i = 1:n;
            theta(i) = atan2(vymean(i), vxmean(i));
            thetad(i) = rad2deg(theta(i));
            thetad2(i) = thetad(i) + heading;
            if (thetad2(i) >= 0) && (thetad2(i) <= 360);
                horz(i) = thetad2(i);
            end
            if thetad2(i) < 0;
                horz(i) = 360+thetad2(i);
            end
            if thetad2(i) > 360;
                horz(i) = thetad2(i)-360;
            end
        end
    end
    if flagrotate == 1 || flagrotate == 2;
        for i = 1:n;
            thetad2(i) = thetad(i) + heading;
            if (thetad2(i) >= 0) && (thetad2(i) <= 360);
                horz(i) = thetad2(i);
            end
            if thetad2(i) < 0;
                horz(i) = 360+thetad2(i);
            end
            if thetad2(i) > 360;
                horz(i) = thetad2(i)-360;
            end
        end
    end

    %% Calculate Timeshift for each burst %%
    for i = 1:n;
        shift(i) = floor(dist./vmean(i).*hz+flowlag.*hz);

    end

    %% recalculate o2 timeseries with shifted data
    if flagshift == 1;
        o2old = o2;
        Hold = H;
        o2 = nan(length(time),1);
        H = nan(length(time),1);
        for i = 1:n;
            o2(j1(i):j2(i)) = o2old(j1(i)+shift(i):j2(i)+shift(i));
            H(j1(i):j2(i)) = Hold(j1(i)+shift(i):j2(i)+shift(i));
        end
    end

    % %calculate flux2 - Linear detrending
    vxprime2 = zeros(length(time),1); vyprime2 = zeros(length(time),1); vzprime2 = zeros(length(time),1); o2prime2 = zeros(length(time),1);
    phprime2 = zeros(length(time),1); o2m2 = zeros(length(time)+n,1); phm2 = zeros(length(time)+n,1); o2fit = zeros(n,2);
    phfit = zeros(n,2); vxfit = zeros(n,2); vyfit = zeros(n,2); vzfit = zeros(n,2);

    %flux 2 no rotatation
    if flagrotate == 0;
            for i = 1:n;
                o2fit(i,:) = polyfit(time(j1(i):j2(i)),o2(j1(i):j2(i)),1);
                phfit(i,:) = polyfit(time(j1(i):j2(i)),H(j1(i):j2(i)),1);
                vxfit(i,:) = polyfit(time(j1(i):j2(i)),vx(j1(i):j2(i)),1);
                vyfit(i,:) = polyfit(time(j1(i):j2(i)),vy(j1(i):j2(i)),1);            
                vzfit(i,:) = polyfit(time(j1(i):j2(i)),vz(j1(i):j2(i)),1);
                vxprime2(j1(i):j2(i)) = vx(j1(i):j2(i))-(vxfit(i,1)*time(j1(i):j2(i))+vxfit(i,2));
                vyprime2(j1(i):j2(i)) = vy(j1(i):j2(i))-(vyfit(i,1)*time(j1(i):j2(i))+vyfit(i,2));
                vzprime2(j1(i):j2(i)) = vz(j1(i):j2(i))-(vzfit(i,1)*time(j1(i):j2(i))+vzfit(i,2));
                o2prime2(j1(i):j2(i)) = o2(j1(i):j2(i))-(o2fit(i,1).*time(j1(i):j2(i))+o2fit(i,2));
                phprime2(j1(i):j2(i)) = H(j1(i):j2(i))-(phfit(i,1).*time(j1(i):j2(i))+phfit(i,2));
                flux2primeo2 = vzprime2.*o2prime2;
                flux2primeph = vzprime2.*phprime2;
                flux2o2(i) = mean(flux2primeo2(j1(i):j2(i)))*60*60;
                flux2ph(i) = mean(flux2primeph(j1(i):j2(i)))*60*60;
                TKE2(i) = ((mean(vxprime2(j1(i):j2(i)).^2))+(mean(vyprime2(j1(i):j2(i)).^2))+(mean(vzprime2(j1(i):j2(i)).^2))).^0.5;
                wv2(i)= ((mean((vx(j1(i):j2(i))- vxmean(i)).^2))+(mean((vy(j1(i):j2(i))- vymean(i)).^2))).^0.5;
            end
    end

    %Flux 2 rotated
    if flagrotate == 1 || flagrotate == 2;
          for i = 1:n;
              o2fit(i,:) = polyfit(time(j1(i):j2(i)),o2(j1(i):j2(i)),1);
              phfit(i,:) = polyfit(time(j1(i):j2(i)),H(j1(i):j2(i)),1);
              vxfit(i,:) = polyfit(time(j1(i):j2(i)),vrot2(j1(i):j2(i),1),1);
              vyfit(i,:) = polyfit(time(j1(i):j2(i)),vrot2(j1(i):j2(i),2),1);
              vzfit(i,:) = polyfit(time(j1(i):j2(i)),vrot2(j1(i):j2(i),3),1);
              vxprime2(j1(i):j2(i)) = vrot2(j1(i):j2(i),1)-(vxfit(i,1)*time(j1(i):j2(i))+vxfit(i,2));
              vyprime2(j1(i):j2(i)) = vrot2(j1(i):j2(i),2)-(vyfit(i,1)*time(j1(i):j2(i))+vyfit(i,2));
              vzprime2(j1(i):j2(i)) = vrot2(j1(i):j2(i),3)-(vzfit(i,1)*time(j1(i):j2(i))+vzfit(i,2));
              o2prime2(j1(i):j2(i)) = o2(j1(i):j2(i))-(o2fit(i,1).*time(j1(i):j2(i))+o2fit(i,2));
              phprime2(j1(i):j2(i)) = H(j1(i):j2(i))-(phfit(i,1).*time(j1(i):j2(i))+phfit(i,2));
              flux2primeo2 = vzprime2.*o2prime2;
              flux2primeph = vzprime2.*phprime2;
              flux2o2(i) = mean(flux2primeo2(j1(i):j2(i)))*60*60;
              flux2ph(i) = mean(flux2primeph(j1(i):j2(i)))*60*60;
              TKE2(i) = ((mean(vxprime2(j1(i):j2(i)).^2))+(mean(vyprime2(j1(i):j2(i)).^2))+(mean(vzprime2(j1(i):j2(i)).^2))).^0.5;
              wv2(i)= ((mean((vx(j1(i):j2(i))- vxmean(i)).^2))+(mean((vy(j1(i):j2(i))- vymean(i)).^2))).^0.5;
          end
    end


    %calculate flux3
    vxprime3 = zeros(length(time),1); vyprime3 = zeros(length(time),1); vzprime3 = zeros(length(time),1);
    o2prime3 = zeros(length(time),1); phprime3 = zeros(length(time),1);

    if flagrotate == 0; 
          for i = 1:n; 
                vxprime3(j1(i):j2(i)) = vx(j1(i):j2(i))-smooth(vx(j1(i):j2(i)),windowsize);
                vyprime3(j1(i):j2(i)) = vy(j1(i):j2(i))-smooth(vy(j1(i):j2(i)),windowsize);
                vzprime3(j1(i):j2(i)) = vz(j1(i):j2(i))-smooth(vz(j1(i):j2(i)),windowsize);
                o2prime3(j1(i):j2(i)) = o2(j1(i):j2(i))-smooth(o2(j1(i):j2(i)),windowsize);
                phprime3(j1(i):j2(i)) = (H(j1(i):j2(i))-smooth(H(j1(i):j2(i)),windowsize));
                flux3primeo2 = vzprime3.*o2prime3;
                flux3primeph = vzprime3.*phprime3;
                flux3o2(i) = mean(flux3primeo2(j1(i):j2(i)))*60*60;
                flux3ph(i) = mean(flux3primeph(j1(i):j2(i)))*60*60;
                flux3o2stor(i) = flux3o2(i)-(mean(o2storage(j1(i):(j1(i)+(hz*60))))-mean(o2storage((j2(i)-(hz*60)):j2(i))))*B*mheight;
                flux3phstor(i) = flux3ph(i)-(mean(phstorage(j1(i):(j1(i)+(hz*60))))-mean(phstorage((j2(i)-(hz*60)):j2(i))))*B*mheight;
                TKE3(i) = ((mean(vxprime3(j1(i):j2(i)).^2))+(mean(vyprime3(j1(i):j2(i)).^2))+(mean(vzprime3(j1(i):j2(i)).^2))).^0.5;
                wv3(i)= ((mean((vx(j1(i):j2(i))- vxmean(i)).^2))+(mean((vy(j1(i):j2(i))- vymean(i)).^2))).^0.5;
                            
                [Cxy,freqxy]=cpsd(vxprime3(j1(i):j2(i)),vzprime3(j1(i):j2(i)),window,[],[],hz);%compute co-spectrum of uprime wprime
                df=diff(freqxy(1:2));%frequency interval
                Flux_cpsd_xy(i)=sum(real(Cxy)*df);%integrate co-spectrum (real part only)
                filt=find(freqxy<1/Td);%find frequencies lower than a wave with 3 second period
                Flux_cpsd_xy_low(i)=sum(real(Cxy(filt))*df);%integrate only low frequency co-spectrum
                Ustar(i) = ((Flux_cpsd_xy_low(i)^2)^0.5)^0.5;
                Zzero(i) = mheight./(exp((vmean(i).*0.41)./Ustar(i))); 
                EddyDiff(i) = 0.41.*Ustar(i).*mheight;
                if  vmean(i) < flagwv % removes Cd values that are contaminated due to waves
                    Cd(i) = NaN;
                    phi2(i) = NaN;
                else
                    Cd(i) = ((Flux_cpsd_xy_low(i)^2)^0.5)/(vmean(i).^2);
                    phi2(i) = phid(i);
                end    
                                  
                [Co2vz,freqo2vz]=cpsd(o2prime3(j1(i):j2(i)),vzprime3(j1(i):j2(i)),window,[],[],hz);%compute co-spectrum of uprime wprime
                df=diff(freqo2vz(1:2));%frequency interval
                Flux3cpsd(i)=(sum(real(Co2vz)*df))*3600;%integrate co-spectrum (real part only)
                low=find(freqo2vz<1/Td);%find frequencies lower than a wave with 3 second period
                Flux3cpsdLow(i)=(sum(real(Co2vz(low))*df))*3600;%integrate only low frequency co-spectrum
                Flux3cpsdLowStor(i) = Flux3cpsdLow(i)-(mean(o2storage(j1(i):(j1(i)+(hz*60))))-mean(o2storage((j2(i)-(hz*60)):j2(i))))*B*mheight;
                
                [Cphvz,freqphvz]=cpsd(phprime3(j1(i):j2(i)),vzprime3(j1(i):j2(i)),window,[],[],hz);%compute co-spectrum of uprime wprime
                df=diff(freqphvz(1:2));%frequency interval
                Flux3cpsdph(i)=(sum(real(Cphvz)*df))*3600;%integrate co-spectrum (real part only)
                low=find(freqphvz<1/Td);%find frequencies lower than a wave with 3 second period
                Flux3cpsdphLow(i)=(sum(real(Cphvz(low))*df))*3600;%integrate only low frequency co-spectrum    
                Flux3cpsdphLowStor(i) = Flux3cpsdphLow(i)-(mean(phstorage(j1(i):(j1(i)+(hz*60))))-mean(phstorage((j2(i)-(hz*60)):j2(i))))*B*mheight;
            end
    end

    if flagrotate == 1 || flagrotate == 2;
            for i = 1:n;
                vxprime3(j1(i):j2(i)) = vrot2(j1(i):j2(i),1)-smooth(vrot2(j1(i):j2(i),1),windowsize);
                vyprime3(j1(i):j2(i)) = vrot2(j1(i):j2(i),2)-smooth(vrot2(j1(i):j2(i),2),windowsize);
                vzprime3(j1(i):j2(i)) = vrot2(j1(i):j2(i),3)-smooth(vrot2(j1(i):j2(i),3),windowsize);
                o2prime3(j1(i):j2(i)) = o2(j1(i):j2(i))-smooth(o2(j1(i):j2(i)),windowsize);
                phprime3(j1(i):j2(i)) = (H(j1(i):j2(i))-smooth(H(j1(i):j2(i)),windowsize));%.*-5583-17.11668;%%%%%%%%%%%%%%
                flux3primeo2 = vzprime3.*o2prime3;
                flux3primeph = vzprime3.*phprime3;
                flux3o2(i) = mean(flux3primeo2(j1(i):j2(i)))*60*60;
                flux3ph(i) = mean(flux3primeph(j1(i):j2(i)))*60*60;
                flux3o2stor(i) = flux3o2(i)-(mean(o2storage(j1(i):(j1(i)+(hz*60))))-mean(o2storage((j2(i)-(hz*60)):j2(i))))*B*mheight;
                flux3phstor(i) = flux3ph(i)-(mean(phstorage(j1(i):(j1(i)+(hz*60))))-mean(phstorage((j2(i)-(hz*60)):j2(i))))*B*mheight;
                TKE3(i) = ((mean(vxprime3(j1(i):j2(i)).^2))+(mean(vyprime3(j1(i):j2(i)).^2))+(mean(vzprime3(j1(i):j2(i)).^2))).^0.5;
                wv3(i)= ((mean((vx(j1(i):j2(i))- vxmean(i)).^2))+(mean((vy(j1(i):j2(i))- vymean(i)).^2))).^0.5;
                            
                [Cxy,freqxy]=cpsd(vxprime3(j1(i):j2(i)),vzprime3(j1(i):j2(i)),window,[],[],hz);%compute co-spectrum of uprime wprime
                df=diff(freqxy(1:2));%frequency interval
                Flux_cpsd_xy(i)=sum(real(Cxy)*df);%integrate co-spectrum (real part only)
                filt=find(freqxy<1/Td);%find frequencies lower than a wave with 3 second period
                Flux_cpsd_xy_low(i)=sum(real(Cxy(filt))*df);%integrate only low frequency co-spectrum
                Ustar(i) = ((Flux_cpsd_xy_low(i)^2)^0.5)^0.5;
                Zzero(i) = mheight./(exp((vmean(i).*0.41)./Ustar(i))); 
                EddyDiff(i) = 0.41.*Ustar(i).*mheight;
                if  vmean(i) < flagwv % removes Cd values that are contaminated due to waves
                    Cd(i) = NaN;
                    phi2(i) = NaN;
                else
                    Cd(i) = ((Flux_cpsd_xy_low(i)^2)^0.5)/(vmean(i).^2);
                    phi2(i) = phid(i);
                end               
                
                                  
                [Co2vz,freqo2vz]=cpsd(o2prime3(j1(i):j2(i)),vzprime3(j1(i):j2(i)),window,[],[],hz);%compute co-spectrum of uprime wprime
                df=diff(freqo2vz(1:2));%frequency interval
                Flux3cpsd(i)=(sum(real(Co2vz)*df))*3600;%integrate co-spectrum (real part only)
                low=find(freqo2vz<1/Td);%find frequencies lower than a wave with 3 second period
                Flux3cpsdLow(i)=(sum(real(Co2vz(low))*df))*3600;%integrate only low frequency co-spectrum
                Flux3cpsdLowStor(i) = Flux3cpsdLow(i)-(mean(o2storage(j1(i):(j1(i)+(hz*60))))-mean(o2storage((j2(i)-(hz*60)):j2(i))))*B*mheight;
                
                [Cphvz,freqphvz]=cpsd(phprime3(j1(i):j2(i)),vzprime3(j1(i):j2(i)),window,[],[],hz);%compute co-spectrum of uprime wprime
                df=diff(freqphvz(1:2));%frequency interval
                Flux3cpsdph(i)=(sum(real(Cphvz)*df))*3600;%integrate co-spectrum (real part only)
                low=find(freqphvz<1/Td);%find frequencies lower than a wave with 3 second period
                Flux3cpsdphLow(i)=(sum(real(Cphvz(low))*df))*3600;%integrate only low frequency co-spectrum    
                Flux3cpsdphLowStor(i) = Flux3cpsdphLow(i)-(mean(phstorage(j1(i):(j1(i)+(hz*60))))-mean(phstorage((j2(i)-(hz*60)):j2(i))))*B*mheight;
            end
    end

    %calculate cumulative fluxes
    disp('Caluclating cumulative fluxes');
    cum2o2 = zeros(length(time)+n,1); cum2ph = zeros(length(time)+n,1); cum3o2 = zeros(length(time)+n,1); cum3ph = zeros(length(time)+n,1);time1 = zeros(length(time)+n,1);

    if flagwrite == 1;
        k = 1;
        for i = 1:n;
            for j = j1(i):j2(i);
                if j == j1(i);
                    cum2o2(k) = flux2primeo2(j);
                    cum2ph(k) = flux2primeph(j);
                    cum3o2(k) = flux3primeo2(j);
                    cum3ph(k) = flux3primeph(j);
                    time1(k) = time(j);
                    k = k+1;
                else
                time1(k) = time(j);
                sum2o2 = cum2o2(k-1)+flux2primeo2(j);
                sum2ph = cum2ph(k-1)+flux2primeph(j);
                cum2o2(k) = sum2o2;
                cum2ph(k) = sum2ph;
                sum3o2 = cum3o2(k-1)+flux3primeo2(j);
                sum3ph = cum3ph(k-1)+flux3primeph(j);
                cum3o2(k) = sum3o2;
                cum3ph(k) = sum3ph;
                k = k+1;
                end
                if j == j2(i);
                    cum2o2(k) = NaN;
                    cum2ph(k) = NaN;
                    cum3o2(k) = NaN;
                    cum3ph(k) = NaN;
                    time1(k) = NaN;
                    k = k+1;
                end
            end
        end
    end
    %Convert uMol m L-1 s-1 at sampling frequency(hz) to mmol m-2 h-1
    cum2o2b = cum2o2./hz.*B;
    cum2phb = cum2ph./hz.*B;
    cum3o2b = cum3o2./hz.*B;
    cum3phb = cum3ph./hz.*B;

    %write data to files
    disp('Writing data to files');
    y1 = [timemean, vmean, o2mean, phmean, presmean, Hs, TKE2, TKE3, flux2o2, flux3o2, flux3o2stor, Flux3cpsdLow, Flux3cpsdLowStor, flux2ph, flux3ph, flux3phstor, Flux3cpsdphLow, Flux3cpsdphLowStor, SNRmean, correlationmean, vxmean, vymean, vzmean, vxrot, vyrot, vzrot, thetad, phid, phi2, horz, wv2, wv3, Ustar, Cd, Zzero, EddyDiff, Flux_cpsd_xy_low];
    y1 = y1';
    if flagrotate == 0;
            fid1 = fopen(outfile2, 'w');
    end
    if flagrotate == 1 || flagrotate == 2;
            fid1 = fopen(outfile4, 'w');
    end

    fprintf(fid1, 'timemean vmean o2mean phmean presmean Hs TKE2 TKE3 flux2o2 flux3o2 flux3o2stor Flux3cpsdLow Flux3cpsdLowStor flux2ph flux3ph flux3phstor Flux3cpsdphLow Flux3cpsdphLowStor SNRmean correlationmean vxmean vymean vzmean vxrot vyrot vzrot theta phid phi2 horz wv2 wv3 Ustar Cd Zzero EddyDiff Flux_cpsd_xy_low \r\n');
    fprintf(fid1,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %e %e %e %e %e\r\n', y1);
    fclose(fid1);

    if flagwrite == 1;
        y2 = [time1, cum2o2b, cum2phb, cum3o2b, cum3phb];
        k = find(time1,1,'last');
        y2 = y2(1:k,:);
        y2 = y2';
        if flagrotate == 0;
            fid2 = fopen(outfile3, 'w');
        end
        if flagrotate == 1 || flagrotate == 2;
            fid2 = fopen(outfile5, 'w');
        end

        fprintf(fid2, 'time1 cum2o2b cum2phb cum3o2b cum3phb \r\n');
        fprintf(fid2,'%f %f %E %f %E \r\n', y2);
        fclose(fid2);
    end
    disp('Done');

  catch ME 
    fprintf('Error: %s\n', ME.message);
    exit(1);
end
