%{
MIT License

Copyright (c) [2016] [Mallory Ann Jensen, jensenma@alum.mit.edu]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
%}

function compare_QSSPC_lifetime(deltanQSSPC,tauQSSPC,deltanLifetime,tauLifetime,filenameSintonCircle,saveStart)
%compare_QSSPC_lifetime(deltanQSSPC,tauQSSPC,deltanLifetime,tauLifetime,filenameSintonCircle).
%This function takes as input a single QSSPC curve and a single pair of
%spatially resolved lifetime/injection level maps. The additional input is
%the filename for the Sinton circle - this script currently assumes the
%maps are generated from PC-PL. The units of the lifetime should be in
%microseconds in both cases. A figure is generated comparing the
%injection-dependent lifetime curve the both the arithmetic and harmonic
%means around the measured area. 

%We need to get the spatially-resolved lifetime averages for comparison

%If PC-PL, use the image of the Sinton circle to choose the averaging area
figure; 
%Define delimiter and header in order to read txt files
delimiterIn = ',';
headerlinesIn = 0;
PLsensor = importdata(filenameSintonCircle,delimiterIn,headerlinesIn);
figure;
imagesc(PLsensor);
axis('image');

circletest = 0;
while circletest == 0
    
    PLsensor_circletest = PLsensor;
    
    disp('Click the center of the circle ...');
    [x,y] = ginput(1);
    disp('Click the edge of the circle ...');
    [xr,yr] = ginput(1);
    radius = sqrt((x-xr)^2 + (y-yr)^2);
    
    %PL maps have 1024 x 1024 pixels
    %set all pixels outside of the sensor region to NaN
    for i = 1:1024
        for j = 1:1024
            if abs(sqrt(abs(i-y)^2+abs(j-x)^2))>radius
                PLsensor_circletest(i,j) = NaN;
            end
        end
    end
    
    %Plot the new map with the sensor region only
    figure;
    imagesc(PLsensor_circletest);
    axis('image');
    
    promptcorrect = 'Is the mask acceptable?';
    str = input(promptcorrect,'s');
    if isempty(str)
        str = 'Y';
    end
    
    if str == 'Y'
        circletest = 1;
    end
    
end

lifetime = figure;
h(1)=loglog(deltanQSSPC,tauQSSPC,'LineWidth',2); 
hold all;
count = 2; 

for i = 1:length(deltanLifetime)
    deltan_now = deltanLifetime{i};
    tau_now = tauLifetime{i};
    %Mask the spatially-resolved lifetime and injection maps
    [m,n] = size(deltan_now); 
    for i = 1:m
        for j = 1:n
            if abs(sqrt(abs(i-y)^2+abs(j-x)^2))>radius
                deltan_now(i,j) = NaN;
                tau_now(i,j) = NaN;
            end
        end
    end

    %Plot the results 
    figure;
    imagesc(deltan_now); 
    axis('image');
    title('Masked injection level');

    figure;
    imagesc(tau_now);
    axis('image');
    title('Masked lifetime'); 

    [deltanLinear,deltan_Mp] = tau_averages(deltan_now);
    [tauLinear,tau_Mp] = tau_averages(tau_now);

    figure(lifetime); 
    h(count) = loglog(deltanLinear,tauLinear,'o','MarkerSize',8);
    hold all;
    count = count+1;
    h(count) = loglog(deltan_Mp,tau_Mp,'s','MarkerSize',8);
    hold all;
    count = count+1;
    hold all; 
    loglog(deltan_now,tau_now,'k.');
    hold all;
end

xlabel('Excess carrier density (cm^-^3)','FontSize',30);
ylabel('Lifetime (\mus)','FontSize',30);
legend(h,'QSSPC','Arithmetic mean','Harmonic mean');
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
hgsave(lifetime,[saveStart 'QSSPC_PCPL']);
print(lifetime,'-dpng','-r0',[saveStart 'QSSPC_PCPL.png']);

