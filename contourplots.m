sThrust = zeros(9,13,21,6);
sFC = 999999.*ones(9,13,21,6);

for t_max = 1400:50:1800 %9 loops
    for pi_c = 16:2:40 %13 loops
        for beta = 0:0.5:10 %21 loops
            for pi_f = 1:0.2:2 %6 loops
                if isreal(ST(t_max,pi_c,beta,pi_f)) && ST(t_max,pi_c,beta,pi_f) > 0
                    sThrust((t_max-1400)/50+1,(pi_c-16)/2+1,beta/0.5+1,round((pi_f-1)/0.2+1)) = ST(t_max,pi_c,beta,pi_f).*(1+beta);
                end
               
                if isreal(TSFC(t_max,pi_c,beta,pi_f)) && TSFC(t_max,pi_c,beta,pi_f) > 0
                    sFC((t_max-1400)/50+1,(pi_c-16)/2+1,beta/0.5+1,round((pi_f-1)/0.2+1)) = TSFC(t_max,pi_c,beta,pi_f);                
                end
            end
        end
        
    end
end

invTSFC = 1./sFC;
normInvTSFC = invTSFC./(max(invTSFC,[],'all')); %normalized reciprocal of TSFC data, 1.00 is highest value

%fpr contour from 1 to 2
%bypass ratio contour from 0 to 5
%to4 = 1700, pic = 22
testrangeST = zeros(701,101); %rows are constant bypass ratio, columns constant FPR
testrangeTSFC = zeros(701,101);
for var1 = 1:701
    for var2 = 1:101
        testrangeST(var1,var2) = ST(1700,22,(var1-1)/100,1+(var2-1)/100);
        testrangeTSFC(var1,var2) = 1000.*TSFC(1700,22,(var1-1)/100,1+(var2-1)/100);
    end
end

[C,h] = contour(1:0.01:2,0:0.01:7,testrangeST,'red');
clabel(C,h);
xlabel("Fan Pressure Ratio");
ylabel("Bypass Ratio β");
title("Specific Thrust (^{N}/_{kg/s}) at T_0_4 = 1700K, π_c = 22");

figure();
[C,h] = contour(1:0.01:2,0:0.01:7,testrangeTSFC,'red');
clabel(C,h);
xlabel("Fan Pressure Ratio");
ylabel("Bypass Ratio β");
title("Thrust Specific Fuel Consumption (^{kg/s}/_{kN}) at T_0_4 = 1700K, π_c = 22");

Ta = 216.65; %kelvin
Pa = 7231.355; %pascals
M = 1.7;
uIn = M*sqrt(Ta*287*1.4); %inlet velocity
B = 1.5; %bypass ratio
Qr = 45*10^6; %J/kg

%To4 from 1400 to 1800
%compression ratio from 16 to 40
%fpr = 2, bypass ratio = 1.5 carpet plot
testrangeST = zeros(9,13); %rows are constant To4, columns constant pi_c
testrangeTSFC = zeros(9,13);
for var1 = 0:8
    for var2 = 0:12
        testrangeST(var1+1,var2+1) = ST(1400+var1*50,16+var2*2,1.5,2);
        testrangeTSFC(var1+1,var2+1) = 1000.*TSFC(1400+var1*50,16+var2*2,1.5,2);
    end
end

figure();
plot(testrangeST(1,:),testrangeTSFC(1,:),'red');
text(testrangeST(1,1)-15,testrangeTSFC(1,1)+.0001, 'T_0_4 = 1400K');
hold on
for temploop = 2:9
    plot(testrangeST(temploop,:),testrangeTSFC(temploop,:),'red');
    text(testrangeST(temploop,1)-15,testrangeTSFC(temploop,1)+.0001, append('T_0_4 = ',num2str(1400+(temploop-1)*50),'K'));
end
for presloop = 1:13
    plot(testrangeST(:,presloop),testrangeTSFC(:,presloop),'blue');
    text(testrangeST(9,presloop),testrangeTSFC(9,presloop)-.0001, append('π_c = ',num2str(16+(presloop-1)*2)));
end
xlim([60 360]);
ylim([.012 .028]);
xlabel("Specific Thrust (^{N}/_{kg/s})");
ylabel("Thrust Specific Fuel Consumption (^{kg/s}/_{kN})");
title("Carpet Plot of TSFC and ST vs T_0_4 and π_c for π_f = 2 and β = 1.5");

density = Pa/(287*Ta);
diameter = 1.6;
area = pi*(diameter/2)^2;
STmin = 80000/((density*uIn*area)*(1+B));
xline(STmin,"--k",'Minimum Specific Thrust');

LD = 7.5;
fuelFraction = 0.4;
TSFCmax = (log(1/(1-fuelFraction))*LD*uIn/9.81)/8000;
yline(TSFCmax,"--k", {'Maximum','TSFC'});
hold off

%contours in pi_c T04 plane
testrangeST = zeros(401,2401); %rows are constant temp, columns constant compression ratio
testrangeTSFC = zeros(401,2401); %rows are constant temp, columns constant compression ratio
testrangeThermEff = zeros(401,2401); %rows are constant temp, columns constant compression ratio
testrangePropEff = zeros(401,2401); %rows are constant temp, columns constant compression ratio
testrangeOverallEff = zeros(401,2401); %rows are constant temp, columns constant compression ratio

for var1 = 1:401
    for var2 = 1:2401
        testrangeST(var1,var2) = ST(1400+var1-1,16+(var2-1)/100,1.5,2);
        testrangeTSFC(var1,var2) = 1000.*TSFC(1400+var1-1,16+(var2-1)/100,1.5,2);
            [uCold,uHot] = exitVel(1400+var1-1,16+(var2-1)/100,1.5,2);
            f = fuelair(1400+var1-1,16+(var2-1)/100,1.5,2);
            deltaKE = 0.5*((1+f)*uHot^2+B*uCold^2-(1+B)*uIn^2);
            energyConsump = f*Qr;
        testrangeThermEff(var1,var2) = deltaKE/energyConsump;
            thrust_ma = testrangeST(var1,var2);
            thrustPower = thrust_ma*uIn;
        testrangePropEff(var1,var2) = thrustPower/deltaKE;
        testrangeOverallEff(var1,var2) = thrustPower/energyConsump;
    end
end

figure();
[C,h] = contour(16:0.01:40,1400:1:1800,testrangeST,'red');
clabel(C,h);
xlabel("Compression Ratio π_c");
ylabel("Max Temperature T_0_4");
title("Specific Thrust (^{N}/_{kg/s}) for β = 1.5, FPR = 2");

figure();
[C,h] = contour(16:0.01:40,1400:1:1800,testrangeTSFC,'red');
clabel(C,h);
xlabel("Compression Ratio π_c");
ylabel("Max Temperature T_0_4");
title("Thrust Specific Fuel Consumption (^{kg/s}/_{kN}) for β = 1.5, FPR = 2");

figure();
[C,h] = contour(16:0.01:40,1400:1:1800,testrangeThermEff,'red');
clabel(C,h);
xlabel("Compression Ratio π_c");
ylabel("Max Temperature T_0_4");
title("Thermal Efficiency for β = 1.5, FPR = 2");

figure();
[C,h] = contour(16:0.01:40,1400:1:1800,testrangePropEff,'red');
clabel(C,h);
xlabel("Compression Ratio π_c");
ylabel("Max Temperature T_0_4");
title("Propulsion Efficiency for β = 1.5, FPR = 2");

figure();
[C,h] = contour(16:0.01:40,1400:1:1800,testrangeOverallEff,'red');
clabel(C,h);
xlabel("Compression Ratio π_c");
ylabel("Max Temperature T_0_4");
title("Overall Efficiency for β = 1.5, FPR = 2");







