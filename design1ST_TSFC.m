clear; 
clc; 
 
%all givens
nd = 0.95;
yd = 1.4;
nc = 0.9;
yc = 1.37;
nb = 0.97;
yb = 1.35;
nt = 0.92;
yt = 1.33;
nn = 0.98;
yn = 1.36;
nf = 0.92;
yf = 1.4;
nfn = 0.99;
 
syms pic;
syms B;
 
QR = 45*10^6;
R = 287;
Ta = 216.65;
Tmax = 1400:50:1800;
mach = 1.7;
pic = 16:1:40;
Pa = 7231.355;
 
%conditions of ambient (state a)
Toa = Ta*(1+(yd-1)/2*mach^2);
Poa = Pa*(1+(yd-1)/2*mach^2)^(yd/(yd-1));
u = mach*sqrt(R*yd*Ta);
 
%conditions after diffuser, before fan (state 2)
To2s = Toa;
To2 = (To2s-Ta)/nd+Ta;
Po2 = Pa*(To2s/Ta)^(yd/(yd-1));
 
%conditions after fan, before fan nozzle (state 8) 
Po8 = Po2*prf;
To8 = To2*(1+1/nf*(prf^((yf-1)/yf)-1));
 
%conditions after fan nozzle, before compressor
uef = sqrt(2*nfn*yf/(yf-1)*R*To8*(1-(Pa/Po8)^((yf-1)/yf)));
 
%conditions after compressor, before burner (state 3)
Po3 = pic*Po8;
To3s = To2*(Po3/Po8)^((yc-1)/yc);
To3 = (To3s-To8)/nc+To8;
cpc = yc/(yc-1)*R;
 
%conditions after burner, before turbine (state 4)
To4 = Tmax;
cp4 = yb/(yb-1)*R;
Po4 = Po3;
f = (To4/To3-1)/(QR/(cp4*To3)-To4/To3);
 
%conditions after turbine, before nozzle (state 5)
cpt = yt/(yt-1)*R;
To5 = (cpc*(To2 - To3) - B*cpc*(To8 - Toa) + To4*cpt*(f + 1))/(cpt*(f + 1)); 
    %cpt*(To4-To5)*(1+f) == cpc*(To3-To2) + B*cpc*(To8-Toa) solved for To5
To5s = To4-(To4-To5)/nt;
Po5 = Po4*(To5s/To4)^(yt/(yt-1));
 
%conditions after nozzle (state 6)
To6 = To5;
Po6 = Po5;
 
%conditions of exit flow (state 7)
P7 = Pa;
Po7 = Po6;
To7s = To6;
T7s = To7s*(P7/Po7)^((yn-1)/yn);
cp6 = yn/(yn-1)*R;
ue = sqrt(2*nn*cp6*(To7s-T7s));
%ue = sqrt(2*nn*yn/(yn-1)*R*To6*(1-(Pa/Po6)^((yn-1)/yn)));
 
 
 
%question(A): ST = thrust/ma
ST = (1+f)*ue+B*uef-(1+B)*u;
 
%question(B): TSFC = mf/thrust = (1/ST)*(mf/ma) = f/ST;
TSFC = f/ST;
 
%question(C): nth = 
% (kinetic energy of exhaust - kinetic energy of intake)/chemical energy of burned fuel
nth = ((1+f)*(ue^2)/2+B*(uef^2)/2-(1+B)*(u^2)/2)/(f*QR);
 
%question(D): np = 
% thrust power/(rate of change of kinetic energy of exhaust - intake)
np = ST*u/((1+f)*(ue^2)/2+B*(uef^2)/2-(1+B)*(u^2)/2);
 
%question(E): no = np*nth
no = nth*np;
 
 
 

