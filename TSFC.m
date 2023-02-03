function thrustspecificfuelconsumption = TSFC(To_max,pic,B,prf)
 
%all givens
nd = 0.95;
yd = 1.4;
nc = 0.9;
yc = 1.37;
nb = 0.97;
yb = 1.35;
pib = 0.95;
nt = 0.92;
yt = 1.33;
nn = 0.98;
yn = 1.36;

nf = 0.92;
yf = 1.4;
nfn = 0.99;
yfn = 1.4;

 
QR = 45*10^6;
R = 287;
Ta = 216.65;
mach = 1.7;
Pa = 7231.355;
 
%conditions of ambient (state a)
Toa = Ta*(1+(yd-1)/2*mach^2);
Poa = Pa*(1+(yd-1)/2*mach^2)^(yd/(yd-1));
u = mach*sqrt(R*yd*Ta);
 
%conditions after diffuser, before fan (state 2)
To2 = Toa;
To2s = (To2-Ta)*nd+Ta;
Po2 = Pa*(To2s/Ta)^(yd/(yd-1));
 
%conditions after fan, before fan nozzle (state 8) 
Po8 = Po2*prf;
To8s = To2*prf^((yf-1)/yf);
To8 = (To8s-To2)/nf+To2;
cpf = yf/(yf-1)*R;
 
%conditions after fan nozzle, before compressor
cpfn = yfn/(yfn-1)*R;
uef = sqrt(2*nfn*cpfn*To8*(1-(Pa/Po8)^((yfn-1)/yfn)));
 
%conditions after compressor, before burner (state 3)
Po3 = pic*Po8;
To3 = To8*(pic^((yc-1)/(yc*nc)));
cpc = yc/(yc-1)*R;
 
%conditions after burner, before turbine (state 4)
To4 = To_max;
cpb = yb/(yb-1)*R;
Po4 = pib*Po3;
f = (To4-To3)/(nb*QR/cpb-To4);
 
%conditions after turbine, before nozzle (state 5)
cpt = yt/(yt-1)*R;
To5 = To4 - (cpc*(To3-To8)+(1+B)*cpf*(To8-To2))/((1+f)*cpt);
    %(1+f)*cpt*(To4-To5) = cpc*(To3-To8) + (1+B)*cpf*(To8-To2) solved for To5
Po5 = Po4*(To5/To4)^(yt*nt/(yt-1));
 
%conditions after nozzle (state 6)
P6 = Pa;
Po6 = Po5;
To6s = To5;
T6s = To6s*(P6/Po6)^((yn-1)/yn);
cpn = yn/(yn-1)*R;
ue = sqrt(2*nn*cpn*(To6s-T6s));
 
 

st_bare = (1+f)*ue+B*uef-(1+B)*u;
sthrust = st_bare/(1.04+0.01*B^1.2);

thrustspecificfuelconsumption = f/sthrust;
%specificthrust = sthrust/(1+B);
%fuelairratio = f;
%vcold = uef;
%vhot = ue;
end
 
 
 

