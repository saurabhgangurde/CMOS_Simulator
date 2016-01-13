
tic
%****************initial definition!*********************
%
clear;             		% Clear memory and print header
cflag=1;           		% 0/1 quantum/classical
imax=5000;			% maximum number of newton iterations
%
%****************INPUT PARAMETERS!!!*********************
%
fid1=fopen('cvdata','w');	% Initialize output file; first arg is filename
T = 300;            		% temperature (K)
Ns =-2.5E19;    			% Substrate doping concentration (cm^-3)
				% (-/+ for n/p-type);
Tox=10.8E-7;         		% Oxide thickness (cm)
pd = 0;            		% consider poly depletion? 1: yes 0: no
Npoly=-5e20;        		% Poly doping (- means n-type; + means p-type)
Vfb=-0.88;				% 99 means autocalculate, else uses value
				% autocalculate uses Npoly and Ns (see below)
				% For metal gate, please input Vfb here
Dvc=  1e-4;       		% small signal voltage for capacitance
np=sign(Ns);			% (do not edit this line) n-sub or p-sub
Vstart=(0.2)*np*1;       	% Silicon Voltage: start (accumulation)
Vend=np*(-0.2);         		% Silicon Voltage: end   (inversion)
Nvstep=2;         		% number of voltage steps 
%
cimp0=0;			% Following parameters for implant
%cimp0=3.0e17;			% Gaussian implant: peak conc. (+: p-type)
%ximp0=0.02e-4;			% Gaussian implant: Rp
%simp0=0.05e-4;			% Gaussian implant: dRp (straggle)
%
%*********DEFINE THE RANGE FOR CALCULATION******************
%
if Nvstep>1;
  Vs0=linspace(Vstart,Vend,Nvstep);
else
  Vs0(1)=Vstart;
end
%
%*********calculate Vfb*****************
%
Ns=abs(Ns);
if Vfb==99
  if np*Npoly<0;
    Vfb=-np*0.026*log(abs(Npoly*Ns)/2.1045e20); 	% n+/nmos or p+/pmos
  else
    Vfb=np*0.026*log(abs(Npoly/Ns)); 	%  n+/nsub or p+/psub
  end
end
%
%****************CONSTANTS**********************
%
pi=3.14;
k = 8.61735E-5;    % Boltzmann constant in eV/K
beta=1/k/T;        % kT/q (1/eV) 
vt=0.0259;
hb= 6.58215e-16;   % Plank's constant in eV-s
eps0 = 8.86E-14;   % Permittivity of free space (F/cm)
eps1 = 11.7;       % Relative permittivity of Si
eps1ox = 3.9;      % Relative permittivity of Si
q0 = 1.602E-19;    % electron charge (C)
au = 0.5262E-8;    % atomic unit in cm
Eg = 1.12;         % Bandgap of Si in eV
Nv = 1.02E19;      % effective density of states of valence band
Nc = 2.8E19;       % effective density of states of conduction band 
ni = 1.45E10;      % density of intrinsic carrier 
pai= 3.14159;      
m0 = 9.109534E-31; % electron mass in Kg
m1 = 0.916;        % Electron effective mass in z direction for lower(3,6)
m2 = 0.19;         % Electron effective mass in z direction for higher(1,4,2,5)
md1 = 0.19;        % Density of state mass low energy valley (3,6)
md2 = 0.417;       % Density of state mass high energy valley (1,4,2,5)
mc1 = 0.19;        % Conductivity     mass low energy valley (3,6)
mc2 = 0.315;       % Conductivity     mass high energy valley (1,4,2,5)
me=m0;
mp=1.67E-27;
%
%****************NUMERICAL CONSTANTS******************
%
aregion=2000e-8; % size of classically allowed region (cm)
fregion=100e-8;	% size of classically forbidden region (cm)
N1 =50;		%number of mesh points in the fregion for poisson solution
N2 = 50;	%number of mesh points in the aregion for poisson solution
NS =50;		%number of mesh points in the aregion for Shrodinger eqn.
N=N1+N2;	% total number of mesh points for the poisson solution
dx0= fregion/real(N1);   		% mesh size in fregion
dx1= (aregion-fregion)/real(N2-1);	% mesh size in aregion
%
%**************** MESH GENERATION    **************** 
%
dx=zeros(N,1);			% initialize spacing array
xscale=zeros(N,1);		% initialize space array
dx(1:N1)=ones(N1,1).*dx0;    	% fregion spacing
xscale(2:N1+1)=dx0:dx0:fregion;	% fregion space array
dx(N1+1:N)=ones(N2,1).*dx1;  	% aregion spacing
xscale(N1+1:N)=fregion:dx1:aregion;		% aregion space array
%
%**************** IMPURITY PROFILE SETUP ************************        
%
Nad= Ns*ones(N,1);  		% constant substrate doping
%
if cimp0~=0;			% Gaussian implant
  Nad=Nad+np*cimp0.*exp(-(xscale-ximp0).^2./(simp0^2));
end
%
%****************CALCULATED PARAMETERS**************** 
%
if sign(np)==1
  Nep0=0.5*(Ns+sqrt(Ns*Ns+4*ni*ni));    % hole density at equilibrium
  Nen0=ni*ni/Nep0;                      % electron density at equilibrium
else
  Nen0=0.5*(Ns+sqrt(Ns*Ns+4*ni*ni)); % electron density at equilibrium
  Nep0=ni*ni/Nen0;                      % hole density at equilibrium
end
Ef=+k*T*log(Nep0/ni);                 % definition of Fermi
%
fprintf('Matrix for Poisson is %g by %g square \n',N,N);
Nen= zeros(N,1);  
Nep= zeros(N,1);  
Ne = zeros(N,1);  
Rho = zeros(N,1); 
V0 = zeros(N,1);
correctionN=zeros(N,1);
correctionP=zeros(N,1);
VS = zeros(N,1); 
A = zeros(N,N);  % Matrix for 2nd differential operator 
A(1,1)=1/dx0^2;	 %bondary condition Vsurface=Vs

gamma=1;
be=gamma*hb^2*q0*q0/(4*pi^2*me*6);
bp=0;
%bp=gamma*hb^2*q0/(4*pi^2*mp*6);
count=0;
%----------------------------------------------------
for j=2:N-1
    avgdx=(dx(j-1)+dx(j))/2;
    A(j,j-1) = 1/dx(j-1)/avgdx;
    A(j,j) = -(1/avgdx)*(1/dx(j-1)+1/dx(j));
    A(j,j+1) = 1/dx(j)/avgdx;
end;
A(N,N)=1/dx(N-1)^2;     
%
%************** BEGIN CALCULATIONS *************************
%
fprintf(fid1,'Vs          Vg          Cac         Lndc        Lnac        Lpdc        Lpac        Es          Qtot        Qinv        \n');
%
%**************LOOP FOR EACH VOLTAGE STEP*************************
%
for ivl=1:Nvstep
  Vss(1)=Vs0(ivl)-Dvc; % for capacitance calculation purpose (dV)
  Vss(2)=Vs0(ivl);
%
%**************LOOP FOR CAPACITANCE CALCULATION*******************
%
  for icc=1:2
    Vs=Vss(icc);
%
%**************POISSON EQUATION SETUP*****************************
%
    Rho(1)=Vs/dx0^2;  %bondary condition on the surface	 
    Rho(N)=0;  %bondary condition V(N)=0
    V0=A\Rho; %inital guess
    deltaNe = zeros(N,1);  
    deltaRho = zeros(N,1); 
    R = zeros(N,1); 
%
%***************** NEWTON LOOP ********************************** 
%
    for i=1:imax;  
        
      Nep=+Nep0*exp(-beta*V0); 

      % Electron
        Nen=+ni^2./Nep;
      
      %density gradient part
      %{
          for k=2:N-1
          correctionN(k)=-1*be*(sqrt(Nen(k+1))+sqrt(Nen(k-1))-2*sqrt(Nen(k)))/(dx(k)^2*sqrt(Nen(j)));
          correctionP(k)=-1*bp*(sqrt(Nep(k+1))+sqrt(Nep(k-1))-2*sqrt(Nep(k)))/(dx(k)^2*sqrt(Nep(j)));
      end;
      %}     
      Ne=-sign(np)*Nad-Nen+Nep; 	% Net charge density
      deltaNe=-beta*Nep-beta*Nen; 	% gradient for NR method

      Rho=-q0*Ne/eps0/eps1;
      deltaRho=-q0*deltaNe/eps0/eps1;
      
      sqrtNen=sqrt(Nen);
      sqrtNep=sqrt(Nep);
      correctionN=-1*be*((A*sqrtNen)./sqrt(Nen));
      correctionP=-1*bp*((A*sqrtNep)./sqrt(Nep));
      
      B=zeros(N,1);
      
      for row=1:N
          b=0;
          for col=1:N
              b=b+A(row,col)*sqrt(Nen0)*exp(V0(col)/(2*vt));
          end
          B(row)=b;%/sqrt(Nen0)*exp^(V0(row)/(2*vt));
      end
              
        
       Bdash=zeros(N,N);
       for row=1:N
          for col=1:N
              if row==col
                  Bdash(row,row)=-be*(0.5/(sqrt(Nen0)*vt))*(B(row)-A(row,row)*sqrt(Nen0)*exp(V0(row)/(2*vt)))*exp(-0.5*V0(row)/vt);
              else
                  Bdash(row,col)=-be*(0.5/(sqrt(Nen0)*vt))*A(row,col)*exp(0.5*(V0(col)-V0(row))/vt);
              end
          end
         
      end
      
    
      %Set up the Newton Raphson matrix
      NR=A;
      for j=2:N-1
        NR(j,j)=NR(j,j)-deltaRho(j);
      end
      NR=NR+A*Bdash;
      R=-A*(V0+correctionN+correctionP)+Rho;	
%
      R(1)=0;
      R(N)=0;
      deltaV0=NR\R;  % Potential in eV
      V0=V0+deltaV0; % Update the V0
%
      % convergence test
      dsort=max(abs(deltaV0));
      if(dsort<1.0e-12)
        break
      end 
      
    end 				% Newton loop
  end
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++
%plot(xscale,V0);
toc