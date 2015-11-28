Ns =-2.5E21;    			% Substrate doping concentration (cm^-3)
Tox=10.8E-9;         		% Oxide thickness (cm)
q=1.6E-19;
Eg = 1.12*q;         % Bandgap of Si in eV
Nv = 1.02E19;      % effective density of states of valence band
Nc = 2.8E19;       % effective density of states of conduction band 
ni = 1.45E16;      % density of intrinsic carrier 
eps0 = 8.85E-12;   % Permittivity of free space (F/cm)
eps1 = 11.7;       % Relative permittivity of Si
eps1ox = 3.9;      % Relative permittivity of Si


Na=10e22;
%Na=0;
vt=0.0259;
%nx=1000; %Number of steps in space(x)
%niter=1000;                      %Number of iterations 
%dx=1e-9/nx;                     %Width of space step(x)
%x=1e-9*(1:dx:nx);                        %Range of x(0,2) and specifying the grid points

L=1000e-6;
nx=1000;
dx=L/nx;
dx2=dx^2;
x=0:dx:L;
Vend=0.1;
Vstrt=1; 
n=1;  %No. of Iterations
%B=zeros(1,n);
V=linspace(Vstrt,Vend,nx);     %Voltage vector
V=V';

%defining A
A=zeros(998,998);
A(1,1)=-2/dx2;
A(1,2)=1/dx2;
for i=2:997
    A(i,i)=-2/dx2;
    A(i,i-1)=1/dx2;
    A(i,i+1)=1/dx2;
end
A(998,998)=-2/dx2;
A(998,997)=1/dx2;
%A=defined;

%defining F
F=zeros(998,1);
for i=1:998
    F(i,1)=(-q/eps0)*(Na*(exp(-V(i,1)/vt)-2)-(ni^2/Na)*((exp(V(i,1)/vt)-1)));
end
Fdash=zeros(998,998);
for i=2:998
    Fdash(i,i)=(q/(eps0*vt))*(Na*(exp(-V(i,1)/vt))+(ni^2/Na)*((exp(V(i,1)/vt))));
end


for j=0:0
    display('sdsd');
    V(2:999)=V(2:999)+((Fdash-A)^-1)*(F-A*V(2:999));
    for i=2:998
    F(i,1)=(-q/eps0)*(Na*(exp(-V(i,1)/vt)-2)-(ni^2/Na)*((exp(V(i,1)/vt)-1)));
    end
    for i=2:998
         Fdash=(q/(eps0*vt))*(Na*(exp(-V(i,1)/vt))+(ni^2/Na)*((exp(V(i,1)/vt))));
    end

end

     plot(x(1:end-1),V(1:end));