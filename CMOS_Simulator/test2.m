Ns =-2.5E21;    			% Substrate doping concentration (cm^-3)
Tox=10.8E-9;         		% Oxide thickness (cm)
q=1.6E-19;
Eg = 1.12*q;         % Bandgap of Si in eV
Nv = 1.02E19;      % effective density of states of valence band
Nc = 2.8E19;       % effective density of states of conduction band 
ni = 1E10;      % density of intrinsic carrier 
eps0 = 8.85E-14;   % Permittivity of free space (F/cm)
eps1 = 11.7;       % Relative permittivity of Si
eps1ox = 3.9;      % Relative permittivity of Si


Na=10e16;
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
Vend=0.4;
Vstrt=0; 
n=1;  %No. of Iterations
V=linspace(Vstrt,Vend,nx);     %Voltage vector
V(1)=0.2;
V(1000)=0.4;
V=V';


%defining A
A=zeros(1000,1000);
A(1,1)=1;
%A(1,2)=1;
for i=2:999
    A(i,i)=-2;
    A(i,i-1)=1;
    A(i,i+1)=1;
end
%A(1000,1000)=-2/dx2;
A(1000,1000)=1;
%A=defined;

%defining F

F=zeros(1000,1);
F(1,1)=0.2;
F(1000,1)=0.4;
for i=2:999
    F(i,1)=dx2*(-q/eps0)*(Na*(exp(-V(i,1)/vt)-1)-(ni^2/Na)*((exp(V(i,1)/vt))));
end
Fdash=zeros(1000,1000);
Fdash(1,1)=0;
F(1000,1)=0;
for i=2:999
    Fdash(i,i)=dx2*(q/(eps0*vt))*(Na*(exp(-V(i,1)/vt))+(ni^2/Na)*((exp(V(i,1)/vt))));
end

for j=0:1
    
    for i=2:999
    F(i,1)=dx2*(-q/eps0)*(Na*(exp(-V(i,1)/vt)-1)-(ni^2/Na)*((exp(V(i,1)/vt))));
    end
    for i=2:999
         Fdash(i,i)=dx2*(q/(eps0*vt))*(Na*(exp(-V(i,1)/vt))+(ni^2/Na)*((exp(V(i,1)/vt))));
    end
    V=V-(Fdash-A)^-1*(F-A*V);

end

%{
for i=2:999
    for j=1:15
        F=(q/(eps0*vt))*(Na*(exp(-V(i)/vt)-1)-(ni^2/Na)*((exp(V(i)/vt))));
        Fdash=(q/(eps0*vt))*(Na*(exp(-V(i)/vt))+ni^2/Na*(exp(V(i)/vt)));
        V(i)=V(i)+(V(i+1)+V(i-1)-2*V(i)-F*dx2)/(2+Fdash*dx2);
    end
end
%}
     plot(x(1:end-1),V(1:end));