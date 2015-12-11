
h=(6.626*10^-34)/(2*3.14);
m=9.1*10^-31;
number=1000;
e=1.6*10^-19;
E=e*5;
length=10*10^-9;
dx=length/number;
wavefunc=zeros(number-1,2);
x=linspace(0,length,number);
x=x';

%potential=zeros(number,1);
potential=linspace(0,10,1000);
potential=potential';
potential_mean=zeros(number-1,1);
for i=1: number-1
    potential_mean(i,1)=(potential(i,1)+potential(i+1,1))/2;
end    
V=potential_mean*e;

k=zeros(number-1,1);

for i=1:number-1
    k(i,1)=sqrt(2*m*(E-V(i,1)))/h;
end

wavefunc1(1,1)=1;
wavefunc2=wavefunc1;
T1=zeros(2,2);
T2=zeros(2,2);

j=1;
T1(1,1)=exp(1i*k(j,1)*x(j+1,1));
T1(1,2)=exp(-1i*k(j,1)*x(j+1,1));
T1(2,1)=k(j,1)*exp(1i*k(j,1)*x(j+1,1));
T1(2,2)=-k(j,1)*exp(-1i*k(j,1)*x(j+1,1));
    
T2(1,1)=exp(1i*k(j+1,1)*x(j+1,1));
T2(1,2)=exp(-1i*k(j+1,1)*x(j+1,1));
T2(2,1)=k(j+1,1)*exp(1i*k(j+1,1)*x(j+1,1));
T2(2,2)=-k(j+1,1)*exp(-1i*k(j+1,1)*x(j+1,1));
T=T2\T1;
%{
syms a2;syms b2;syms b1;
S=solve([a2 == T(1,1) + T(1,2) * b1,b2 == T(2,1) + T(2,2) * b1,a2^2 + b1^2 == 1 + b2^2],[b1,a2,b2]);
wavefunc1(1,2)=S.b1(1);
wavefunc1(2,1)=S.a2(1);
wavefunc1(2,2)=S.b2(1);

wavefunc2(1,2)=S.b1(2);
wavefunc2(2,1)=S.a2(2);
wavefunc2(2,2)=S.b2(2);

%}
p=T(1,2)^2+1-T(2,2)^2;
q=(2*T(1,1)*T(1,2))-(2*T(2,1)*T(2,2));
r=T(1,1)^2-1-T(2,1)^2;

if(p==0)
    wavefunc1(1,2)=-q/r;
    wavefunc2(1,2)=-q/r;
    
else
    wavefunc1(1,2)=(-q+sqrt(q^2-(4*p*r)))/(2*p);
    wavefunc2(1,2)=(-q-sqrt(q^2-(4*p*r)))/(2*p);
end
wavefunc1(2,1)=T(1,1) +T(1,2)*wavefunc1(1,2);
wavefunc1(2,2)=T(2,1) +T(2,2)*wavefunc1(1,2);

wavefunc2(2,1)=T(1,1) +T(1,2)*wavefunc2(1,2);
wavefunc2(2,2)=T(2,1) +T(2,2)*wavefunc2(1,2);


for j=2:number-2
    
    T1(1,1)=exp(1i*k(j,1)*x(j+1,1));
    T1(1,2)=exp(-1i*k(j,1)*x(j+1,1));
    T1(2,1)=k(j,1)*exp(1i*k(j,1)*x(j+1,1));
    T1(2,2)=-k(j,1)*exp(-1i*k(j,1)*x(j+1,1));
    
    a11=exp(1i*k(j+1,1)*x(j+1,1));
    a12=exp(-1i*k(j+1,1)*x(j+1,1));
    a21=k(j+1,1)*exp(1i*k(j+1,1)*x(j+1,1));
    a22=-k(j+1,1)*exp(-1i*k(j+1,1)*x(j+1,1));
     
    T2(1,1)=a22;
    T2(1,2)=-a12;
    T2(2,1)=-a21;
    T2(1,2)=a11;
    
    T=(T2*T1)/(a11*a22-a21*a12);
    wavefunc1(j+1,:)=wavefunc1(j,:)*T;
    wavefunc2(j+1,:)=wavefunc2(j,:)*T;
    
end
%{
func1=zeros(number,1);
func2=zeros(number,1);

func1(1,1)=wavefunc1(1,1)*exp(1i*k(1,1)*x(1,1))+wavefunc1(1,2)*exp(-1i*k(1,1)*x(1,1));
func2(1,1)=wavefunc2(1,1)*exp(1i*k(1,1)*x(1,1))+wavefunc2(1,2)*exp(-1i*k(1,1)*x(1,1));

for j=1:number-1;
 
    func1(j+1,1)=wavefunc1(j,1)*exp(1i*k(j,1)*x(j+1,1))+wavefunc1(j,2)*exp(-1i*k(j,1)*x(j+1,1));
    func2(j+1,1)=wavefunc2(j,1)*exp(1i*k(j,1)*x(j+1,1))+wavefunc2(j,2)*exp(-1i*k(j,1)*x(j+1,1));

    
end    


A=zeros(number,number);
for i=1:number-1
    A(i,i)=-1/dx;
    A(i,i+1)=1/dx;
end
A(number,number-1)=-1/dx;
A(number,number)=1/dx;

J1=(1/(2*m))*(real(transpose(conj(func1))*(-1i*h)*A*func1));
J2=(1/(2*m))*(real(transpose(conj(func2))*(-1i*h)*A*func2));
%}
%{
hold on;
plot(x,J1,'r');
plot(x,J2,'b');

%}
display('done');

