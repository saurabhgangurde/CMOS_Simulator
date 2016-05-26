function N= DG(V,N0,H,Ni)
bn=1e-8;
dim=length(V);
ktbyq=0.025;

display('in function');

for j=1:1000
    %j
    X=V-ktbyq*log(N0/Ni);
    G0= bn*H*sqrt(N0)+ sqrt(N0).*X;
 

    first_term=bn*H*0.5*diag(1./sqrt(N0));
    sqrt(N0);

    temp_vector=(0.5./sqrt(N0)).*V;
    second_term=diag(temp_vector);
    for i=1:dim-1
        second_term(i,i)=second_term(i,i)+sqrt(N0(i))*(V(i+1)-V(i))/(N0(i+1)-N0(i));
    end
    second_term(dim,dim)=second_term(dim,dim)+sqrt(N0(dim))*V(dim)/N0(dim);

    temp_vector=(ktbyq*0.5./sqrt(N0)).*log(N0/Ni);
    third_term=diag(temp_vector);
    for i=1:dim
        third_term(i,i)=third_term(i,i)+Ni/sqrt(N0(i));
    end

    %third_term
    deltaG= first_term+second_term+third_term;


    deltaN=-(deltaG)^-1*G0;
    N0=N0+deltaN;

   %dsort=max(abs(deltaN));
    %      if(dsort<norm(N0)/(1000*sqrt(length(N0))))
     %       break
      %    end
end
N=N0;
end