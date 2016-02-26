function [z,ncg]=cgp(A,b,del)
    r=-b;
    d=b;
    z=zeros(length(b),1);
    ncg=0;
    while norm(r)>1e-7
        ncg=ncg+1;
        if d'*A*d<=0
            dn=norm(d)^2;
            zn=norm(z)^2;
            dz=z'*d;
            tau=(sqrt(dz^2+dn*(del^2-zn))-dz)/dn;
            tau1=(-sqrt(dz^2+dn*(del^2-zn))-dz)/dn;
            z=z+tau*d;
            z1=z+tau1*d;
            if -b*(z-z1)+0.5*z'*A*z-0.5*z1'*A*z1>0
                z=z1;
            end
            break;
        end
        alpha=r'*r/(d'*A*d);
        z0=z;
        z=z+alpha*d;
        if norm(z)>=del
            dn=norm(d)^2;
            zn=norm(z0)^2;
            dz=z0'*d;
            tau=(sqrt(dz^2+dn*(del^2-zn))-dz)/dn;
            z=z0+tau*d;
            break;
        end
        r0=r;
        r=r+alpha*A*d;
        beta=r'*r/(r0'*r0);
        d=-r+beta*d;
    end
end