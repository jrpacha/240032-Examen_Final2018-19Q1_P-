function [w,pt] = gaussValues1D(n)
%FUNCTION [W,PT]= GAUSSVALUES1D(N)
% points and weights for the Gauss-Legendre quadrature 
% formulas.

switch (n)
    case 1
        w=2; pt=0;
    case 2
        w=[1,1]; pt=[-1/sqrt(3), 1/sqrt(3)];
    case 3
        w=[5/9, 8/9, 5/9]; pt=[-sqrt(3/5), 0, sqrt(3/5)];
    case 4
        a=(18-sqrt(30))/36;
        b=(18+sqrt(30))/36;
        w=[a,b,b,a]; 
        a=sqrt(3/7+2*sqrt(6/5)/7);
        b=sqrt(3/7-2*sqrt(6/5)/7);
        pt=[-a,-b,b,a];
    case 5
        a=(322-13*sqrt(70))/900;
        b=(322+13*sqrt(70))/900;
        w=[a,b,128/225,b,a];
        a=sqrt(5+2*sqrt(10/7))/3;
        b=sqrt(5-2*sqrt(10/7))/3;
        pt=[-a,-b,0,b,a];
    otherwise
        error('No data are defined for this value');
end

end

