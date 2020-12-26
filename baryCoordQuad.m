function [alphas, isInside] = baryCoordQuad(vertices,point)
% FUNCTION BARYCOORDQUAD
% Computes the barycentric coordinates of a point 
% w.r.t. a given quadrilateral of vertexs v1=[x1,y1], 
% v2=[x2,y2], v3=[x3,y3], v4 = [x4,y4]. 
%
% INPUT
%   vertices: vertices of the quad as a 4x2 matrix
%             vertices = [v1,v2,v3,v4]'.
%      point: row matrix with the coordinates (x,y) of
%             the point, i.e.: point=[x,y]. 
% OUTPUT
%     alphas: row matrix with the barycentric coordina-
%             tes of the point.
%   isInside: equals 1 if the point belongs to the 
%             quadrilater, 0 otherwise. 

isInside = 1;

% define the previous notation

v1 = vertices(1,:);
v2 = vertices(2,:);
v3 = vertices(3,:);
v4 = vertices(4,:);

a=(v1-point);
b=v2-v1;
c=v4-v1;
d=v1-v2-v4+v3;

if (norm(d) <1.d-14) %this not a general quadrilateral,
                     %it is the case of rectangular meshes!!!
                     %the main equation becomes a+b*lambda+c*mu=0
                     %that can be solved directly as a 2x2 system
    M=[b',c'];
    t=-a';
    xx=M\t;
    w1=xx(1);
    w2=w1;
    u1=xx(2);
    u2=u1;
else
    % compute 2on order equation A mu^2 + B mu + C=0
    % as the vertices are 2D, we add a zero third component
    % to compute cross products.
    A=cross([c,0],[d,0]); %must be 3D vectors
    B=cross([c,0],[b,0])+cross([a,0],[d,0]);
    C=cross([a,0],[b,0]);
    % Only third component is needed (the other two are zero)
    A=A(3);
    B=B(3);
    C=C(3);
    %
    % Check for unique solutions
    %
    if (abs(A)<1.e-14)
        u1= -C/B;
        u2=u1;
    else
     %
     % Check for non complex solutions
     %
        if (B^2-4*A*C >0)
            u1=(-B+sqrt(B^2-4*A*C))/(2*A);
            u2=(-B-sqrt(B^2-4*A*C))/(2*A);
        else %complex solution
            u1=-1000;
            u2=u1;
        end
    end
    %
    % compute 2on order equation A lambda^2 + B lambda + C=0
    A=cross([b,0],[d,0]); %must be 3D vectors
    B=cross([b,0],[c,0])+cross([a,0],[d,0]);
    C=cross([a,0],[c,0]);
    % Only third component is needed (the other two are zero)
    A=A(3);
    B=B(3);
    C=C(3);
    %
    % Check for unique solutions
    %
    if (abs(A)<1.e-14)
        w1= -C/B;
        w2=w1;
    else
        %
        % Check for non complex solutions
        %
        if (B^2-4*A*C >0)
            w1=(-B+sqrt(B^2-4*A*C))/(2*A);
            w2=(-B-sqrt(B^2-4*A*C))/(2*A);
        else %complex solution
            w1=-1000;
            w2=w1;
        end
    end
end
%
mu=-10000; %stupid value small enough
if (u1>=0 && u1<=1)
    mu=u1;
end
if (u2>=0 && u2<=1)
  mu=u2;
end
%
lambda=-10000; %stupid value small enough
if(w1>=0 && w1<=1)
    lambda=w1;
end
if(w2>=0 && w2<=1)
    lambda=w2;
end
%[mu,lambda] %parameters
% Barycentric coordinates
alpha1=(1-mu)*(1-lambda);
alpha2=lambda*(1-mu);
alpha3=mu*lambda;
alpha4=(1-lambda)*mu;
alphas=[alpha1,alpha2,alpha3,alpha4];

% Finally, check if the point is inside the quadrilateral
if min(alphas) < 0
    isInside = 0;
end

end

