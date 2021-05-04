% thermalBilinearConvec.m
% MN-ExFinal1Q-2018-19P4A 
%Solving 2D thermal problems using quadrilateral 
% elements with convection boundary conditions.   

%%
clearvars
close all

f=0.0;

%Perm.A
perm='I';
kc=0.7;
TLeft=100.0;
TRight=200.0;
TinfTop=20.0;
betaTop=2.0;
TinfBot=10.0;
betaBot=3.0;
p = [2.3,2.2];

eval('meshAleta2DQuad')

numNod=size(nodes,1);
numElem=size(elem,1);


indTop=find(5*nodes(:,2) - 2*nodes(:,1) -10 > -0.01);
indBot=find(nodes(:,2)<0.01);
indLeft=find(nodes(:,1)<0.01);
indRight=find(nodes(:,1)>4.99); 

figure();
numbering = 0;
plotElementsOld(nodes, elem, numbering);
hold on

plot(nodes(indTop,1),nodes(indTop,2),'or','markerFaceColor','r');
plot(nodes(indBot,1),nodes(indBot,2),'og','markerFaceColor','g');
plot(nodes(indRight,1),nodes(indRight,2),'om','markerFaceColor','m');
plot(nodes(indLeft,1),nodes(indLeft,2),'ob','markerFaceColor','b');

hold off

%%

a11=kc;
a12=0.0;
a21=a12;
a22=a11;
a00=0;

coeff=[a11,a12,a21,a22,a00,f];

K=zeros(numNod);
Q=zeros(numNod,1);
F=zeros(numNod,1);
for e=1:numElem
    [Ke,Fe]=bilinearQuadElement(coeff,nodes,elem,e);
    rows=elem(e,:);
    cols=rows;
    K(rows,cols)=K(rows,cols)+Ke;
    if coeff(6)~=0
        F(rows)=F(rows)+Fe;
    end
end
Fini=F;
Kini=K;

fixedNodes=[indLeft',indRight'];
freeNodes=setdiff(1:numNod,fixedNodes);

%BC
u = zeros(numNod,1);
%----------- Natural B.C.: convection
indCV=indTop';
[K,Q]=applyConvQuad(indCV,betaTop,TinfTop,K,Q,nodes,elem);
indCV=indBot';
[K,Q]=applyConvQuad(indCV,betaBot,TinfBot,K,Q,nodes,elem);

%----------- Essential BC
u(indLeft,1)=TLeft;
u(indRight,1)=TRight;

% Modify F
%F = F - K(:,fixedNodes)*u(fixedNodes);

%Reduced system
Im = F(freeNodes,1) + Q(freeNodes,1) - K(freeNodes,fixedNodes)*u(fixedNodes,1);
Km=K(freeNodes,freeNodes);

%Solve the Reduced System
um = Km\Im;
u(freeNodes,1)=um;

%Post Process: compute Q's:
Q = Kini*u-Fini;

%Numerical output
% R=[(1:numNod)',nodes,u,Q];
% fprintf(1,'\n%5s%8s%14s%14s%14s\n','Nod.','X','Y','T','Q')
% fprintf(1,'%4d%14.5e%14.5e%14.5e%14.5e\n',R')

%Graphical output
titol='Equation solution';
colorScale='jet';
plotContourSolution(nodes,elem,u,titol,colorScale);

bar(1:52) = '-';
fprintf('%s',bar)
%fprintf('\nPERM.%s\n',perm)
%fprintf('%s',bar)
[minTemp,nodMinTemp]=min(u);
fprintf('\nPART (A)\n')  %% Fancy output. Dont try this at exams!!!
fprintf('Compute the minimum temperature achieved at the\n') 
fprintf('nodes and the coordinates (x,y) of the node with\n')
fprintf('minimum temperature:\n\n')
fprintf('Min Temperature: %.5e\n',minTemp)
fprintf('Node Coordinates: (%.5e, %.5e)\n',nodes(nodMinTemp,:))
fprintf('Hint: the temperature of node 156 is: %.5e\n',u(156))

%Compute the temperature for point p=[2.3, 2.2].
fprintf(1,'%s',bar);
fprintf('\nPART (B)\n')
fprintf('Using the point p = [%.1f,%.1f], find thenumber of\n',p)
fprintf('the element to which p belongs, its second bary-\n')
fprintf('centric coordinate and the interpolated temperature\n')
fprintf('Tp of the point p.\n\n')

for e = 1:numElem
    n1 = elem(e,1);
    n2 = elem(e,2);
    n3 = elem(e,3);
    n4 = elem(e,4);
    v1 = nodes(n1,:);
    v2 = nodes(n2,:);
    v3 = nodes(n3,:);
    v4 = nodes(n4,:);
    vertexs = [v1;v2;v3;v4];
    [alphas,isInside] = baryCoordQuad(vertexs,p);
    if (isInside >= 1)
        break;
    end
end

tempInterp = alphas(1)*u(n1) + alphas(2)*u(n2) + alphas(3)*u(n3) + ...
    alphas(4)*u(n4);
fprintf('Element number = %d\n',e);
fprintf('        alpha2 = %.5e\n',alphas(2));
fprintf('         tempP = %.5e\n',tempInterp);
% fprintf(1,...
%     '\n\nExercise 1. Temperature at point P = (%.5g,%5g):\n',...
%     p(1,1),p(1,2))
% fprintf('Element: %d\nNods.: [%d,%d,%d,%d]\nInterp.Temp: %.5e\n\n',...
%     e,n1,n2,n3,n4,tempInterp);

fprintf(1,'%s',bar)
fprintf('\nPART (C)\n')
fprintf('The number of nodes whose temperature differs from\n')
fprintf('Tp by less than 5 degrees.\n\n')
nods1 = find(abs(tempInterp-u) < 5.0);
fprintf('Number of Nodes: %d\n',size(nods1,1))
nods2 = find(abs(tempInterp-u) < 2.0);
fprintf('Hint: The number of points for a difference of 2\n')
fprintf('degrees is: %d\n',size(nods2,1))
fprintf(1,'%s\n',bar);