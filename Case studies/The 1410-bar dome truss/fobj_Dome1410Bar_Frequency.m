function [W,PFIT,sol] = fobj_Dome1410Bar_Frequency(X,se1,se2)

E = 20e10;
Density = 7850;
Mass = 100;
minfrq = [7 9];
FF = zeros(117,117);
MM = zeros(117,117);
W = 0;
AA = [X';X'];
NUM_S = 30;
NUM_NODE = 13;
THETA = 2*pi/NUM_S;
X_COORDS = [1 3 5 7 9 11 13 1.989 3.978 5.967 7.956 9.945 11.934];
Y_COORDS = [0 0 0 0 0 0 0 0.209 0.418 0.627 0.836 1.0453 1.2543];
Z_COORDS = [4 3.75 3.25 2.75 2 1.25 0 3 2.75 2.25 1.75 1 -0.5];
coorda = [X_COORDS;Y_COORDS];

nodes = [];
for ns = 1:3
    trans = [cos((ns-1)*THETA) -sin((ns-1)*THETA);
        sin((ns-1)*THETA) cos((ns-1)*THETA)];
    nodes = [nodes;[trans*coorda;Z_COORDS]'];
end

B = [];
for ns = 1:2
    firstNodeCur = (ns-1)*NUM_NODE + 1;
    firstNodeNext = firstNodeCur+NUM_NODE;
    B = [B;[firstNodeCur:firstNodeCur+5;firstNodeCur+1:firstNodeCur+6]'];
    B = [B;[firstNodeCur+7:firstNodeCur+11;firstNodeCur+8:firstNodeCur+12]'];
    B = [B;[firstNodeCur+1:firstNodeCur+6;firstNodeCur+7:firstNodeCur+12]'];
    B = [B;[firstNodeCur+1:firstNodeCur+5;firstNodeCur+8:firstNodeCur+12]'];
    B = [B;[firstNodeCur,firstNodeCur+7]];
    B = [B;[firstNodeCur:firstNodeCur+5;firstNodeNext:firstNodeNext+5]'];
    B = [B;[firstNodeCur+7:firstNodeCur+12;firstNodeNext+7:firstNodeNext+12]'];
    B = [B;[firstNodeCur+7:firstNodeCur+12;firstNodeNext:firstNodeNext+5]'];
    B = [B;[firstNodeCur+7:firstNodeCur+12;firstNodeNext+1:firstNodeNext+6]'];
end

Constr = zeros(3,39);
for j = 1:3
    Constr(:,13*(j-1)+7) = 1;
end

Ur = 1-Constr;
f = find(Ur);

for i = 1:94
    XTOOL = nodes(B(i,2),1)-nodes(B(i,1),1);
    YTOOL = nodes(B(i,2),2)-nodes(B(i,1),2);
    ZTOOL = nodes(B(i,2),3)-nodes(B(i,1),3);
    ELEMTOOL = sqrt(XTOOL^2+YTOOL^2+ZTOOL^2);
    roi = sqrt(nodes(B(i,1),1)^2+nodes(B(i,1),2)^2);
    roj = sqrt(nodes(B(i,2),1)^2+nodes(B(i,2),2)^2);
    loi = nodes(B(i,1),1)/roi;
    moi = nodes(B(i,1),2)/roi;
    loj = nodes(B(i,2),1)/roj;
    moj = nodes(B(i,2),2)/roj;
    Roi = [loi -moi 0; moi loi 0; 0 0 1];
    Roj = [loj -moj 0; moj loj 0; 0 0 1];
    R = [Roi zeros(3,3);zeros(3,3) Roj];
    XT = XTOOL/ELEMTOOL;
    YT = YTOOL/ELEMTOOL;
    ZT = ZTOOL/ELEMTOOL;
    se = [XT^2,XT*YT,XT*ZT;XT*YT,YT^2,YT*ZT;XT*ZT,YT*ZT,ZT^2];
    I2 = (B(i,1)-1)*3;
    H2 = (B(i,2)-1)*3;
    EAL = AA(i,1)*E/ELEMTOOL;
    APL6 = AA(i,1)*Density*ELEMTOOL/6;
    Weight = AA(i,1)*ELEMTOOL*Density;
    W = W+15*Weight;
    SE = EAL*[se -se;-se se];
    K_cylindrical_coordinate = R'*SE*R;
    elementdof = [I2+1 I2+2 I2+3 H2+1 H2+2 H2+3];
    FF(elementdof,elementdof) = FF(elementdof,elementdof)+K_cylindrical_coordinate;
    ME = APL6*[2 0 0 1 0 0;0 2 0 0 1 0;0 0 2 0 0 1;1 0 0 2 0 0;0 1 0 0 2 0;0 0 1 0 0 2];
    M_cylindrical_coordinate = R'*ME*R;
    MM(elementdof,elementdof) = MM(elementdof,elementdof)+M_cylindrical_coordinate;
end

for j = 1:39
    MM(3*j-2,3*j-2) = MM(3*j-2,3*j-2)+Mass;
    MM(3*j-1,3*j-1) = MM(3*j-1,3*j-1)+Mass;
    MM(3*j,3*j) = MM(3*j,3*j)+Mass;
end

MM = MM(f,f);
FF = FF(f,f);
AK = FF(37:72,37:72);
AM = MM(37:72,37:72);
BK = FF(1:36,37:72);
BM = MM(1:36,37:72);
sw = [];

for J = 1:14
    LANDA1 = complex(cos(2*J*pi/30),sin(2*J*pi/30));
    LANDA2 = complex(cos(2*J*pi/30),-sin(2*J*pi/30));
    kj = AK+LANDA1*BK+LANDA2*BK';
    mj = AM+LANDA1*BM+LANDA2*BM';
    sw = [sw;sqrt(real(eig(kj,mj)))/(2*pi)'];
end

sw = [sw;sw];
LANDA1 = -1;
kj = AK+LANDA1*BK+LANDA1*BK';
mj = AM+LANDA1*BM+LANDA1*BM';
sw = [sw;sqrt(real(eig(kj,mj)))/(2*pi)'];
LANDA1 = 1;
kj = AK+LANDA1*BK+LANDA1*BK';
mj = AM+LANDA1*BM+LANDA1*BM';
sw = [sw;sqrt(real(eig(kj,mj)))/(2*pi)'];
Fr = sort(sw);
Frequencies = Fr(1:5,1);
penalty = 0;
penalty = penalty+max(0,(1-(Frequencies(1,1)/minfrq(1,1))))+max(0,(1-(Frequencies(3,1)/minfrq(1,2))));
e1 = 2.5;
e2 = 0.5;
e3 = 0.5;

if se1>0.95*se2
    e1 = 200;
end

ppenalty = (1+e1*penalty)^(e2+e3*(((se1-1)/(se2-1))));
PFIT = W*ppenalty;

%% SOL
sol.Penalized_weight=PFIT;% kg
sol.weight=W;
sol.penalty=penalty;
sol.penalty_coefficient=ppenalty;
sol.A=AA'.*1e4;% cm^2
sol.x=X;
sol.omega=Frequencies;
sol.e1=e1;
sol.e2=e2;
sol.e3=e3;
sol.connectivity_1410bar=B;
sol.Coordinate_1410bar=nodes;
sol.load_1410bar=f;
end