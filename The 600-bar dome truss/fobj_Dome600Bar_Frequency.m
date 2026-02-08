function [W,PFIT,sol] = fobj_Dome600Bar_Frequency(X,se1,se2)

E = 20e10;
Density = 7850;
Mass = 100;
minfrq = [5 7];
FF = zeros(81, 81);
MM = zeros(81, 81);
W = 0; 
AA = [X';X'];
NUM_S = 24;
NUM_NODE = 9;
THETA = 2*pi/NUM_S;
RADIUS = [1 1 3 5 7 9 11 13 14];
Z_COORDS = [7 7.5 7.25 6.75 6 5 3.5 1.5 0];

nodes = [];
for ns = 1:3
    nodes = [nodes; [RADIUS * cos((ns - 1) * THETA); RADIUS * sin((ns - 1) * THETA); Z_COORDS]'];
end

B = [];
for ns = 1:2
    firstNodeCur = (ns - 1) * NUM_NODE + 1;
    firstNodeNext = firstNodeCur + NUM_NODE;
    B = [B; [firstNodeCur:firstNodeCur+NUM_NODE-2;firstNodeCur+1:firstNodeCur+NUM_NODE-1]'];
    B = [B; [firstNodeCur, firstNodeCur + 2]];
    B = [B; [firstNodeCur:firstNodeCur+NUM_NODE-2;firstNodeNext:firstNodeNext+NUM_NODE-2]'];
    B = [B; [firstNodeCur+2:firstNodeCur+NUM_NODE-1;firstNodeNext+1:firstNodeNext+NUM_NODE-2]'];
    B = [B; [firstNodeCur, firstNodeNext + 1]];
end

Constr = zeros(3, 27);
for j = 1:3
    Constr(:, 9*j) = 1;
end

Ur = 1-Constr;
f = find(Ur);

for i = 1:50
    XTOOL = nodes(B(i, 2), 1) - nodes(B(i, 1), 1);
    YTOOL = nodes(B(i, 2), 2) - nodes(B(i, 1), 2);
    ZTOOL = nodes(B(i, 2), 3) - nodes(B(i, 1), 3);
    ELEMTOOL = sqrt(XTOOL^2 + YTOOL^2 + ZTOOL^2);
    roi = sqrt(nodes(B(i, 1), 1)^2 + nodes(B(i, 1), 2)^2);
    roj = sqrt(nodes(B(i,2),1)^2+nodes(B(i,2),2)^2);
    loi = nodes(B(i,1),1)/roi;
    moi = nodes(B(i,1),2)/roi;
    loj = nodes(B(i,2),1)/roj;
    moj = nodes(B(i,2),2)/roj;
    Roi = [loi -moi 0;moi loi 0;0 0 1];
    Roj = [loj -moj 0;moj loj 0;0 0 1];
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
    W = W+12*Weight;
    
    SE = EAL*[se -se;-se se];
    K_cylindrical_coordinate = R'*SE*R;
    elementdof = [I2+1 I2+2 I2+3 H2+1 H2+2 H2+3];
    FF(elementdof,elementdof) = FF(elementdof,elementdof)+K_cylindrical_coordinate;
    ME = APL6*[2 0 0 1 0 0;0 2 0 0 1 0;0 0 2 0 0 1;1 0 0 2 0 0;0 1 0 0 2 0;0 0 1 0 0 2];
    M_cylindrical_coordinate = R'*ME*R;
    MM(elementdof,elementdof) = MM(elementdof,elementdof)+M_cylindrical_coordinate;
end

for j = 1:27
    MM(3*j-2,3*j-2) = MM(3*j-2,3*j-2)+Mass;
    MM(3*j-1,3*j-1) = MM(3*j-1,3*j-1)+Mass;
    MM(3*j,3*j) = MM(3*j,3*j)+Mass;
end

MM = MM(f,f);
FF = FF(f,f);
AK = FF(25:48,25:48);
AM = MM(25:48,25:48);
BK = FF(1:24,25:48);
BM = MM(1:24,25:48);
sw = [];

for J = 1:11
    LANDA1 = complex(cos(2*J*pi/24),sin(2*J*pi/24));
    LANDA2 = complex(cos(2*J*pi/24),-sin(2*J*pi/24));
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
Frequencies = Fr(1:6,1);
penalty = 0;
penalty = penalty+max(0,(1-(Frequencies(1,1)/minfrq(1,1))))+max(0,(1-(Frequencies(3,1)/minfrq(1,2))));
if se1<=0.95*se2
e1 = 1;
e2 = 1.5;
e3 = 4.5;
else
e1 = 1;
e2 = 1.5;
e3 = 15;
end
ppenalty = (1+e1*penalty)^(e2+e3*(((se1-1)/(se2-1))));
PFIT = W*ppenalty;

%% SOL
sol.Penalized_weight = PFIT;
sol.weight = W;
sol.penalty = penalty;
sol.penalty_coefficient = ppenalty;
sol.A = AA'.*1e4;% cm^2
sol.x = X;
sol.omega = Frequencies;
sol.e1 = e1;
sol.e2 = e2;
sol.e3 = e3;
sol.connectivity_600bar = B;
sol.Coordinate_600bar = nodes;
sol.load_600bar = f;

end