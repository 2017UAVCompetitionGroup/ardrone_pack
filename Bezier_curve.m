clear;
clc;

%输入测试用给定点矩阵%
k = 1:0.5:20;
T0 = 10;
Velocity_mam = 0.2;
X=k.^2/10;
Y=(k.^(2/3) + k.^(1/3))/10;
P=cat(1,X,Y);
%plot(X,Y);
%通过路径点反求控制点%
Q(:,1)=P(:,1);
Q(:,39)=P(:,39);
for i = 2:1:38
    Q(:,i)=( P(:,i-1) + 4*P(:,i) + P(:,i+1))/6/10;     %这里除以10为测试数据实际要去除
end
Q;

%通过控制点拟合B样条曲线并参数化,计算总路径长度与总用时%
Length_of_Curve = 0;
for i = 1:1:36
   syms t
   Bx_= 3*(-Q(1,i) + 3*Q(1,i+1) - 3*Q(1,i+2) + Q(1,i+3)).*t.^2 + 6*(Q(1,i) - 2*Q(1,i+1) + Q(1,i+2)).*t + 3*(-Q(1,i) + Q(1,i+1) + Q(1,i+2));
   By_= 3*(-Q(2,i) + 3*Q(2,i+1) - 3*Q(2,i+2) + Q(2,i+3)).*t.^2 + 6*(Q(2,i) - 2*Q(2,i+1) + Q(2,i+2)).*t + 3*(-Q(2,i) + Q(2,i+1) + Q(2,i+2));
   BL_= ((Bx_).^2 + (By_).^2).^(1/2);
   BL = int(BL_, 0, 1);
   Length_of_Curve = Length_of_Curve + BL;
end
double(Length_of_Curve)
T =double(Length_of_Curve / Velocity_mam - T0);

%自定义时间间隔各时刻飞机x，y方向分速度%
internal = 0.2;%自定义间隔时刻采样点
Terminal = double((T+2*T0)/internal);
for k = 1:1:Terminal
    L(1,k) = ( k-1) *internal;
    Velocity_Rate(1,k) = ( k-1) *internal;
    if (k-1)*internal < T0
       L(2,k) = ((k-1)*internal)^2*Velocity_mam/2/T0 ;
       Velocity_Rate(1,k) = Velocity_mam/T0 *( k-1) * internal;
      else if (k-1)*internal >= T0 && (k-1)*internal < T + T0
              L(2,k) = Velocity_mam*T0/2 + Velocity_mam*((k-1)*internal-T0);
              Velocity_Rate(1,k) = Velocity_mam;
        else if (k-1)*internal >= T + T0
                L(2,k) = (T + T0)*Velocity_mam - (T + 2*T0 - (k-1)*internal)*Velocity_mam/2/T0;
                Velocity_Rate(1,k) =Velocity_mam - Velocity_mam/T0 *(( k-1) * internal - T - T0);
            end
          end
    end
end
L;
Length_of_Curve = 0;

for n=1:Terminal
    while i < 37
        syms t
        Bx_= 3*(-Q(1,i) + 3*Q(1,i+1) - 3*Q(1,i+2) + Q(1,i+3)).*t.^2 + 6*(Q(1,i) - 2*Q(1,i+1) + Q(1,i+2)).*t + 3*(-Q(1,i) + Q(1,i+1) + Q(1,i+2));
        By_= 3*(-Q(2,i) + 3*Q(2,i+1) - 3*Q(2,i+2) + Q(2,i+3)).*t.^2 + 6*(Q(2,i) - 2*Q(2,i+1) + Q(2,i+2)).*t + 3*(-Q(2,i) + Q(2,i+1) + Q(2,i+2));
        BL_= ((Bx_).^2 + (By_).^2).^(1/2);
        BL = int(BL_, 0, 1);
        Length_of_Curve = Length_of_Curve + BL;
        
        if double(Length_of_Curve) >= L(2,n)
            syms v
            pos = solve(int(BL_,0,v)-double(L(2,n))+double(Length_of_Curve) + BL == 0 ,v);
            Dx_pos = 3*(-Q(1,i) + 3*Q(1,i+1) - 3*Q(1,i+2) + Q(1,i+3)).*double(pos).^2 + 6*(Q(1,i) - 2*Q(1,i+1) + Q(1,i+2)).*double(pos) + 3*(-Q(1,i) + Q(1,i+1) + Q(1,i+2));
            Dy_pos = 3*(-Q(2,i) + 3*Q(2,i+1) - 3*Q(2,i+2) + Q(2,i+3)).*double(pos).^2 + 6*(Q(2,i) - 2*Q(2,i+1) + Q(2,i+2)).*double(pos) + 3*(-Q(2,i) + Q(2,i+1) + Q(2,i+2));
            Velocity(1,n) = Velocity_Rate(1,n) * cos(atan(Dy_pos/Dx_pos));
            Velocity(2,n) = Velocity_Rate(1,n) * sin(atan(Dy_pos/Dx_pos));
                               
            break;
        else
            i = i+1;
        end
    end
end 
