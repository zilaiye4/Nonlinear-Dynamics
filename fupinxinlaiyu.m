% 定义参数
clear;clc;close all;
format long;
epsilon = 0.1;
eta = 0.02;
gama=0.02;
w0 = 2;
omega0=2;
alpha =4;
K = 1;
p = 0.6;
f = 1;
F2 =30;
EE=[];
FF=[];
EEE=[];
FFF=[];
EEE1=[];
EEE2=[];
EEE3=[];
EEE4=[];
EEE5=[];
EEE6=[];
FFF2=[];
FFF3=[];
FFF4=[];
SIGMA=[];

% 定义频率范围和步长
sigma_values = linspace(0,4 , 99); % 调整范围和点数以获得更好的分辨率
omega1_values = omega0 + epsilon * sigma_values;
omega2_values = 3 * omega1_values;

% 初始化振幅和相位数组
a_values = zeros(size(sigma_values));
phi_values = zeros(size(sigma_values));
for  i=1:length(sigma_values)
     
            SIGMA=[SIGMA;sigma_values(i)];                                       %将对应的W，存入xx
           
end
% 计算稳态振幅和相位
for i = 1:length(sigma_values)
    stabA=[];
unstabA=[];
uuu=1;
uuuu=1;
u1=1;
    sigma = sigma_values(i)
    omega1 = omega1_values(i);
    omega2 = omega2_values(i);
    
    B = F2 / (2 * (omega0^2 - omega2^2));
    
    % 定义方程函数（输入为 [a; phi]）
    fun = @(x) [-f*sin(x(2))/(2*w0)+3*B*alpha*x(1)^2*sin(3*x(2))/(4*w0)+eta*x(1)/2-27*B^2*gama*x(1)*w0^2-3*gama*x(1)^3*w0^2/8+9*B*gama*x(1)^2*cos(3*x(2))*w0^2/4-K/2*x(1)*sin(p*pi/2)*w0^(p-1);  
        -sigma * x(1)-f*cos(x(2))/(2*w0)+3*B^2*alpha*x(1)/w0+3*alpha*x(1)^3/(8*w0)+3*B*alpha*x(1)^2*cos(3*x(2))/(8*w0)-9*B*gama*x(1)^2*w0^2*sin(3*x(2))/4+K/2*x(1)*cos(p*pi/2)*w0^(p-1) ];
    a1=0.1:0.1:2;
b1=0:0.5:2*pi;
    for ii=1:length(a1)
a=a1(ii);
for iii=1:length(b1)
    b=b1(iii);
    % 初始猜测值 [a; phi]
    x_guess = [a; b]; % 调整初始猜测值以匹配实际物理意义
    mmm=nan;
 
   options = optimoptions('fsolve', ...
    'FunctionTolerance', 1e-10, ...  % 函数值容差
    'StepTolerance', 1e-10, ...      % 步长容差
    'OptimalityTolerance', 1e-10, ... % 最优性容差
    'MaxIterations', 1000, 'MaxFunctionEvaluations',5000,'Display', 'off');             % 显示迭代过程
    sol= myve(fun, x_guess, options);
    if sol(1)>0&&sol(2)>=0&&sol(2)<=2*pi
    a_values(i) = sol(1);
    phi_values(i) = sol(2);
    mmm=sol(1);
    end
   rho11=eta/2 - 27*B^2*gama*w0^2 - (9*sol(1)^2*gama*w0^2)/8 - (K*w0^(p - 1)*sin((pi*p)/2))/2 + (9*B*sol(1)*gama*w0^2*cos(3*sol(2)))/2 + (3*B*sol(1)*alpha*sin(3*sol(2)))/(2*w0);
rho12=(9*B*sol(1)^2*alpha*cos(3*sol(2)))/(4*w0) - (f*cos(sol(2)))/(2*w0) - (27*B*sol(1)^2*gama*w0^2*sin(3*sol(2)))/4;
rho21=(3*B^2*alpha)/w0 - sigma + (9*sol(1)^2*alpha)/(8*w0) + (K*w0^(p - 1)*cos((pi*p)/2))/2 + (3*B*sol(1)*alpha*cos(3*sol(2)))/(4*w0) - (9*B*sol(1)*gama*w0^2*sin(3*sol(2)))/2;
rho22= (f*sin(sol(2)))/(2*w0) - (27*B*sol(1)^2*gama*w0^2*cos(3*sol(2)))/4 - (9*B*sol(1)^2*alpha*sin(3*sol(2)))/(8*w0);
eta1=-rho11-rho22;
eta2=rho11*rho22-rho12*rho21;
   if  eta1>1e-4&&eta2>1e-4
stabA(uuu)=mmm;
uuu=uuu+1;
    elseif (eta1<-1e-1)||(eta2<-1e-1)
   unstabA(uuuu)=mmm;
uuuu=uuuu+1;
end
end   
    

    end
    
    C=[];
D=[];
C=uniquetol(stabA,1e-4);
D=uniquetol(unstabA,1e-4);
C=C(~isnan(C));
D=D(~isnan(D));
C
D;
if length(C)>0
EE(i,:)=C(1);
 EEE=[EEE;[SIGMA(i),EE(i)]];
end
if length(C)>1
EE1(i,:)=C(2);
EEE1=[EEE1;[SIGMA(i),EE1(i)]];
end
if length(C)>2
EE2(i,:)=C(3);
EEE2=[EEE2;[SIGMA(i),EE2(i)]];
end
if length(C)>3
EE3(i,:)=C(4);
EEE3=[EEE3;[SIGMA(i),EE3(i)]];
end
if length(C)>4
EE4(i,:)=C(5);
EEE4=[EEE4;[SIGMA(i),EE4(i)]];
end
if length(C)>5
EE5(i,:)=C(6);
EEE5=[EEE5;[SIGMA(i),EE5(i)]];
end
if length(C)>6
EE6(i,:)=C(7);
EEE6=[EEE6;[SIGMA(i),EE6(i)]];
end
if length(D)>0
FF(i,:)=D(1);
FFF=[FFF;[SIGMA(i),FF(i)]];
end
if length(D)>1
FF1(i,:)=D(2);
FFF2=[FFF2;[SIGMA(i),FF1(i)]];
end
if length(D)>1
FF1(i,:)=D(2);
FFF2=[FFF2;[SIGMA(i),FF1(i)]];
end
if length(D)>2
FF2(i,:)=D(3);
FFF3=[FFF3;[SIGMA(i),FF2(i)]];
end
if length(D)>3
FF3(i,:)=D(4);
FFF4=[FFF4;[SIGMA(i),FF3(i)]];
end
  end
figure(1)
if length(FFF)
plot(FFF(:,1),FFF(:,2),'k*','MarkerSize',4);
hold on;
end
plot(EEE(:,1),EEE(:,2),'ko','MarkerSize',4);
hold on;
if length(EEE1)
plot(EEE1(:,1),EEE1(:,2),'ko','MarkerSize',4);
hold on
end
if length(EEE2)
plot(EEE2(:,1),EEE2(:,2),'ko','MarkerSize',4);
hold on
end
if length(EEE3)
plot(EEE3(:,1),EEE3(:,2),'ko','MarkerSize',4);
hold on
end
if length(EEE4)
plot(EEE4(:,1),EEE4(:,2),'ko','MarkerSize',4);
hold on
end
if length(EEE5)
plot(EEE5(:,1),EEE5(:,2),'ko','MarkerSize',4);
hold on
end
if length(EEE6)
plot(EEE6(:,1),EEE6(:,2),'ko','MarkerSize',4);
hold on
end
if length(FFF2)
plot(FFF2(:,1),FFF2(:,2),'k*','MarkerSize',4);
hold on;
end
if length(FFF3)
plot(FFF3(:,1),FFF3(:,2),'k*','MarkerSize',4);
hold on;
end
if length(FFF4)
plot(FFF4(:,1),FFF4(:,2),'k*','MarkerSize',4);
hold on;
end

legend('unstable solution','stable solution')
  xlabel('$\sigma$','interpreter','latex','fontsize',18)
 ylabel('$a$', 'interpreter', 'latex','fontsize',18)
set(gca, 'Fontname' ,'Times New Roman','FontSize',18,'LineWidth',1);

function [x, fval, exitflag, output] = myve(fun, x0, options)
   
    %   fun     - 函数句柄，返回方程组值F(x)
    %   x0      - 初始点
    %   options - 选项结构体，包含容差、最大迭代次数等
    % 输出:
    %   x       - 求解结果
    %   fval    - 方程组在解处的值
    %   exitflag - 退出标志
    %   output  - 输出信息结构体
    
    % 默认参数设置
    if nargin < 3
        options = struct();
    end
    if ~isfield(options, 'TolX')
        options.TolX = 1e-6;
    end
    if ~isfield(options, 'TolFun')
        options.TolFun = 1e-6;
    end
    if ~isfield(options, 'MaxIter')
        options.MaxIter = 100;
    end
    
    % 初始化
    x = x0;
    iter = 0;
    exitflag = 0;
    output.iterations = 0;
    
    % 计算初始函数值
    F = feval(fun, x);
    fval = F;
    f_norm = norm(F);
    
    % 检查初始点是否已满足条件
    if f_norm < options.TolFun
        exitflag = 1;
        output.message = '初始点已满足收敛条件';
        return;
    end
    
    % 初始信赖域半径
    delta = 0.1;
    if length(x0) > 1
        delta = min(delta, norm(x0) * 0.1);
    end
    
    % 主迭代循环
    while iter < options.MaxIter
        iter = iter + 1;
        
        % 计算雅可比矩阵（使用数值微分）
        J = jacobian(fun, x);
        
        % 求解信赖域子问题
        s = trust_region_step(J, F, delta);
        
        % 计算新点和函数值
        x_new = x + s;
        F_new = feval(fun, x_new);
        f_new_norm = norm(F_new);
        
        % 计算实际下降和预测下降
        actual_reduction = f_norm^2 - f_new_norm^2;
        predicted_reduction = compute_predicted_reduction(J, F, s);
        
        % 计算rho
        if predicted_reduction > 0
            rho = actual_reduction / predicted_reduction;
        else
            rho = 0;
        end
        
        % 更新信赖域半径
        if rho > 0.75
            delta = min(2*delta, 1e6);  % 限制最大半径
        elseif rho < 0.25
            delta = 0.5*delta;
        end
        
        % 接受或拒绝步长
        if rho > 1e-4
            x = x_new;
            F = F_new;
            f_norm = f_new_norm;
        end
        
        % 检查收敛条件
        if f_norm < options.TolFun
            exitflag = 1;
            break;
        end
        
        if norm(s) < options.TolX*(norm(x) + options.TolX)
            exitflag = 2;
            break;
        end
    end
    
    % 输出信息
    output.iterations = iter;
    output.fval = f_norm;
    
    if exitflag == 0
        if iter == options.MaxIter
            exitflag = 0;
            output.message = '已达到最大迭代次数';
        else
            output.message = '未满足收敛条件';
        end
    else
        output.message = '成功收敛';
    end
end

function J = jacobian(fun, x)
    % 数值方法计算雅可比矩阵
    h = 1e-8;  % 步长
    n = length(x);
    F0 = feval(fun, x);
    m = length(F0);
    J = zeros(m, n);
    
    for i = 1:n
        x1 = x;
        x1(i) = x1(i) + h;
        F1 = feval(fun, x1);
        J(:, i) = (F1 - F0) / h;
    end
end

function s = trust_region_step(J, F, delta)
    % 求解信赖域子问题：min ||F + J*s||^2  subject to ||s|| <= delta
    A = J' * J;
    g = J' * F;
    
    % 无约束解
    s_uncon = -A\g;
    
    if norm(s_uncon) <= delta
        s = s_uncon;
        return;
    end
    
    % 约束解，使用二分法寻找最优lambda
    lambda = 0;
    lambda_max = 1e6;
    
    % 二分法迭代
    for i = 1:50
        lambda_mid = (lambda + lambda_max) / 2;
        s_mid = -(A + lambda_mid * eye(size(A))) \ g;
        norm_s = norm(s_mid);
        
        if norm_s < delta
            lambda_max = lambda_mid;
        else
            lambda = lambda_mid;
        end
        
        if abs(norm_s - delta) < 1e-8
            break;
        end
    end
    
    s = -(A + lambda_max * eye(size(A))) \ g;
end

function pred_red = compute_predicted_reduction(J, F, s)
    % 计算预测下降量
    pred_red = 2 * F' * J * s + s' * J' * J * s;
    pred_red = -pred_red;  % 确保为正值
end
