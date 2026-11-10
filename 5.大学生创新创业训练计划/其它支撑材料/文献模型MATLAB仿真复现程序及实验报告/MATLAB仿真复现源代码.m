%1.列出已知信息
%1.1 时间
dt=0.01;
time=0:dt:100; %下面仿真进行循环时需要设t即每一小份时间，故将这里改成time，且在第一部分中无法定义延迟，因为延迟需要用到t，只能在循环中定义
n=4; %4个无人机
N=length(time); %将0到100s的时间分成10000份

%1.2 编队拓扑
A=[0 1 0 0;
    0 0 1 0;
    0 0 0 1;
    1 0 0 0]; %邻接矩阵
%r_w=[-50 50;
%    50 50 ;
%    50 -50;
%   -50 -50]; %各无人机相对中心位置的偏移量，起先我设置成四行二列，但与下面设置的矩阵不匹配
r_w=[-50 50 50 -50;
    50 50 -50 -50];
alpha=zeros(n,1); %设置领导者，先设n行1列的零矩阵
alpha(1)=1; %由于我们设UAV1为领导者，故除了第一个为1，其余都为0

%%2.初始条件
r_h=zeros(2,n,N); %实际历史位置,先初始化为一个三维的矩阵（演示一下）
r_h_d=zeros(2,n,N); %期望历史位置
r_c_d=[0;0]; %路径中心初始位置，注意与上面的偏移量形式相同

for i=1:n
    r_h_d(:,i,1)=r_c_d+r_w(:,i);%中心位置加上偏移量，即每个UAV（i）初始位置 
    %这就是为什么上面的偏移量不按以往的习惯设置成四行二列的，r_h_d对于每个时刻是二行四列的
end

%%3.主仿真开始，即开始循环
for k=1:N  %(N-1还是N)对于时间进行循环
    t=time(k); %第k份时的具体时间是第几时刻

    %3.1期望路径计算
    v_0=10; %初始时刻的速度
    if t<=50
        theta=60*(pi/180); %在前50s，无人机前进的角度是60°，这里转换成弧度制，下面同理
    else
        theta=-60*(pi/180);
    end
    v_d=[v_0*cos(theta);
        v_0*sin(theta)]; %前进过程中的速度，注意这里也是两行一列（向量格式）
    %下面更新中心期望位置
    r_c_d=r_c_d+v_d*dt;
    %然后就可以计算每个无人机的期望位置与速度
    r_d=r_c_d+r_w;
    v_d=repmat(v_d,1,n); %%复制矩阵v_d，行上复制1次，列上复制4次(4个无人机)，使每个无人机都使用相同的期望速度向量。
    %由控制器u(:,i) = v_d(:,i) + sum - alpha_term，需要为每个无人机(i)指定一个速度向量v_d(:,i)，所以这里需要复制
    u=zeros(2,n); %式8的控制器

    %3.2时变延迟函数定义
    %通信延迟
    tau12=0.3*abs(sin(t))+0.1;
    tau23=0.2*abs(sin(t))+0.2;
    tau34=0.1*abs(sin(t))+0.3;
    tau14=0.1*abs(sin(t))+0.2;
    tau=[NaN tau12 NaN NaN;
        NaN NaN tau23 NaN;
        NaN NaN NaN tau34;
        tau14 NaN NaN NaN]; %设置成矩阵，方便下面索引（i，j），自身延迟同理
    %自身延迟
    T=[0.20*abs(sin(t));
        0.30*abs(sin(t));
        0.25*abs(sin(t));
        0.15*abs(sin(t))];

    for i=1:n
        %3.3求和项的计算
        sum=[0;0]; %初始化求和项
        neighbors=find(A(i,:)==1); %返回A中每行值为1的列索引，即找到i的邻居
        for j=neighbors
            %下面获得延迟状态（在历史记录中找最近时间点）
            k_delay=max(1,k-round(tau(i,j)/dt)); 
            %将延迟换成时间步数后进行四舍五入（round），之后用当前步数减去延迟，计算延迟后的时间步索引，max保证索引不小于1
            %然后获得邻居j在通信延迟时刻的实际位置与期望位置
            r_j_delay=r_h(:,j,k_delay);
            r_d_j_delay=r_h_d(:,j,k_delay);

            %对自身延迟重复上面的操作
            k_delay_s=max(1,k-round(T(i)/dt));
            r_i_delay=r_h(:,i,k_delay_s);
            r_d_i_delay=r_h_d(:,i,k_delay_s);

            %按照公式8，计算一致性协议的误差
            error_j=r_j_delay-r_d_j_delay;
            error_i=r_i_delay-r_d_i_delay;

            sum=sum+A(i,j)*(error_j-error_i); %更新求和项
        end
        %3.4领导者跟随项
        k_delay_s=max(1,k-round(T(i)/dt));
        r_i_delay=r_h(:,i,k_delay_s);
        r_d_i_delay=r_h_d(:,i,k_delay_s);
        error_i_a=r_i_delay-r_d_i_delay;
        alpha_term=alpha(i)*error_i_a;

        %3.5最终控制器（公式8）
        u(:,i)=v_d(:,i)+sum-alpha_term;
    end

    %3.6更新状态
    r_c=r_h(:,:,k); %获取当前位置
    r_n=r_c+u(:,:)*dt; %更新位置r(k+1)=r(k)+u*dt,其实这里的u既然用（:,:),可以直接用u，不用写后面那些

    r_h(:,:,k+1)=r_n;
    r_h_d(:,:,k+1)=r_d;
end

%%4.画图
% 提取历史数据便于绘图
x_hist = squeeze(r_h(1, :, :)); % 所有无人机的x坐标 (4×N)
y_hist = squeeze(r_h(2, :, :));
xd_hist = squeeze(r_h_d(1, :, :)); % 所有无人机的期望x坐标
yd_hist = squeeze(r_h_d(2, :, :));

%编队轨迹
figure('Name', '图3：折线路径编队飞行');
hold on;
colors = ['r', 'g', 'b', 'm']; % 无人机颜色 (修改为与论文一致的颜色)
for i = 1:n
    plot(x_hist(i, :), y_hist(i, :), colors(i), 'LineWidth', 1.5);
end
%定期绘制期望编队形状
plot_interval = 10 / dt; %每10秒绘制一次
for k = 1:plot_interval:N
    %直接使用存储的期望位置
    plot(xd_hist(:, k), yd_hist(:, k), ...
         '--ks', 'MarkerFaceColor', 'k', 'LineWidth', 0.5);
    
    %绘制中心点
    center_x = mean(xd_hist(:, k));
    center_y = mean(yd_hist(:, k));
    plot(center_x, center_y, 'k*');
end

%标记路径转折点（50秒处）
k50 = find(time >= 50, 1); %找到50秒对应的索引
turn_point = [mean(xd_hist(:, k50)), mean(yd_hist(:, k50))];
plot(turn_point(1), turn_point(2), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
text(turn_point(1)+50, turn_point(2), '转折点 (t=50s)', 'FontSize', 10);

title('图3：折线路径编队飞行');
xlabel('x(米)');
ylabel('y(米)');
legend('UAV1', 'UAV2', 'UAV3', 'UAV4', '期望编队', '中心点', 'Location', 'best');
axis equal;
grid on;
hold off;

%速度曲线
figure('Name', '图4：速度曲线');
hold on;

%计算实际速度大小 (使用实际位置数据)
speed_hist = zeros(n, N);
for i = 1:n
    for k = 1:N-1
        dx = x_hist(i, k+1) - x_hist(i, k);
        dy = y_hist(i, k+1) - y_hist(i, k);
        speed_hist(i, k) = sqrt(dx^2 + dy^2) / dt;
    end
    %最后一时刻使用前一时刻的速度
    speed_hist(i, end) = speed_hist(i, end-1);
end

%绘制速度曲线
for i = 1:n
    plot(time, speed_hist(i, :), colors(i), 'LineWidth', 1);
end

%标记转折点
line([50,50], [0, 15], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1); %修改y范围
text(52, 14, 't=50s', 'FontSize', 10); %调整文本位置

title('图4：速度曲线');
xlabel('时间(秒)');
ylabel('速度(米/秒)');
legend('UAV1', 'UAV2', 'UAV3', 'UAV4', 'Location', 'best');
grid on;
ylim([0, 150]); %调整y轴范围匹配实际速度
hold off;

%跟踪误差
figure('Name', '图5：跟踪误差');
hold on;
error_hist = zeros(n, N);
for i = 1:n
    for k = 1:N
        dx = x_hist(i, k) - xd_hist(i, k);
        dy = y_hist(i, k) - yd_hist(i, k);
        error_hist(i, k) = sqrt(dx^2 + dy^2);
    end
end
%绘制误差曲线
for i = 1:n
    plot(time, error_hist(i, :), colors(i), 'LineWidth', 1.5);
end
%标记转折点
line([50,50], [0, max(error_hist(:))], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
text(52, max(error_hist(:))-5, 't=50s', 'FontSize', 10);

title('图5：跟踪误差');
xlabel('时间(秒)');
ylabel('位置误差(米)');
legend('UAV1', 'UAV2', 'UAV3', 'UAV4', 'Location', 'best');
grid on;
hold off;