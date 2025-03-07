# zikuzyuchisuitei-hampu-up
hanpu up

clear; clc;

%% ===============================
%% 1. シミュレーション設定
%% ===============================
m1 = 1900;  m2 = 50;  m3 = 50;   % 車体、車軸の質量
k1 = 33e3;  k2 = 33e3;          % 前後サスペンション剛性 [N/m]
k3 = 260e3; k4 = 260e3;         % タイヤ剛性 [N/m]
c1 = 300;   c2 = 300;           % サスペンション減衰 [Ns/m]
l1 = 1.34;  l2 = 1.46;          % 重心から前・後車軸までの距離 [m]
J  = 3000;                      % 車体慣性モーメント [kg*m^2]

% 走行条件
v = 20;           % 車両速度 [m/s]
t_end = 10;       % シミュレーション時間 [s]
dt = 0.001;      
t = 0:dt:t_end;

% 状態変数の初期化
x1 = zeros(size(t));  dx1    = zeros(size(t));
theta = zeros(size(t)); dtheta = zeros(size(t));
x2 = zeros(size(t));  dx2    = zeros(size(t));
x3 = zeros(size(t));  dx3    = zeros(size(t));

Ff = zeros(size(t));  
Fr = zeros(size(t));
Fy1_array = zeros(size(t)); 
Fy2_array = zeros(size(t));

%% ===============================
%% 2. 台形型ハンプの設定
%% ===============================
% ハンプ開始位置を s=0 とする
hump_start      = 0;    % ハンプ開始 [m]
ramp_up_length  = 10;   % 上昇部の長さ [m]
top_length      = 15;   % 上辺の長さ [m]
ramp_down_length= 10;   % 下降部の長さ [m]
hump_height     = 0.20;  % ハンプの最大高さ [m]

% 台形型ハンプの路面高さ関数
road_hump = @(s) ...
    (s < hump_start).*0 + ...
    (s >= hump_start & s < hump_start + ramp_up_length) .* ...
        ((s - hump_start)/ramp_up_length * hump_height) + ...
    (s >= hump_start + ramp_up_length & ...
     s < hump_start + ramp_up_length + top_length) .* ...
        hump_height + ...
    (s >= hump_start + ramp_up_length + top_length & ...
     s < hump_start + ramp_up_length + top_length + ramp_down_length) .* ...
        (hump_height - ...
         ( (s - (hump_start + ramp_up_length + top_length)) / ramp_down_length * hump_height)) + ...
    (s >= hump_start + ramp_up_length + top_length + ramp_down_length).*0;

%% ===============================
%% 3. シミュレーションループ
%% ===============================
for i = 1:length(t)-1
    % 車両の現在進行距離 s [m]
    s = v * t(i);
    
    % 路面高さ (前後同一と仮定)
    road_height = road_hump(s);
    y1 = road_height;
    y2 = road_height;
    
    % サスペンション力(前: Ff, 後: Fr)
    Ff(i) = -k1 * ((x1(i) - l1*theta(i)) - x2(i)) ...
            - c1 * ((dx1(i) - l1*dtheta(i)) - dx2(i));
    Fr(i) = -k2 * ((x1(i) + l2*theta(i)) - x3(i)) ...
            - c2 * ((dx1(i) + l2*dtheta(i)) - dx3(i));
    
    % タイヤ力 (前: Fy1, 後: Fy2)
    Fy1 = -k3 * (x2(i) - y1);
    Fy2 = -k4 * (x3(i) - y2);
    Fy1_array(i) = Fy1;
    Fy2_array(i) = Fy2;
    
    % 車体 (m1, J) の運動方程式
    ddx1    = 2*(Ff(i) + Fr(i)) / m1;
    ddtheta = 2*(-Ff(i)*l1 + Fr(i)*l2) / J;
    
    % 前軸 (m2), 後軸 (m3) の運動方程式
    ddx2 = ( -k3*(x2(i)-y1) - Ff(i) ) / m2;
    ddx3 = ( -k4*(x3(i)-y2) - Fr(i) ) / m3;
    
    % オイラー法で状態更新
    dx1(i+1)   = dx1(i)    + ddx1*dt;
    x1(i+1)    = x1(i)     + dx1(i)*dt;
    
    dtheta(i+1)= dtheta(i) + ddtheta*dt;
    theta(i+1) = theta(i)  + dtheta(i)*dt;
    
    dx2(i+1)   = dx2(i)    + ddx2*dt;
    x2(i+1)    = x2(i)     + dx2(i)*dt;
    
    dx3(i+1)   = dx3(i)    + ddx3*dt;
    x3(i+1)    = x3(i)     + dx3(i)*dt;
end

% 最終ステップの力を補完
Ff(end) = -k1*((x1(end)-l1*theta(end)) - x2(end)) ...
          - c1*((dx1(end)-l1*dtheta(end)) - dx2(end));
Fr(end) = -k2*((x1(end)+l2*theta(end)) - x3(end)) ...
          - c2*((dx1(end)+l2*dtheta(end)) - dx3(end));

final_s = v * t(end);
final_road = road_hump(final_s);
Fy1_array(end) = -k3*( x2(end) - final_road );
Fy2_array(end) = -k4*( x3(end) - final_road );

%% ===============================
%% 4. タイヤ力の可視化 (全区間)
%% ===============================
figure('Color','w');
plot(t, Fy1_array, 'b', 'LineWidth',1.5); hold on;
plot(t, Fy2_array, 'r', 'LineWidth',1.5);
xlabel('時間 [s]', 'FontSize', 16);
ylabel('垂直方向力 [N]', 'FontSize', 16);
legend({'前車軸タイヤ力','後車軸タイヤ力'}, 'FontSize', 16, 'Location','best');
grid on; set(gca,'FontSize',16);

%% ===============================
%% 5. ハンプ上辺区間のみのデータを抽出
%% ===============================
% 走行距離 s = v * t
s_all = v * t;


top_start = 5; 
top_end   = 10; 

% インデックス抽出
top_idx = find(s_all >= top_start & s_all <= top_end);

% 上辺区間のみのデータ
t_top   = t(top_idx);
Fy1_top = Fy1_array(top_idx);
Fy2_top = Fy2_array(top_idx);

% 上辺の路面高さ 
u_top   = arrayfun(@(ss) road_hump(ss), s_all(top_idx));

%% ===============================
%% 6. ARXモデルによる推定 (前車軸, 上辺区間のみ)
%% ===============================
% iddata形式に変換 (入力: 台形ハンプ高さ, 出力: 前車軸タイヤ力)
data_top_front = iddata(Fy1_top', u_top', dt);

% ARXモデル
na = 5; nb = 0; nk = 1;
model_arx_front = arx(data_top_front, [na nb nk]);

% ARXモデル出力
[simOut_front, fitVal_front, ~] = compare(data_top_front, model_arx_front);
y_pred_top_front = simOut_front.OutputData;

% 誤差評価
err_top_front   = Fy1_top' - y_pred_top_front;
RMSE_top_front  = sqrt(mean(err_top_front.^2));
maxErr_top_front= max(abs(err_top_front));

%% 結果表示 (前車軸)
fprintf('=== ハンプ上辺のみ ARXモデル(前車軸) 適用結果 ===\n');
disp(model_arx_front);
fprintf('Fit率         = %.2f %%\n', fitVal_front);
fprintf('RMSE          = %.4f [N]\n', RMSE_top_front);
fprintf('最大ピーク誤差 = %.4f [N]\n', maxErr_top_front);

%% 前車軸タイヤ力の予測比較プロット
figure('Color','w');
plot(t_top, Fy1_top, 'b-', 'LineWidth',1.5); hold on;
plot(t_top, y_pred_top_front, 'm--', 'LineWidth',2.5);
xlabel('時間 [s]', 'FontSize', 16);
ylabel('垂直方向力 [N]', 'FontSize', 16);
legend({'モデルによる実測出力','モデルによる出力予測'}, ...
       'FontSize', 16, 'Location','best');
grid on; set(gca,'FontSize',16);

%% ===============================
%% 7. ARXモデルによる推定 (後車軸, 上辺区間のみ)
%% ===============================
% 後車軸用に iddata を作成
data_top_rear = iddata(Fy2_top', u_top', dt);

% 前車軸と同じ次数で ARX モデルを作成
model_arx_rear = arx(data_top_rear, [na nb nk]);

% ARXモデル出力
[simOut_rear, fitVal_rear, ~] = compare(data_top_rear, model_arx_rear);
y_pred_top_rear = simOut_rear.OutputData;

% 誤差評価
err_top_rear   = Fy2_top' - y_pred_top_rear;
RMSE_top_rear  = sqrt(mean(err_top_rear.^2));
maxErr_top_rear= max(abs(err_top_rear));

%% 結果表示 (後車軸) 
fprintf('\n=== ハンプ上辺のみ ARXモデル(後車軸) 適用結果 ===\n');
disp(model_arx_rear);
fprintf('Fit率         = %.2f %%\n', fitVal_rear);
fprintf('RMSE          = %.4f [N]\n', RMSE_top_rear);
fprintf('最大ピーク誤差 = %.4f [N]\n', maxErr_top_rear);
