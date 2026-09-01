clc; clear; close all;

simTime = 5;
dt = 0.01;

kOS = 5; %overshoot penalty coefficient
Gs = tf(1, [1, 1, -6]);
Rtf = @(t) t.^2+sin(pi*t);


nParticles = 50;
steps = 20;
wMax = 0.5;    %momentum
wMin = 0.2;
c1 = 0.6;   %cognitive confidence
c2 = 0.6;   %group confidence

%initialize the random state
stateMax = [200, 200, 200]';
stateMin = [0, 0, 0]';

PState = stateMin + rand(3, nParticles).*(stateMax - stateMin);
v = zeros(size(PState));

timePState = zeros([3, nParticles, steps + 1]);
timePState(:, :, 1) = PState;

pbest = PState;
bestCost = ones(nParticles, 1)*inf;
gbestCost = inf;
for i = 1:steps
    w = wMax - (wMax - wMin)*i/steps;
    cost = zeros(nParticles, 1);
    parfor j = 1:nParticles
        Kp = PState(1, j); Kd = PState(2, j); Ki = PState(3, j);
        num = [Kd, Kp, Ki]; den = [1, 0];

        Cs = tf(num, den);
        YsRs = Gs*Cs/(1 + Gs*Cs);
        
        t = 0:dt:simTime;
        Rt = ones(size(t)).*Rtf(t);

        [y, tout] = lsim(YsRs, Rt, t);
        
        e = Rt' - y;
        cost(j) = trapz(tout,  tout.*abs(e));
        if ~all(e > 0)
            c = trapz(t(e < 0), abs(e(e < 0)).^2);
            cost(j) = cost(j) + c + kOS*min(e)^2;
        end

        if(cost(j) < bestCost(j))
            bestCost(j) = cost(j);
            pbest(:, j) = PState(:, j);
        end
    end

    if(min(cost) < gbestCost)
        [gbestCost, idx] = min(cost);
        gbest = PState(:, idx);
    end

    v = w*v + c1*(pbest - PState).*rand(size(PState)) + c2*(gbest - PState).*rand(size(PState));
    PState = PState + v;
    timePState(:, :, i+1) = PState;
end

disp("Global Best: ");
disp(gbest);

figure(1);
hold on; grid on; grid minor;
for i = 1:nParticles
    plot3(timePState(1, i, 1), timePState(2, i, 1), timePState(3, i, 1), "Marker", "s", "MarkerSize", 8, "LineWidth", 1.5);
    plot3(timePState(1, i, end), timePState(2, i, end), timePState(3, i, end), "Marker", "^", "MarkerSize", 8, "LineWidth", 1.5);

    plot3(squeeze(timePState(1, i, :)), squeeze(timePState(2, i, :)), squeeze(timePState(3, i, :)), "LineWidth", 1.5);
end
xlabel("K_p");
ylabel("K_d");
zlabel("K_i");
hold off;


t = 0:dt:simTime;
Rt = ones(size(t)).*Rtf(t);
figure(2);
hold on; grid on; grid minor;
Cs = tf([gbest(2), gbest(1), gbest(3)], [1, 0]);
YsRs = Gs*Cs/(1 + Gs*Cs);
[y, tout] = lsim(YsRs, Rt, t);
plot(tout, Rt', "LineStyle", ":", "LineWidth", 1.5);
plot(tout, y, "LineWidth", 1.5);
xlabel("Response");
ylabel("Time");
hold off;
