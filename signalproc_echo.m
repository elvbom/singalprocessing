close all;
clear;

%load signals
[s, fs_v]  = audioread('vocal.wav'); %far end signal
[x, fs_d] = audioread('drumloop.wav'); %near end signal

%make signals same length
x(length(x):length(s)) = 0;

%play signals
% sound(s, fs_v)
% sound(x, fs_d)

%generate echo, play and plot against near end signal
delay = fs_d*[0.04 0.08 0.1]; %fs*s
c = [0.5 0.3 0.1];

%echo paths - start at first delay, end when last one is heard
[y1, y2, y3] = deal(zeros(length(s) + delay(3), 1));
for i = delay(1)+1:length(s)+delay(1) 
    y1(i) = c(1)*x(i-delay(1));
end

for i = delay(2)+1:length(s)+delay(2)
    y2(i) = c(2)*x(i-delay(2));
end

for i = delay(3)+1:length(s)+delay(3)
    y3(i) = c(3)*x(i-delay(3));
end

%sum echos
y = y1+y2+y3;
%sound(y, fs_d);

figure;
subplot(2, 1, 1)
plot(x, 'g')
title('Near end signal');
xlabel('Time [s]');
ylabel('Amplitude');

subplot(2, 1, 2)
plot(y, 'b')
title('Echo');
xlabel('Time [s]');
ylabel('Amplitude');

%generate noise contamined far end signal
variance = 0.01;
v = sqrt(variance)*randn(length(s), 1);
u = v+s;
%sound(s, fs_d)
figure;
subplot(2, 1, 1)
plot(s, 'r')
title('Far end signal');
xlabel('Time [s]');
ylabel('Amplitude');

subplot(2, 1, 2)
plot(u, 'y')
title('Noise');
xlabel('Time [s]');
ylabel('Amplitude');

%make u same lenght as y, generate mixed signal sr
u(length(u):length(y)) = 0;
sr = y+u;
% sound(sr, fs_v);
figure;
plot(sr);
title('Mixed signal');
xlabel('Time [s]');
ylabel('Amplitude');

%estimate echo cofficients c using LMS and RLS

%LMS w/ known delay
M = length(c); %number of coefficents
my = 0.001; %step size
[c_lms, c_lms_help] = deal(zeros(M, 1));
c_lms_ev = zeros(M, length(x));
[y_lms, e_lms] = deal(zeros(length(x), 1));
x(length(x):length(sr)) = 0; %for whole loop to run
tic;
for i = 1:length(sr)
    for j = 1:M
        if i-delay(j) > 0
            c_lms_help(j) = x(i-delay(j));
        end
    end
    y_lms(i) = c_lms'*c_lms_help;
    e_lms(i) = sr(i)-y_lms(i);
    c_lms = c_lms+2*my*c_lms_help*e_lms(i);
    c_lms_ev(:,i) = c_lms;
end
timer_lms = toc;
fprintf('LMS known: %1.3f, %1.3f, %1.3f\n', c_lms(1), c_lms(2), c_lms(3));

%RLS w/ known delay
[c_rls, c_rls_help] = deal(zeros(M, 1));
c_rls_ev = zeros(M, length(sr));
rho = 1;
P = rho*eye(M); %eye = identity matrix
lambda = 1;
e_rls = zeros(length(x), 1);
tic;
for i = 1:length(sr)
    for j = 1:M
        if i-delay(j) > 0
            c_rls_help(j) = x(floor(i-delay(j)));
        end
    end
    cPc = c_rls_help'*P*c_rls_help;
    K_tilde = P*c_rls_help/(lambda+cPc);
    P = 1/lambda*(P-P*(c_rls_help*c_rls_help')*P/(lambda+cPc));
    e_rls(i) = sr(i)-c_rls'*c_rls_help;
    c_rls = c_rls+K_tilde*e_rls(i);
    c_rls_ev(:,i) = c_rls;    
end
timer_rls = toc;
fprintf('RLS known: %1.3f, %1.3f, %1.3f\n', c_rls(1), c_rls(2), c_rls(3));

%LMS w/ unknown delay
N_m = 5;
delay_u = fs_d*(1:1:N_m);
M_u = length(delay_u);
[c_lmsu, c_lmsu_help] = deal(zeros(M_u, 1));
[y_lmsu, e_lmsu] = deal(zeros(length(x), 1));
tic;
for i = 1:length(x)
    for j = 1:M_u
        if i-delay_u(j) > 0
            c_lmsu_help(j) = x(i-delay_u(j));
        end
    end
    y_lmsu(i) = c_lmsu'*c_lmsu_help;
    e_lmsu(i) = sr(i)-y_lmsu(i);
    c_lmsu = c_lmsu+2*my*c_lmsu_help*e_lmsu(i);
end
timer_lmsu = toc;
fprintf('LMS unknown:\n');
disp(c_lmsu);

%RLS w/ unknown delay
[c_rlsu, c_rlsu_help] = deal(zeros(M_u, 1));
P_u = rho*eye(M_u);
e_rlsu = zeros(length(x), 1);
tic;
for i = 1:length(x)
    for j = 1:M_u
        if i-delay_u(j) > 0
            c_rlsu_help(j) = x(i-delay_u(j));
        end
    end
    cPcu = c_rlsu_help'*P_u*c_rlsu_help;
    K_tilde = P_u*c_rlsu_help/(lambda+cPcu);
    P_u = 1/lambda*(P_u-P_u*(c_rlsu_help*c_rlsu_help')*P_u/(lambda+cPcu));
    e_rlsu(i) = sr(i)-c_rlsu'*c_rlsu_help;
    c_rlsu = c_rlsu+K_tilde*e_rlsu(i);   
end
timer_rlsu = toc;
fprintf('RLS unknown:\n');
disp(c_rlsu);

%remove echo from signal: LMS w/ known and unknown delay
y_lms_rev = zeros(M_u, length(x));
y_lmsu_rev = zeros(M_u, length(x));
for i = 1:length(sr)
    for j = 1:M
        if i-delay(j) > 0
            y_lms_rev(j, i) = c_lms(j)*x(i-delay(j));
        end
    end
    for j = 1:M_u
        if i-delay_u(j) > 0
            y_lms_rev(j, i) = c_lmsu(j)*x(i-delay_u(j));
        end
    end
end
y_lms_rev_sum = sum(y_lms_rev)'; %sum echo paths
y_lmsu_rev_sum = sum(y_lmsu_rev)';
s_lms_rev = sr-y_lms_rev_sum; %remove echo paths
s_lmsu_rev = sr-y_lmsu_rev_sum;
error_lms = (u-s_lms_rev).^2; %squared error
error_lmsu = (u-s_lmsu_rev).^2;

%remove echo from signal: RLS w/ known and unknown delay
y_rls_rev = zeros(M_u, length(x));
y_rlsu_rev = zeros(M_u, length(x));
for i = 1:length(x)
    for j = 1:M
        if i-delay(j) > 0
            y_rls_rev(j, i) = c_rls(j)*x(i-delay(j));
        end
    end
    for j = 1:M_u
        if i-delay_u(j) > 0
            y_rlsu_rev(j, i) = c_rlsu(j)'*x(i-delay_u(j));
        end
    end
end
y_rls_rev_sum = sum(y_rls_rev)';
y_rlsu_rev_sum = sum(y_rlsu_rev)';
s_rls_rev = sr-y_rls_rev_sum;
s_rlsu_rev = sr-y_rlsu_rev_sum;
error_rls = (u-s_rls_rev).^2;
error_rlsu = (u-s_rlsu_rev).^2;

fprintf('Average sqared error lms: %1.10f\n', mean(error_lms));
fprintf('Average sqared error lmsu: %1.10f\n', mean(error_lmsu));
fprintf('Average sqared error rls: %1.10f\n', mean(error_rls));
fprintf('Average sqared error rlsu: %1.10f\n', mean(error_rlsu));


figure;
subplot(2, 1, 1)
plot(error_lms);
title('Squared error LMS with known delay');
xlabel('Time [s]');
ylabel('Squared error');
subplot(2, 1, 2)
plot(error_lmsu);
title('Squared error LMS with unknown delay');
xlabel('Time [s]');
ylabel('Squared error');

figure;
subplot(2, 1, 1)
plot(error_rls);
title('Squared error RLS with known delay');
xlabel('Time [s]');
ylabel('Squared error');
subplot(2, 1, 2)
plot(error_rlsu);
title('Squared error RLS with unknown delay');
xlabel('Time [s]');
ylabel('Squared error');

%plot recovered signals
figure;
subplot(2, 1, 1);
plot(s_lms_rev);
title('Signal recovered using LMS with known delay');
xlabel('Time [s]');
ylabel('Amplitude');
sound(s_lms_rev, fs_d)
subplot(2, 1, 2);
plot(s_lmsu_rev);
title('Signal recovered using LMS with unknown delay');
xlabel('Time [s]');
ylabel('Amplitude');

figure;
subplot(2, 1, 1);
plot(s_rls_rev);
title('Signal recovered using RLS with known delay');
xlabel('Time [s]');
ylabel('Amplitude');
subplot(2, 1, 2);
plot(s_rlsu_rev);
title('Signal recovered using RLS with unknown delay');
xlabel('Time [s]');
ylabel('Amplitude');

%plot squared error between y and y_hat
y_sqrd_lms = (y-y_lms_rev_sum).^2;
y_sqrd_lmsu = (y-y_lmsu_rev_sum).^2;
y_sqrd_rls = (y-y_rls_rev_sum).^2;
y_sqrd_rlsu = (y-y_rlsu_rev_sum).^2;

figure;
title('Squared errors between echo and its estimate LMS');
subplot(2, 1, 1);
plot(y_sqrd_lms);
title('Squared errors between echo and its estimate LMS with known delays');
xlabel('Time [s]');
ylabel('Squared error');
subplot(2, 1, 2)
plot(y_sqrd_lmsu);
title('Squared errors between echo and its estimate LMS with unknown delay');
xlabel('Time [s]');
ylabel('Squared error');

figure;
title('Squared errors between echo and its estimate RLS');
subplot(2, 1, 1);
plot(y_sqrd_rls);
title('Squared errors between echo and its estimate RLS with known delays');
xlabel('Time [s]');
ylabel('Squared error');
subplot(2, 1, 2);
plot(y_sqrd_rlsu);
title('Squared errors between echo and its estimate RLS with unknown delay');
xlabel('Time [s]');
ylabel('Squared error');

%time-averaged squared error y-y_hat
y_avg_lms = mean(y_sqrd_lms);
fprintf('Time-average sqared error y-y_hat_lms: %1.10f\n', y_avg_lms);
y_avg_lmsu = mean(y_sqrd_lmsu);
fprintf('Time-average sqared error y-y_hat_lmsu: %1.10f\n', y_avg_lmsu);
y_avg_rls = mean(y_sqrd_rls);
fprintf('Time-average sqared error y-y_hat_rls: %1.10f\n', y_avg_rls);
y_avg_rlsu = mean(y_sqrd_rlsu);
fprintf('Time-average sqared error y-y_hat_rlsu: %1.10f\n', y_avg_rlsu);

%plot evolution of 3 estimated echo-path coefficients w/ true values
x_ev = 1:length(c_lms_ev);
c_plot = c'*ones(1, length(c_lms_ev));
figure;
subplot(2, 1, 1);
plot(x_ev, c_lms_ev(1,:), x_ev, c_lms_ev(2,:), x_ev, c_lms_ev(3,:), x_ev, c_plot)
title('Evolution of echo paths compared to set values using LMS')
xlabel('No iterations');
ylabel('Value');

subplot(2, 1, 2);
plot(x_ev, c_rls_ev(1,:), x_ev, c_rls_ev(2,:), x_ev, c_rls_ev(3,:), x_ev, c_plot)
title('Evolution of echo paths compared to set values using RLS')
xlabel('No iterations');
ylabel('Value');

%averaged square error parameter estimates for RLS and LMS
av_lms = sum(mean((c-c_lms_ev').^2));
fprintf('Average sqared error parameter estimates LMS: %1.10f\n', av_lms);
av_rls = sum(mean((c-c_rls_ev').^2));
fprintf('Average sqared error parameter estimates RLS: %1.10f\n', av_rls);

%CPU time different algorithms
fprintf('CPU time LMS known delays: %1.10f\n', timer_lms);
fprintf('CPU time RLS known delays: %1.10f\n', timer_rls);
fprintf('CPU time LMS unknown delays: %1.10f\n', timer_lmsu);
fprintf('CPU time LMS unknown delays: %1.10f\n', timer_rlsu);