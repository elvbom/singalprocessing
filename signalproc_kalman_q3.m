close all
clear

%Q3

sigman = 0.25;
n = eye(2)*sigman;


%Method 1 w/ noise
%params
T = 0.1;
samples = 100;

q = 0.2;
sigma1 = 0.2;
sigma2 = 0.3;

A = [1 T; 0 1];
Q = [T^4/4 T^3/2; T^3/2 T^2]*q;
[C_y1, C_y2] = deal([1 0]);
C = [C_y1; C_y2];
R = [sigma1 0; 0 sigma2]+n;

%apri apost
[x, y, x_pri, x_post, E] = deal(zeros(2, samples));
[p_pri, p_post, K] = deal(zeros([2 2 samples]));
[p_pri(:,:,1), p_pri(:,:,1)] = deal(eye(2));

[trace_pri11, trace_post1] = deal(length(samples));

%initial condition/guesses
x(:,1) = A*[0 0]' + sqrt(Q)*randn(2, 1);
y(:,1) = C*x(:,1) + sqrt(R)*randn(2, 1);

x_pri(:,1) = (mvnrnd([0 0]', Q))'; 
E(:,1) = y(:,1) - C*x_pri(:,1);
p_pri(:,:,1) = Q;

K(:,:,1) = p_pri(:,:,1)*C'*(C*p_pri(:,:,1)*C' + R)^(-1); 
p_post(:,:,1) = p_pri(:,:,1)*(1 - C*K(:,:,1));
x_post(:,1) = x_pri(:,1) + K(:,:,1)*E(:,1);

for i = 2:samples
    %actual
    x(:,i) = A*x(:,i-1) + sqrt(Q)*randn(2, 1);
    y(:,i) = C*x(:,i) + sqrt(R)*randn(2, 1);
    
    %estimates
    %a priori prediction
    x_pri(:,i) = A*x_post(:,i-1);
    p_pri(:,:,i) = A*p_post(:,:,i-1)*A' + Q;
    
    %a posteriori correction
    K(:,:,i) = p_pri(:,:,i)*C'*(C*p_pri(:,:,i)*C' + R)^(-1);
    E(:,i) = y(:,i) - C*x_pri(:,i);
    x_post(:,i) = x_pri(:,i) + K(:,:,i)*E(:,i);
    p_post(:,:,i) = (eye(2) - K(:,:,i)*C)*p_pri(:,:,i);
    
    %traces
    trace_pri11(i) = trace(p_post(:,:,i));
    trace_post1(i) = trace(p_pri(:,:,i));
end

k = 1:samples;

%1.c) plot position
figure;
subplot(2, 1, 1);
plot(k, x(1,1:samples),'k', k , x_pri(1,1:samples), 'r--', k, x_post(1,1:samples),'b--')
legend('True position', 'Apriori estimate', 'Aposteriori estimate')
ylabel('Position')
title('Noise corrupted signal - method 1: position and velocity estimates')

%1.d) plot velocity
subplot(2, 1, 2), 
plot(k, x(2,1:samples), 'k', k, x_pri(2,1:samples), 'r--', k, x_post(2,1:samples), 'b--')
ylabel('Velocity')
legend('True velocity', 'Apriori estimate', 'Aposteriori estimate')
xlabel('Time')

%1.e) error covariance maatrix
figure;
plot(k, trace_pri11, 'r', k, trace_post1, 'b')
title('Noise corrupted signal - method 1: traces error covariance matrix')
xlabel('Time')
legend('Trace error covariance matrix apriori', 'Trace error covariance matrix aposteriori')
axis([2 100 0 0.1]);

%1.f) What happens to x?1k|k?1 as k ? ??
%more more acurate

%1.g) What happens to Pk|k?1 as k ? ?? 
%Write down the limiting value of the covariance matrix Pk|k?1, if it exists.
fprintf('Lim value Pk|k-1: %1.10f\n', trace(p_pri(:,:,samples)));


%Method 2 w/ noise
%params
C = C_y1 + C_y2;
R = sigma1 + sigma2 + sigman;

%apri apost
[E,  y] = deal(zeros(1, samples));
[x, x_pri, x_post, K] = deal(zeros(2, samples));
[p_pri, p_post] = deal(zeros([2 2 samples]));
[p_pri(:,:,1), p_pri(:,:,1)] = deal(eye(2));

[trace_pri12, trace_post1] = deal(length(samples));

%initial condition/guesses
x(:,1) = A*[0 0]' + sqrt(Q)*randn(2, 1);
y(1) = C*x(:,1) + sqrt(R)*randn;

x_pri(:,1) = (mvnrnd([0 0]', Q))'; 
E(1) = y(1) - C*x_pri(:,1);
p_pri(:,:,1) = Q;

K(:,1) = p_pri(:,:,1)*C'/(C*p_pri(:,:,1)*C' + R); 
p_post(:,:,1) = p_pri(:,:,1)*(1 - C*K(:,1));
x_post(:,1) = x_pri(:,1) + K(:,1)*E(1);

for i = 2:samples
    %actual
    x(:,i) = A*x(:,i-1) + sqrt(Q)*randn(2, 1);
    y(i) = C*x(:,i) + sqrt(R)*randn;
    
    %estimates
    %a priori prediction
    x_pri(:,i) = A*x_post(:,i-1);
    p_pri(:,:,i) = A*p_post(:,:,i-1)*A' + Q;
    
    %a posteriori correction
    K(:,i) = p_pri(:,:,i)*C'*(C*p_pri(:,:,i)*C' + R)^(-1);
    E(i) = y(i) - C*x_pri(:,i);
    p_post(:,:,i) = p_pri(:,:,i)*(1 - C*K(:,i));
    x_post(:,i) = x_pri(:,i) + K(:,i)*E(i);
    
    %traces
    trace_pri12(i) = trace(p_post(:,:,i));
    trace_post1(i) = trace(p_pri(:,:,i));
end

k = 1:samples;

%plot position
figure;
subplot(2, 1, 1);
plot(k, x(1,1:samples),'k', k , x_pri(1,1:samples), 'r--', k, x_post(1,1:samples),'b--')
legend('True position', 'Apriori estimate', 'Aposteriori estimate')
ylabel('Position')
title('Noise corrupted signal - method 2: position and velocity estimates')

%plot velocity
subplot(2, 1, 2), 
plot(k, x(2,1:samples), 'k', k, x_pri(2,1:samples), 'r--', k, x_post(2,1:samples), 'b--')
ylabel('Velocity')
legend('True velocity', 'Apriori estimate', 'Aposteriori estimate')
xlabel('Time')

%error covariance matrix
figure;
plot(k, trace_pri12, 'r', k, trace_post1, 'b')
title('Noise corrupted signal - method 2: traces error covariance matrix')
xlabel('Time')
legend('Trace error covariance matrix apriori', 'Trace error covariance matrix aposteriori')
axis([2 100 0 0.1]);


%1.f) What happens to x?1k|k?1 as k ? ??
%more more acurate

%1.g) What happens to Pk|k?1 as k ? ?? 
%Write down the limiting value of the covariance matrix Pk|k?1, if it exists.
fprintf('Lim value Pk|k-1: %1.10f\n', trace(p_pri(:,:,samples)));

%Q3.a)
figure;
plot(k, trace_pri11,'b', k, trace_pri12, 'r' )
title('Noise corrupted signal - comparison error covariance matrices for method 1 and method 2')
legend('Method 1' , 'Method 2')
xlabel('Time')
axis([2 100 0 0.1]);