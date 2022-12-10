% Golnoosh Aghamohammadi


close all;
clear;
clc;


a = 0.1; 
b = 0.2;
omega = 1; 
error = 1e-9;
t = linspace(0, 1, 101);
X = [pi; 0]; 



for i=1:length(t)
    phi = pi/6 + omega * t(i);
    F = @(X) [a * cos(phi) + b * cos(X(1)) - X(2)
         a * sin(phi) - b * sin(X(1))]; % Constraint Matrix
    J = @(X) [-b * sin(X(1)), -1
         -b * cos(X(1)), 0]; %Jacobian Matrix

    [X, iteration_counter] = NewtonRaphson_method(F, J, X, error);

     dFt = [-a * omega * sin(omega * t(i))
             a * omega * cos(omega * t(i))];

    dFq = [-b * sin(X(1)), -1
           -b * cos(X(1)), 0];

    d_X = - dFt \ dFq;
    
    theta(i,1) = X(1);
    d(i,1) = X(2);
    d_theta(i,1) = d_X(1);
    d_d(i,1) = d_X(2);

    
end

figure(1)

subplot(2,2,1)
plot(t, theta, 'LineWidth',2);
ylabel('\Theta');
xlabel('Time(s)');

subplot(2,2,2);
plot(t, d,'LineWidth',2);
ylabel('X');
xlabel('Time(s)');

subplot(2,2,3)
plot(t, d_theta,'LineWidth',2);
ylabel('d_\Theta');
xlabel('Time(s)');

subplot(2,2,4);
plot(t, d_d,'LineWidth',2);
ylabel('d_X');
xlabel('Time(s)');


% Newton-Raphson method
function [x, iteration_counter] = NewtonRaphson_method(F, J, x, tol)
    F_value = F(x);
    F_norm = norm(F_value);
    iteration_counter = 0;

    while F_norm > tol && iteration_counter < 100
        delta = J(x)\-F_value;
        x = x + delta;
        F_value = F(x);
        F_norm = norm(F_value);
        iteration_counter = iteration_counter + 1;
    end

    if F_norm > tol
        iteration_counter = -1;
    end
end


