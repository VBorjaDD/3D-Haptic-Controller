clear all;
close all;


addpath('/Users/louisjouret/Documents/GitHub/ADCS_LouisJouret/casadi-matlabR2015a-v3.3.0')
import casadi.*

model = stlread('Cellulo concept 1 v47.stl');
Points = model.Points;
Points = Points/1000;
CList = model.ConnectivityList;


%% physical characteristics of the drone and constants
I1 = 4.5;
I2 = 4.3;
I3 = 6.0;

kf = 1;
km = 1;
a = 0.2;
m = 1;
g = 9.81;

%% diff equations of the system

x = MX.sym('x',12);
u = MX.sym('u',4);

 ode =   [      x(4);
                x(5);
                x(6);
                (1/I1)*(-x(6)*x(5)*(I2-I3) + kf*a*(u(2) - u(4)));
                (1/I2)*(-x(4)*x(6)*(I3-I1) + kf*a*(u(3) - u(1)));
                (1/I3)*(-x(4)*x(5)*(I1-I2) + km*(u(1) - u(2) + u(3) - u(4)));
                x(10);
                x(11);
                x(12);
                (kf/m)*(u(1) + u(2) + u(3) + u(4))*x(2);
                (kf/m)*(u(1) + u(2) + u(3) + u(4))*x(1);
                -g + (kf/m)*(u(1) + u(2) + u(3) + u(4))];

f = Function('f',{x,u},{ode},{'x','u'},{'ode'});


T   = 1.0; % Time horizon
N   = 15; % Number of control intervals
time= 8.0; %time simulated

% Integrator to discretize the system
intg_options = struct;
intg_options.tf = T/N; %integrate over T/N time
intg_options.number_of_finite_elements = 4; % Rugge-Kutta finite elements

% DAE problem structure
dae    = struct;
dae.x  = x;         % states
dae.p  = u;         % parameters (here the inputs)7
dae.ode= f(x,u);    % Expression for the right-hand side

intg   = integrator('intg','rk',dae,intg_options);
res    = intg('x0',x,'p',u); % Evaluate with symbols
x_next = res.xf;
% map the initial x&u to the next state
F      = Function('F',{x,u},{x_next},{'x','u'},{'x_next'});

V= VideoWriter('Drone2Movie','MPEG-4');
fps = N/T;
V.FrameRate = fps;

%% optimization loop

opti = casadi.Opti();
x    = opti.variable(12,N+1); % Decision variables for state trajectory
u    = opti.variable(4,N);
p    = opti.parameter(12,1);  % Parameter (not optimized over)
ref  = opti.parameter(12,1);

Q = diag([1 1 1 1 1 1 10 10 10 1 1 1]); 

J = 0;

for k=1:N
    opti.subject_to(x(:,k+1)==F(x(:,k),u(:,k)));
    J = J + ((x(:,k)-ref).'*Q*(x(:,k)-ref));
end

opti.minimize(J); %minimize the cost function
opti.subject_to(0 <= u <= 20); %constraints on the input
opti.subject_to(x(:,1) == p); %first state is the previous state
opti.solver('ipopt');

%%
reference =     [0; 0; 0; 0; 0; 0; 0; 1; 0; 0; 0; 0];
x_new =         [0; 0; 0; 0; 0; 0; 0; 1; 0; 0; 0; 0];
Noise =         zeros(12,1);

X = x_new;
U = [];
Sol = [];
Ref = [];


opti.set_value(p,x_new);
opti.set_value(ref, reference);

for i=1:N*time/T
    sol     = opti.solve();
    input   = sol.value(u(:,1));
    Sol = [Sol sol.value(J)];
    x_new   =  F(x_new,input) + Noise;
    x_new   = full(x_new);
    X = [X x_new];
    U = [U input];
    Noise = 0.001*(randn(12,1)-0.5);
    reference(7:9) = rotz(3)*reference(7:9);
    reference(9) = reference(9) + 0.01;
    Ref = [Ref reference];
    opti.set_value(ref, reference);
    opti.set_value(p,x_new); 
    i
end


plot1 = figure;
plot(linspace(0,N*time/T,N*time/T + 1), X(1,:))
hold on;
plot(linspace(0,N*time/T,N*time/T + 1), X(2,:))
hold on;
plot(linspace(0,N*time/T,N*time/T + 1), X(3,:))
xlabel('time [s]') 
ylabel('inclination [rad]') 
legend({'roll','pitch','yaw'})


plot2 = figure;
plot(linspace(0,N*time/T,N*time/T + 1), X(7,:))
hold on;
plot(linspace(0,N*time/T,N*time/T + 1), X(8,:))
hold on;
plot(linspace(0,N*time/T,N*time/T + 1), X(9,:))
xlabel('time [s]') 
ylabel('coordinates [m]') 
legend({'x','y','z'})

plot3 = figure;
plot(Sol);
xlabel('time [s]') 
ylabel('cost function value J') 

movie = figure;
open(V);
G = [];
for i = 1:N*time/T
    Rx = rotx(X(1,i)*180/pi);
    Ry = roty(X(2,i)*180/pi);
    Rz = rotz(X(3,i)*180/pi);
    rotatedPoints = Points*Rx*Ry*Rz;
    G = [G  [X(7,i); X(8,i); X(9,i)]];
    Translation = [X(7,i)*ones(length(Points),1), X(8,i)*ones(length(Points),1), X(9,i)*ones(length(Points),1)];
    translatedPoints = rotatedPoints + Translation;
    newPos = triangulation(CList,translatedPoints);
    trisurf(newPos);
    hold on;
    fnplt(cscvn(G(:,1:end)),'r',2)
    hold on;
    scatter3(Ref(7,i),Ref(8,i),Ref(9,i),'*','g');
    hold on;
    fnplt(cscvn(Ref(7:9,1:i)),'g',2)
    hold off;
    axis(2*[-1,1,-1,1,-1,10])
    frame = getframe;
    writeVideo(V,frame);
    i
end
close(V);


