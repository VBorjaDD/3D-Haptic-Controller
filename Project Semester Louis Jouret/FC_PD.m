clear all;
close all;

time = 100;

V= VideoWriter('PDMovie','MPEG-4');
fps = 10;
V.FrameRate = fps;


model = stlread('Cellulo concept 1 v47.stl');
Points = model.Points;
Points = Points/1000;
CList = model.ConnectivityList;

%x = sym('x',[12 1]);
%u = sym('u',[4 1]);
xz = sym('xz',[15 1]);
% xz = sym('xz',[12 1]);
para = sym('para',[16,1]);
% para = sym('para',[4,1]);


I1 = 4.5;
I2 = 4.3;
I3 = 6.0;

kf = 1;
km = 1;
a = 0.2;
m = 1;
g = 9.81;

%% small angle approx
            
fun = @(t,xz, para) [       xz(4);
                            xz(5);
                            xz(6);
                            (1/I1)*(-xz(6)*xz(5)*(I2-I3) + kf*a*(para(2) - para(4)));
                            (1/I2)*(-xz(4)*xz(6)*(I3-I1) + kf*a*(para(3) - para(1)));
                            (1/I3)*(-xz(4)*xz(5)*(I1-I2) + km*(para(1) - para(2) + para(3) - para(4)));
                            xz(10);
                            xz(11);
                            xz(12);
                            (kf/m)*(para(1) + para(2) + para(3) + para(4))*xz(2);
                            (kf/m)*(para(1) + para(2) + para(3) + para(4))*xz(1);
                            -g + (kf/m)*(para(1) + para(2) + para(3) + para(4))
                            xz(1) - para(5);
                            xz(2) - para(6);
                            xz(3) - para(7);
                            xz(4) - para(8);
                            xz(5) - para(9);
                            xz(6) - para(10);
                            xz(7) - para(11);
                            xz(8) - para(12);
                            xz(9) - para(13);
                            xz(10) - para(14);
                            xz(11) - para(15);
                            xz(12) - para(16);
                            ];

h = 1/fps;

state = [0.1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
input = [0; 0; 0; 0];

a = 60.0;
c = 30.0; %pb
b = 2.5;
d = 5.0;  %pb


Kpd = [ 0   a   -a   0    c     -c   b   0   -b  d   0  -d;
        -a  0   a    -c   0     c    0   -b  -b  0  -d  -d;
        0   -a  -a   0    -c    -c   -b  0   -b  -d  0  -d;
        a   0   a    c    0     c    0   b   -b  0   d  -d];

Ki = Kpd
        
        
tspan = [0 h];

reference =     [0; 0; 0; 0; 0; 0; 0; 1; 0; 0; 0; 0];
NominalInput =  [m*g/4; m*g/4; m*g/4; m*g/4];
intg = zeros(12,1);
XZ = [state;intg];
% XZ = state;
U = NominalInput;
tspan = [0 h];
Ref = reference;
para = zeros(16,1);



for i = 1:fps*time
    input0 = Kpd*(state-reference) + 0.0*Ki*intg;
    input = input0 + NominalInput;
    input(input<0) = 0;
    input(input>20) = 20;
    para = [input; reference];
%     xz = [x ;z];
    [t,sol] = ode45(@(t, xz) fun(t,xz,para), tspan, [state; intg]);
%     [t,sol] = ode45(@(t, xz) fun(t,xz,para), tspan, state);
    intg = sol(end,13:24).';
    oldref = reference;
    reference(7:9) = rotz(1)*reference(7:9);
    reference(9) = reference(9) + 0.01;
    reference(10:12) = (reference(7:9) - oldref(7:9))/h;
    Noise = 0.005*(randn(12,1)-0.5);
    state = sol(end,1:12).' + Noise;
    XZ = [XZ [state;intg]];
%     XZ = [XZ state];
    U = [U input0];
    Ref = [Ref reference];
    i
end


plot1 = figure;
plot(linspace(0,fps*time,fps*time + 1), XZ(1,:))
hold on;
plot(linspace(0,fps*time,fps*time + 1), XZ(2,:))
hold on;
plot(linspace(0,fps*time,fps*time + 1), XZ(3,:))
hold off;
xlabel('time [s]') 
ylabel('inclination [rad]') 
legend({'roll','pitch','yaw'})


plot2 = figure;
plot(linspace(0,fps*time,fps*time + 1), XZ(7,:))
hold on;
plot(linspace(0,fps*time,fps*time + 1), XZ(8,:))
hold on;
plot(linspace(0,fps*time,fps*time + 1), XZ(9,:))
hold on;
plot(linspace(0,fps*time,fps*time + 1), Ref(7,:),'--')
hold on;
plot(linspace(0,fps*time,fps*time + 1), Ref(8,:),'--')
hold on;
plot(linspace(0,fps*time,fps*time + 1), Ref(9,:),'--')
hold off;
xlabel('time [s]') 
ylabel('coordinates [m]') 
legend({'x','y','z'})


% movie = figure;
% open(V);
% G = [];
% for i = 1:fps*time
%     Rx = rotx(XZ(1,i)*180/pi);
%     Ry = roty(XZ(2,i)*180/pi);
%     Rz = rotz(XZ(3,i)*180/pi);
%     rotatedPoints = Points*Rx*Ry*Rz;
%     G = [G  [XZ(7,i); XZ(8,i); XZ(9,i)]];
%     Translation = [XZ(7,i)*ones(length(Points),1), XZ(8,i)*ones(length(Points),1), XZ(9,i)*ones(length(Points),1)];
%     translatedPoints = rotatedPoints + Translation;
%     newPos = triangulation(CList,translatedPoints);
%     trisurf(newPos);
%     hold on;
%     fnplt(cscvn(G(:,1:end)),'r',2)
%     hold on;
%     scatter3(Ref(7,i),Ref(8,i),Ref(9,i),'*','g');
%     hold on;
%     fnplt(cscvn(Ref(7:9,1:i)),'g',2)
%     hold off;
%     axis(2*[-1,1,-1,1,-1,10])
%     frame = getframe;
%     writeVideo(V,frame);
%     i
% end
% close(V);


 
