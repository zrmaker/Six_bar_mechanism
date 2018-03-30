% Author: Renyuan Zhang
% dr.leo.zhang@outlook.com

clc; clear all; close all;
options.Interpreter = 'Latex';

%% Setup

A=7.2; B=32.4; C=17; D=40.5; E=2; F=21.6; G=36.54; H=15.0;

t = 0:.05:10;   % Time series
w = -2*pi/10;   % Angular velocities
theta = w*t;    % Input angle

str1 = 'P2';
str2 = 'P3';
str3 = 'P5';
str4 = 'P6';
str5 = 'P7';

animation = subplot(2,2,1);
axis(animation,'equal');
set(gca,'xlim',[-30 60],'ylim',[-10 50]);
title(animation,'Animation and Trace of P6&P7');
    
dis = subplot(2,2,2);
set(dis,'xlim',[0 max(abs(theta))/pi*180],'ylim',[0 3]);
xlabel(dis, 'Input Angle \theta (\circ)');
ylabel(dis, 'Amplitude (cm)');
title(dis,'Displacement');
grid on;

vel = subplot(2,2,3);
set(vel,'xlim',[0 max(abs(theta))/pi*180],'ylim',[-1 10]);
xlabel(vel, 'Input Angle \theta (\circ)');
ylabel(vel, 'Amplitude (cm/s)');
title(vel,'Speed');
grid on;

acc = subplot(2,2,4);
set(acc,'xlim',[0 max(abs(theta))/pi*180],'ylim',[-10 10]);
xlabel(acc, 'Input Angle \theta (\circ)');
ylabel(acc, 'Amplitude (cm/s^{2})');
title(acc,'Accelaration');
grid on;

%% Points

P1 = [0;0];
P4 = [D;0];
P2 = [A*cos(theta); A*sin(theta)]; 

l24 = sqrt(A^2 + D^2 - 2*A*D*cos(theta));
alpha = asin(A*sin(theta)./l24);
beta = acos((l24.^2+C^2-B^2)./(2*l24*C));
theta2=pi-alpha-beta;
P3 = [D + C*cos(theta2); C*sin(theta2)];

P7 = [D + (C+H)*cos(theta2); (C+H)*sin(theta2)];
P5(1,:)=P2(1,:)-(P3(1,:)-P2(1,:))/B*E;
P5(2,:)=P2(2,:)-(P3(2,:)-P2(2,:))/B*E;

l57=sqrt((P7(2,:)-P5(2,:)).^2+(P7(1,:)-P5(1,:)).^2);

a1=acos((l57.^2+F^2-G^2)./(2*l57*F));
a2=acos((l57.^2+(B+E)^2-H^2)./(2*l57*(B+E)));
a3=atan((P3(2,:)-P5(2,:))./(P3(1,:)-P5(1,:)));
a3(a3<0)=a3(a3<0)+pi;
theta3=a1+a2+a3;
P6(1,:)=P5(1,:)+F*cos(theta3);
P6(2,:)=P5(2,:)+F*sin(theta3);

%% Calculations

P2_d=sqrt((P2(1,:)-P2(1,1)).^2+(P2(2,:)-P2(2,1)).^2);
P2_vx = diff(P2(1,:))./diff(t);
P2_vy = diff(P2(2,:))./diff(t);
P2_v = sqrt(P2_vx.^2 + P2_vy.^2);
P2_a = diff(P2_v(1,:))./diff(t(1:end-1));

P3_d=sqrt((P3(1,:)-P3(1,1)).^2+(P3(2,:)-P3(2,1)).^2);
P3_vx = diff(P3(1,:))./diff(t);
P3_vy = diff(P3(2,:))./diff(t);
P3_v = sqrt(P3_vx.^2 + P3_vy.^2);
P3_a = diff(P3_v(1,:))./diff(t(1:end-1));

P5_d=sqrt((P5(1,:)-P5(1,1)).^2+(P5(2,:)-P5(2,1)).^2);
P5_vx = diff(P5(1,:))./diff(t);
P5_vy = diff(P5(2,:))./diff(t);
P5_v = sqrt(P5_vx.^2 + P5_vy.^2);
P5_a = diff(P5_v(1,:))./diff(t(1:end-1));

P6_d=sqrt((P6(1,:)-P6(1,1)).^2+(P6(2,:)-P6(2,1)).^2);
P6_vx = diff(P6(1,:))./diff(t);
P6_vy = diff(P6(2,:))./diff(t);
P6_v = sqrt(P6_vx.^2 + P6_vy.^2);
P6_a = diff(P6_v(1,:))./diff(t(1:end-1));

P7_d=sqrt((P7(1,:)-P7(1,1)).^2+(P7(2,:)-P7(2,1)).^2);
P7_vx = diff(P7(1,:))./diff(t);
P7_vy = diff(P7(2,:))./diff(t);
P7_v = sqrt(P7_vx.^2 + P7_vy.^2);
P7_a = diff(P7_v(1,:))./diff(t(1:end-1));
%%

for i=1:length(theta)
    
    %% Plot animation
    
    animation = subplot(2,2,1);
    P1_circle = viscircles(P1',0.05);
    P2_circle = viscircles(P2(:,i)',0.05);
    P3_circle = viscircles(P3(:,i)',0.05);
    P4_circle = viscircles(P4',0.05); 
    P5_circle = viscircles(P5(:,i)',0.05);
    P6_circle = viscircles(P6(:,i)',0.05);
    P7_circle = viscircles(P7(:,i)',0.05);
   
    A_bar = line([P1(1) P2(1,i)],[P1(2) P2(2,i)]);
    B_bar = line([P2(1,i) P3(1,i)],[P2(2,i) P3(2,i)]);
    C_bar = line([P3(1,i) P4(1)],[P3(2,i) P4(2)]);
    D_bar = line([P1(1) P4(1)],[P1(2) P4(2)]);
    E_bar = line([P2(1,i) P5(1,i)],[P2(2,i) P5(2,i)]);
    F_bar = line([P5(1,i) P6(1,i)],[P5(2,i) P6(2,i)]);
    G_bar = line([P6(1,i) P7(1,i)],[P6(2,i) P7(2,i)]);
    H_bar = line([P4(1) P7(1,i)],[P4(2) P7(2,i)]);
    
    str6 = ['Input \theta=',num2str(abs(theta(i))/pi*180),'\circ'];
    P2_text = text(P2(1,i),P2(2,i)+4,str1);
    P3_text = text(P3(1,i),P3(2,i)+4,str2);
    P5_text = text(P5(1,i),P5(2,i)-4,str3);
    P6_text = text(P6(1,i),P6(2,i)+4,str4);
    P7_text = text(P7(1,i),P7(2,i)+4,str5);
    angle_text=text(0,40,str6);
    
    pause(0.0001);
    
	%% Plot displacement
    
    dis = subplot(2,2,2);
    hold on;
    plot(dis,abs(theta(1:i))/pi*180,P2_d(1:i)*0.1,'r');
    plot(dis,abs(theta(1:i))/pi*180,P3_d(1:i)*0.1,'b');
    plot(dis,abs(theta(1:i))/pi*180,P5_d(1:i)*0.1,'g');
    plot(dis,abs(theta(1:i))/pi*180,P6_d(1:i)*0.1,'m');
    plot(dis,abs(theta(1:i))/pi*180,P7_d(1:i)*0.1,'k');
    hold off;
    %%
    if i<length(theta)
        
        %% Delete
        
        delete(P1_circle);
        delete(P2_circle);
        delete(P3_circle);
        delete(P4_circle);
        delete(P5_circle);
        %delete(P6_circle);
        %delete(P7_circle);
        delete(A_bar);
        delete(B_bar);
        delete(C_bar);
        delete(D_bar);
        delete(E_bar);
        delete(F_bar);
        delete(G_bar);
        delete(H_bar);
        delete(P2_text);
        delete(P3_text);
        delete(P5_text);
        delete(P6_text);
        delete(P7_text);
        delete(angle_text);
        
        %% Plot velocity
        
        vel = subplot(2,2,3);
        hold on;
        plot(vel,abs(theta(1:i))/pi*180,P2_v(1:i),'r');
        plot(vel,abs(theta(1:i))/pi*180,P3_v(1:i),'b');
        plot(vel,abs(theta(1:i))/pi*180,P5_v(1:i),'g');
        plot(vel,abs(theta(1:i))/pi*180,P6_v(1:i),'m');
        plot(vel,abs(theta(1:i))/pi*180,P7_v(1:i),'k');
        hold off;
        %%
        if i<length(theta)-1
            
            %% Plot Acceleration
            
            acc = subplot(2,2,4);
            hold on;
            plot(acc,abs(theta(1:i))/pi*180,P2_a(1:i),'r');
            plot(acc,abs(theta(1:i))/pi*180,P3_a(1:i),'b');
            plot(acc,abs(theta(1:i))/pi*180,P5_a(1:i),'g');
            plot(acc,abs(theta(1:i))/pi*180,P6_a(1:i),'m');
            plot(acc,abs(theta(1:i))/pi*180,P7_a(1:i),'k');
            hold off;
        end
    end
end

legend(dis,'Point 2','Point 3','Point 5','Point 6','Point 7');
legend(vel,'Point 2','Point 3','Point 5','Point 6','Point 7');
legend(acc,'Point 2','Point 3','Point 5','Point 6','Point 7');