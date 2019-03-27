%%Copyright Ilarion Kvitnevskiy 2018

function [ ] = visualize( params, t, q1, q2, file)

m1 = params(1);
m2 = params(2);
l1 = 1000*params(3);
l2 = 1000*params(4);
c1 = params(5);
c2 = params(6);

s = 0.8;
figure('units','normalized','outerposition',[(1-s)/2 (1-s)/2 s s]);

framerate = 30;
frameperiod = 1/framerate;

curtime = -frameperiod;
i = 1;

x1 = l1*cos(q1);
x2 = x1 + l2*cos(q1+q2);
y1 = l1*sin(q1);
y2 = y1 + l2*sin(q1+q2);

display(max(x2));
display(max(y2));

while i < size(t,1)
    curtime = min(curtime + frameperiod, t(end));
    while curtime > t(i)
        i = i + 1;
        continue
    end
    subplot(2,3,[1,4]);
    
    plot(0,0,'MarkerSize',30,'Marker','.', 'Color', 'cyan');
    title('Weeee');
    a1 = (l1+l2)*1.02;
    axis equal;
    axis([0, a1, 0, a1]);
  
    hold on;
    line(1000*[0.22, 0.22], 1000*[0, 0.22], 'LineWidth', 2, 'Color', 'red')
    line(1000*[0, 0.22], 1000*[0.22, 0.22], 'LineWidth', 2, 'Color', 'red')
    plot(1000*0.1,1000*0.2,'MarkerSize',30,'Marker','.', 'Color', 'green');
    plot(1000*0.2,1000*0.2,'MarkerSize',30,'Marker','.', 'Color', 'green');
    plot(1000*0.2,1000*0.1,'MarkerSize',30,'Marker','.', 'Color', 'green');
    plot(1000*0.1,1000*0.1,'MarkerSize',30,'Marker','.', 'Color', 'green');
    
    plot(x2(1:i), y2(1:i), 'MarkerSize', 10, 'Marker', '.', 'Color', 'red')
    
    line([0, x1(i)], [0 y1(i)], 'LineWidth',0.5*m1+0.1);
    plot(x1(i),y1(i),'MarkerSize',30,'Marker','.', 'Color', 'blue');
    
    line([x1(i), x2(i)], [y1(i), y2(i)], 'LineWidth',0.5*m2+0.1);

    hold off;
    
    subplot(2,3, 2);
    plot(t, q1);
    hold on;
    plot(t(i), q1(i), 'MarkerSize', 20, 'Marker', '.', 'Color', 'red');
    hold off;
    title('State uno');

    subplot(2,3, 5);
    plot(t, q2);
    hold on;
    plot(t(i), q2(i), 'MarkerSize', 20, 'Marker', '.', 'Color', 'red');
    hold off;
    title('State dos');

    subplot(2,3, 3);
    plot(t, x2);
    hold on;
    plot(t(i), x2(i), 'MarkerSize', 20, 'Marker', '.', 'Color', 'red');
    hold off;
    title('Pos X')
    
    subplot(2,3, 6);
    plot(t, y2);
    hold on;
    plot(t(i), y2(i), 'MarkerSize', 20, 'Marker', '.', 'Color', 'red');
    hold off;
    title('Pos Y')    
    
    drawnow;
    if ~strcmp(file, '')
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        dt = 0.05;
        if i == 1
            imwrite(imind,cm,file,'gif', 'Loopcount',inf, 'DelayTime', dt);
        else
            imwrite(imind,cm,file,'gif','WriteMode','append', 'DelayTime', dt);
        end        
    end
    
end
end