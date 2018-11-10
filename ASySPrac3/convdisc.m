function convdisc(x1,x2)
    Tam1 = size(x1(1,:));
    Tam2 = size(x2(1,:));
    T1 = Tam1(2);
    T2 = Tam2(2);
    T3 = T1+T2-1;
    vec1 = zeros(1,T1+2*T2-2);
    vec2 = zeros(1,T1+2*T2-2);
    Tec2 = -(T2-1):(T1+T2-2);
    size(vec2);
    size(Tec2);
    x3=zeros(2,T3);
    x3(1,:)=0:T3-1;
    x3(1,:)=x3(1,:)+x1(1,1)+x2(1,1);
    vec1(1,T2:T2+T1-1)=x1(2,:);
    figure (1)
    hFig = figure(1);
    set(hFig, 'Position', [0 0 1000 1000])
    filename = 'testAnimated.gif';
    ttt=1;
    for m=0:T1+T2-2
        vec2(1,1+m:T2+m)=fliplr(x2(2,:));
        x3(2,m+1)=sum(vec1.*vec2);
        subplot(4,1,4)
        stem(x3(1,:),x3(2,:),'r','LineWidth',3)
        axis([Tec2(1),T1+2*T2-2,  min(x3(2,:))-1 max(x3(2,:))+1])
        ylabel('x_3[n]','FontWeight','bold','FontSize',16);
        xlabel('n','FontWeight','bold','FontSize',16);
        grid on
        subplot(4,1,1)
        stem(x2(1,:),x2(2,:),'k','LineWidth',3)
        ylabel('x_2[m]','FontWeight','bold','FontSize',16);
        xlabel('m','FontWeight','bold','FontSize',16);
        title('CONVOLUCIÛN','FontWeight','bold','FontSize',16);
        grid on
        axis([Tec2(1),T1+2*T2-2, min(x2(2,:))-1 max(x2(2,:))+1])
        subplot(4,1,2)
        stem(x1(1,:),x1(2,:),'b','LineWidth',3)
        ylabel('x_1[m]','FontWeight','bold','FontSize',16);
        xlabel('m','FontWeight','bold','FontSize',16);
        grid on
        axis([Tec2(1),T1+2*T2-2, min(x1(2,:))-1 max(x1(2,:))+1])
        subplot(4,1,3)
        stem(Tec2,vec2(1,:),'g','LineWidth',3)
        axis([Tec2(1),T1+2*T2-2,min(vec2)-1, max(vec2)+1 ])
        ylabel('x_2[n-m]','FontWeight','bold','FontSize',16);
        xlabel('m','FontWeight','bold','FontSize',16);
        grid on
        drawnow;
        pause(0.4)
        vec2 = zeros(1,T1+2*T2-2);
        frame = getframe(hFig); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 
        % Write to the GIF File 
        if ttt == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
        else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
        end
        ttt=ttt+1;
    end
end