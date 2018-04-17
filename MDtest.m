function MDtest()
    day=1:40907;
    daylength=length(day);
    %the probability of mine disaster happened everyday 
    %p(x=1)=lamda.*e^(-lamda)
    lamda=rand(1,6);
    pt=lamda.*exp(-lamda);
    p=ones(40907,1);
    interlength=floor(length(day)/length(pt));
    for i=1:length(pt)
        p((i-1)*interlength+1:i*interlength)=pt(i);
    end
    p(i*interlength+1:end)=0;
    %which days mine disaster happened
    happendays=[];
    for i=1:daylength
        if rand()<p(i)
            happendays=[happendays,i];
        end
    end
    
    lamdare=MineDisaster(daylength,happendays);
    figure;
    subplot(2,1,1);plot(day,p);
    xlabel('t(days)');ylabel('\lambda');
    legend('\lambda');
    xlim([1,daylength]);
    title('Real Parameter');
    subplot(2,1,2);plot(day,lamdare);
    xlabel('t(days)');ylabel('\lambda');
    xlim([1,daylength]);
    legend('\lambda');
    title('Predicted Parameter');

end