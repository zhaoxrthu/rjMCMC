function rjMCMCtest(x1,y1,x2,y2)
    tic
    ytest=rjMCMCSA(x1(1:800,:),y1(1:800,:),x1(801:1000,:));
    deltay=ytest-y1(801:1000,:);
    error=1/300*sum(sum(deltay.*deltay),2)%/sum(sum(y1(701:1000,:).*y1(701:1000,:),2))
    figure;
    hold on;
    title('data1');
    plot(y1(801:1000,1),y1(801:1000,2),'b.');
    plot(ytest(1:200,1),ytest(1:200,2),'r*');
    xlabel('y1(:,1)');ylabel('y1(:,2)');
    legend('Real Position','Predicted Position');
    xlim([-4,4]);ylim([-4,4]);
    toc

    tic
    ytest=rjMCMCSA(x2(1:800,:),y2(1:800,:),x2(801:1000,:));
    deltay=ytest-y2(801:1000,:);
    error=0;
    for i=1:200
        error=error+norm(deltay(i,:))^2;
    end
    error=error/200
    % error=1/300*sum(sum(deltay.*deltay),2)%/sum(sum(y2(701:1000,:).*y2(701:1000,:),2))
    figure;
    title('data2');
    hold on;
    plot(y2(801:1000,1),y2(801:1000,2),'b.');
    plot(ytest(1:200,1),ytest(1:200,2),'r*');
    xlabel('y2(:,1)');ylabel('y2(:,2)');
    legend('Real Position','Predicted Position');
    toc
end