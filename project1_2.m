clc,clear all;
%2(1)
%     tic
%     MDtest();
%     toc


%2(2)
    load('data1.mat');
    x1=x;y1=y;xtest1=xtest;
    load('data2.mat');
    x2=x;y2=y;xtest2=xtest;
%     rjMCMCtest(x1,y1,x2,y2);
%     tic
%     v1=rjMCMCSA(x1,y1,xtest1);
%     toc
%     tic
%     v2=rjMCMCSA(x2,y2,xtest2);
%     toc
%     save('RJMCMC2015011046.mat','v1','v2');
%     figure;
%     plot(v1(:,1),v1(:,2),'.')
%     figure;
%     plot(v2(:,1),v2(:,2),'.')

%     [error,likelihood]=Likelihood(x1(1:800,:),y1(1:800,:),x1(801:1000,:),y1(801:1000,:));
%     save('err_data1.mat','error');
%     save('llh_data1.mat','likelihood');
    load('llh_data1.mat');
    llh1=likelihood;
    load('llh_data2.mat');
    llh2=likelihood; 
%     BICtest(llh1,x1,y1,llh2,x2,y2);

    v1=BIC(llh1,x1,y1,xtest1);
    v2=BIC(llh2,x2,y2,xtest2);
    save('BIC2015011046.mat',vi,v2);