function [err,llh]=Likelihood(inputdata,outputdata,inputtest,outputtest)
    [N,d]=size(inputdata);
    [Ntest,~]=size(inputtest);
	[~,c]=size(outputdata);
    Kmax=150;
    Iiter=4000;
    lamda=0.1;
    errortest=zeros(Kmax,1);
    likelihood=zeros(Kmax,1);
%     for k=1:Kmax
    for k=36:36
        tic
        k
        maxinput=max(inputdata);
        mininput=min(inputdata);
        miu=mininput+(maxinput-mininput).*rand(k,d);
        D=GenerateD(inputdata,N,k,miu,lamda);
        P=eye(N)-D*pinv(D'*D)*D';
        Pcoe=k*(c+1)*log(N)/2;
        for i=1:Iiter
           [miure,Dre,Pre,~]=UpdataMove(inputdata,outputdata,mininput,maxinput,k,miu,N,Pcoe,lamda,D,P);
            miu=miure;
            D=Dre;
            P=Pre;
        end
        %finding the likelihood
        for j=1:c
            likelihood(k)=likelihood(k)+N/2*log(outputdata(:,j)'*P*outputdata(:,j));
        end
        
        %finding the error
        Dtemp=GenerateD(inputdata,N,k,miu,lamda);
        alpha=zeros(1+d+k,c);
        for j=1:c
           alpha(:,j)=pinv(Dtemp'*Dtemp)*Dtemp'*outputdata(:,j);
        end
        Dtemp=GenerateD(inputtest,Ntest,k,miu,lamda);
        outputtestre=Dtemp*alpha;
        for i=1:Ntest
            errortest(k)=errortest(k)+norm(outputtestre(i,:)-outputtest(i,:))^2;
        end
        figure;
        plot(outputtestre(:,1),outputtestre(:,2),'b.',outputtest(:,1),outputtest(:,2),'r*');
        title(k);
        toc
    end
    err=errortest./Ntest;
    llh=likelihood;
end

function D=GenerateD(inputdata,N,k,miu,lamda)
    %Generating the D matrix;    
	fi=zeros(N,k);
    for ncount=1:N
        for kcount=1:k
			distan1=inputdata(ncount,:)-miu(kcount,:);
			fi(ncount,kcount)=exp(-lamda*sum(distan1.*distan1,2));
        end
    end
	D=[ones(N,1),inputdata,fi];
end

function [miure,Dre,Pre,flag]=UpdataMove(inputdata,outputdata,mininput,maxinput,k,miu,N,Pcoe,lamda,D,P)
	[~,c]=size(outputdata);
    [~,d]=size(inputdata);
    flag=false;
	miutemp=mininput+(maxinput-mininput).*rand(1,d);
	updatak=max(ceil(k*rand()),1);
    miunext=miu;
	miunext(updatak,:)=miutemp;
    Dnext=D;
    for j=1:N
        distan1=inputdata(j,:)-miutemp;
        Dnext(j,1+d+updatak)=exp(-lamda*sum(distan1.*distan1,2));
    end
    Pnext=eye(N)-Dnext*pinv(Dnext'*Dnext)*Dnext';
% 	r=exp(-Pcoe);
    r=1;
	for j=1:c
		r=r*power(outputdata(:,j)'*P*outputdata(:,j)/(outputdata(:,j)'*Pnext*outputdata(:,j)),N/2);
	end
	A=min(1,r);
	if rand()<A
		miure=miunext;
        Dre=Dnext;
        Pre=Pnext;
		flag=true;
    else
        miure=miu;
        Dre=D;
        Pre=P;
	end
end