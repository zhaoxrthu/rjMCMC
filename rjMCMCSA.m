function [outputtestre]=rjMCMCSA(inputdata,outputdata,inputtest)
	[N,d]=size(inputdata);
    [Ntest,~]=size(inputtest);
	[~,c]=size(outputdata);
    
	k=70;
	bk=0.2;dk=0.2;sk=0.2;mk=0.2;
	T=1;
    %Ccoe=c+1;    
    Ccoe=(c+1)/2*log(N/2);

    Iiter=1000;
	maxinput=max(inputdata);
	mininput=min(inputdata);
%     ct=1/2.*(maxinput+mininput);
%     maxinput=maxinput+ct;
%     mininput=maxinput-ct;
    zeta=1;
	miu=mininput+(maxinput-mininput).*rand(k,d);  
    ktest=zeros(Iiter,1);
%     outputtest=zeros(Ntest,c);
%     error=zeros(Iiter-Ilimit,1);
    P=GenerateP(inputdata,N,k,miu);
    
	for i=1:Iiter
%         T=T*0.995;
        Pcoe=Ccoe*k;
		%MH
		u=rand();
		if u<=bk
			[kre,miure,Pre,flag]=BirthMove(inputdata,outputdata,mininput,maxinput,k,miu,N,Ccoe,P);
        else
            if u<=(bk+dk)
			[kre,miure,Pre,flag]=DeathMove(inputdata,outputdata,k,miu,N,Ccoe,P);
            else
                if u<=(bk+dk+sk)
                [kre,miure,Pre,flag]=SplitMove(inputdata,outputdata,k,miu,N,Ccoe,zeta,P);
                else
                    if u<=(bk+dk+sk+mk)
                    [kre,miure,Pre,flag]=MergeMove(inputdata,outputdata,k,miu,N,Ccoe,zeta,P);
                    else 
                        [kre,miure,Pre,flag]=UpdataMove(inputdata,outputdata,mininput,maxinput,k,miu,N,Pcoe,P);
                    end
                end
            end
        end
        if flag
        	Arjsa=SimAnneal(k,kre,miure,P,Ccoe,T,inputdata,outputdata);
        	uT=rand();
        	if uT<=Arjsa
        		k=kre;
        		miu=miure;
        	end
        end
%         if mod(i,100)==0
%             i
%         end
        ktest(i)=k;
        k=kre;
        miu=miure;
        P=Pre;

    end
    Dtemp=GenerateD(inputdata,N,k,miu);
    alpha=zeros(1+d+k,c);
    for j=1:c
       alpha(:,j)=pinv(Dtemp'*Dtemp)*Dtemp'*outputdata(:,j);
    end
    Dtemp=GenerateD(inputtest,Ntest,k,miu);
    outputtestre=Dtemp*alpha;

    figure;
    plot(ktest);
    xlabel('Iteration Times');ylabel('k');
    legend('k');
end

function D=GenerateD(inputdata,N,k,miu)
    %Generating the D matrix;    
    lamda=0.10;
	fi=zeros(N,k);
	for ncount=1:N
		for kcount=1:k
			distan1=inputdata(ncount,:)-miu(kcount,:);
			fi(ncount,kcount)=exp(-lamda*sum(distan1.*distan1,2));
		end
	end
%     theta=[2,1;1,2];
%     fi=zeros(N,k);
%     for i = 1:N
%         fi(i,:)=1/(2*pi*sqrt(det(theta)))*(sum((exp(-1/2*(inputdata(i,:)-miu)/theta.*(inputdata(i,:)-miu))),2))';
%     end
	D=[ones(N,1),inputdata,fi];
end

function P=GenerateP(inputdata,N,k,miu)
    %Generating the P matrix;
		D=GenerateD(inputdata,N,k,miu);
		P=eye(N)-D*pinv(D'*D)*D';
end

function [kre,miure,Pre,flag]=BirthMove(inputdata,outputdata,mininput,maxinput,k,miu,N,Ccoe,P)
    [~,c]=size(outputdata);
	[~,d]=size(inputdata);
	flag=false;
	area=prod(max(miu)-min(miu));
	miutemp=mininput+(maxinput-mininput).*rand(1,d);
	miunext=[miu;miutemp];
	Pnext=GenerateP(inputdata,N,k+1,miunext);
	r=area*exp(-Ccoe)/(k+1);
	for j=1:c
		r=r*((outputdata(:,j)'*P*outputdata(:,j))/(outputdata(:,j)'*Pnext*outputdata(:,j)))^(N/2);
	end
	A=min(1,r);
	if rand()<A
		kre=k+1;
		miure=miunext;
        Pre=Pnext;
		flag=true;
	else
		kre=k;
		miure=miu;
        Pre=P;
	end
end

function [kre,miure,Pre,flag]=DeathMove(inputdata,outputdata,k,miu,N,Ccoe,P)
	[~,c]=size(outputdata);
	flag=false;
	area=prod(max(miu)-min(miu));
	deletek=max(ceil(rand()*k),1);
	miunext=miu;
	miunext(deletek,:)=[];
	Pnext=GenerateP(inputdata,N,k-1,miunext);
	r=k*exp(Ccoe)/area;
	for j=1:c
		r=r*power(outputdata(:,j)'*P*outputdata(:,j)/(outputdata(:,j)'*Pnext*outputdata(:,j)),N/2);
	end
	A=min(1,r);
	if rand()<A
		kre=k-1;
		miure=miunext;
        Pre=Pnext;
		flag=true;
	else
		kre=k;
		miure=miu;
        Pre=P;
	end
end

function [kre,miure,Pre,flag]=SplitMove(inputdata,outputdata,k,miu,N,Ccoe,zeta,P)
	[~,c]=size(outputdata);
	flag=false;
	splitk=max(ceil(rand()*k),1);
	splitmiu=miu(splitk,:);
	ums=rand();
	splitmiu1=splitmiu-ums*zeta;
	splitmiu2=splitmiu+ums*zeta;
	miunext=miu;
	miunext(splitk,:)=[];
	miunext=[miunext;splitmiu1;splitmiu2];
	Pnext=GenerateP(inputdata,N,k+1,miunext);
	r=k*zeta*exp(-Ccoe)/(k+1);
		for j=1:c
		r=r*power((outputdata(:,j)'*P*outputdata(:,j))/(outputdata(:,j)'*Pnext*outputdata(:,j)),N/2);
	end
	A=min(1,r);
	if rand()<A
		kre=k+1;
		miure=miunext;
        Pre=Pnext;
		flag=true;
	else
		kre=k;
		miure=miu;
        Pre=P;
	end
end

function [kre,miure,Pre,flag]=MergeMove(inputdata,outputdata,k,miu,N,Ccoe,zeta,P)
	[~,c]=size(outputdata);
	flag=false;
	mergek=max(ceil(rand()*k),1);
	mergemiu=miu(mergek,:);
	miunext=miu;
	miunext(mergek,:)=[];
	distan=sqrt(sum((miunext-mergemiu).*(miunext-mergemiu),2));
	[~,mergek]=min(distan);
	mindismiu=miunext(mergek,:);
	if sqrt(sum((mindismiu-mergemiu).*(mindismiu-mergemiu),2))<(2*zeta)
		miutemp=1/2*(mergemiu+mindismiu);
		miunext(mergek,:)=miutemp;
		Pnext=GenerateP(inputdata,N,k-1,miunext);
		r=k*exp(Ccoe)/(k-1)/zeta;
		for j=1:c
			r=r*power(outputdata(:,j)'*P*outputdata(:,j)/(outputdata(:,j)'*Pnext*outputdata(:,j)),N/2);
		end
		A=min(1,r);
		if rand()<A
			kre=k-1;
			miure=miunext;
            Pre=Pnext;
			flag=true;
		else
			kre=k;
			miure=miu;
            Pre=P;
		end
	else
		kre=k;
		miure=miu;
        Pre=P;
	end
end

function [kre,miure,Pre,flag]=UpdataMove(inputdata,outputdata,mininput,maxinput,k,miu,N,Pcoe,P)
	[~,c]=size(outputdata);
    [~,d]=size(inputdata);
    flag=false;
	miutemp=mininput+(maxinput-mininput).*rand(1,d);
	updatak=max(ceil(k*rand()),1);
    miunext=miu;
	miunext(updatak,:)=miutemp;
	Pnext=GenerateP(inputdata,N,k,miunext);
	r=exp(-Pcoe);
	for j=1:c
		r=r*power(outputdata(:,j)'*P*outputdata(:,j)/(outputdata(:,j)'*Pnext*outputdata(:,j)),N/2);
	end
	A=min(1,r);
	if rand()<A
		kre=k;
		miure=miunext;
        Pre=Pnext;
		flag=true;
	else
		kre=k;
		miure=miu;
        Pre=P;
	end
end

function Arjsa=SimAnneal(k,knext,miunext,P,C,T,inputdata,outputdata)
    [N,~]=size(inputdata);
    [~,c]=size(outputdata);
	Pnext=GenerateP(inputdata,N,knext,miunext);
	r=exp(-C*k);
	rnext=exp(-C*knext);
	for j=1:c
		r=r*power(outputdata(:,j)'*P*outputdata(:,j),N/2);
		rnext=rnext*power(outputdata(:,j)'*Pnext*outputdata(:,j),N/2);
	end
	Arjsa=min(1,power((rnext/r),(1/T-1)));
end