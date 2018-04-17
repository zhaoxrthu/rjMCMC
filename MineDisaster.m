function lamdare=MineDisaster(daylen,happendays)  
    %paraments from paper
    Iiter=20000;
    kmax=30;lamda=3;alpha=1;beta=200;
    %initializing thresholds 
    bk=zeros(kmax,1);dk=zeros(kmax,1);
    for i=1:kmax-1
        bk(i)=min(1,poisspdf(i,lamda)/poisspdf(i-1,lamda));
        dk(i+1)=min(1,poisspdf(i-1,lamda)/poisspdf(i,lamda));
    end
        rate=1/max(bk+dk);
        bk=rate.*bk;    dk=rate.*dk;
        hk=(1-bk-dk)/2;    pk=hk;
        hk(1)=1-bk(1)-dk(1);    pk(1)=0;
    %%setting the initial value
    k=3;
    divide=ceil([daylen/4,daylen/2,3*daylen/4]');
    E=alpha/beta*[1;1;1;1];

    for i=1:Iiter
        u=rand();
        if u<hk(k)
            [Ere,dividere,kre]=height(k,E,divide,happendays,beta,daylen);
        else
            if u<hk(k)+pk(k)
            [Ere,dividere,kre]=position(k,E,divide,happendays,daylen);
            else
                if u<hk(k)+pk(k)+bk(k)
                [Ere,dividere,kre]=birth(divide,E,k,happendays,lamda,alpha,beta,bk,dk,daylen);
            else
                [Ere,dividere,kre]=death(divide,E,k,happendays,lamda,alpha,beta,bk,dk,daylen);
                end
            end
        end
        E=Ere;
        divide=dividere;
        k=kre;  
    end 
    lamdare=getlamda(divide,E,daylen);
end

function llh=LLHR(L,x1,lamda1,x2,lamda2,y)
    lamdanext1=ones(length(y),1);
    if isempty(x1)
        lamdanext1=lamda1*lamdanext1;    
    else
        lamdanext1(y<x1(1))=lamda1(1);    
        for i=1:length(x1)-1
            lamdanext1(y>=x1(i)&y<x1(i+1))=lamda1(i+1);        
        end
        lamdanext1(y>=x1(end))=lamda1(end);    
    end

    lamdanext2=ones(length(y),1);
    if isempty(x2)
        lamdanext2=lamda2*lamdanext2;    
    else
        lamdanext2(y<x2(1))=lamda2(1);    
        for i=1:length(x2)-1
            lamdanext2(y>=x2(i)&y<x2(i+1))=lamda2(i+1);         
        end
        lamdanext2(y>=x2(end))=lamda2(end);   
    end
    if(isempty(x1))
        L1=sum(log(lamdanext1))-sum(diff([0;L]).*lamda1);
    else
        L1=sum(log(lamdanext1))-sum(diff([0;x1;L]).*lamda1);
    end
    if(isempty(x2))
        L2=sum(log(lamdanext2))-sum(diff([0;L]).*lamda2);
    else
        L2=sum(log(lamdanext2))-sum(diff([0;x2;L]).*lamda2);
    end
    llh=exp(L1-L2);
end

function [Ere,dividere,kre]=birth(divide,E,k,y,lamda,alpha,beta,bk,dk,daylen)
    dividenext=1+max(ceil((daylen-3)*rand()),1);
    j=find(dividenext<divide,1)-1;
    if isempty(j)
        j=k;
    end
    dividenew=[divide(1:j);dividenext;divide(j+1:end)];
    time=[1;divide;daylen];
    u=rand();
    lamda1new=E(j+1)*(u/(1-u))^((time(j+2)-dividenext)/(time(j+2)-time(j+1)));
    hr_new=lamda1new*(1-u)/u;
    h_new=[E(1:j);lamda1new;hr_new;E(j+2:end)];
    llh=LLHR(daylen,dividenew,h_new,divide,E,y);
    pllh=2*lamda*(2*k+3)/...
        daylen^2*(dividenext-time(j+1))*(time(j+2)-dividenext)/(time(j+2)-time(j+1))...
        *beta^alpha/gamma(alpha)*(lamda1new*hr_new/E(j+1))^(alpha-1)...
        *exp(-beta*(lamda1new+hr_new-E(j+1)));
    P=dk(k+2)*daylen/bk(k+1)/(k+1);
    J=(lamda1new+hr_new)^2/E(j+1);
    if rand()<llh*pllh*P*J
        Ere=h_new;
        dividere=dividenew;
        kre=k+1;
    else
        Ere=E;
        dividere=divide;
        kre=k;
    end
end

function [Ere,dividere,kre]=death(divide,E,k,y,lamda,alpha,beta,bk,dk,daylen)
    j=max(ceil(k*rand()),1);
    time=[1;divide;daylen];            
    Enext=exp((time(j+1)-time(j))/(time(j+2)-time(j))*log(E(j))...
        +(time(j+2)-time(j+1))/(time(j+2)-time(j))*log(E(j+1)));            
    dividenext=[divide(1:j-1);divide(j+1:end)];
    h_new=[E(1:j-1);Enext;E(j+2:end)];
    llh=LLHR(daylen,dividenext,h_new,divide,E,y);
    pllh=daylen^2/2/lamda/(2*k+1)*(time(j+2)-time(j))...
        /(time(j+1)-time(j))/(time(j+2)-time(j+1))*gamma(alpha)/beta^alpha*...
        (Enext/E(j)/E(j+1))^(alpha-1)*exp(beta*(E(j)+E(j+1)-Enext));
    P=bk(k)*k/dk(k+1)/daylen;
    J=(Enext/(E(j)+E(j+1))^2);
    if rand() <llh*pllh*P*J        
        Ere=h_new;
        dividere=dividenext;
        kre=k-1;
    else
        Ere=E;
        dividere=divide;
        kre=k;
    end
end

function [Ere,dividere,kre]=height(k,E,divide,y,beta,daylen)
    knext=max(ceil(rand()*(k+1)),1);
    egammanex=E(knext)*exp(rand()-0.5);
    llhr=LLHR(daylen,divide,[E(1:knext-1);egammanex;E(knext+1:end)],divide,E,y);
    acceptance=llhr*(egammanex/E(knext))*exp(-beta*(egammanex-E(knext)));
    Ere=E;
    if rand()<acceptance
        Ere(knext)=egammanex;
    end
    dividere=divide;
    kre=k;
end

function [Ere,dividere,kre]=position(k,E,divide,y,daylen)
    dnext=max(ceil(rand()*(k)),1);
    time=[1;divide;daylen];
    s_l=time(dnext);
    s_r=time(dnext+2);
    temp=[s_l+1,s_r-1];
    dividenext=randi(temp);
    ratio=LLHR(daylen,[divide(1:dnext-1);dividenext;divide(dnext+1:end)],E,divide,E,y);
    dividere=divide;
    if rand()<(s_r-dividenext)*(dividenext-s_l)/(s_r-divide(dnext))/(divide(dnext)-s_l)*ratio
        dividere(dnext)=dividenext; 
    end
    Ere=E;
    kre=k;
end

function lamdare=getlamda(divide,E,daylen)
    lamdatemp=zeros(1,daylen);          
    %how many group do days be divided   
    divnum=length(divide);
    for d=1:divnum
        for i=1:daylen
            if(d==1)
                if(i<divide(d))
                    lamdatemp(i)=E(d);
                end
            else
                if(i<divide(d) && i>=divide(d-1))
                    lamdatemp(i)=E(d);
                end
            end
        end
    end
    lamdatemp(divide(divnum):end)=E(divnum+1);
    lamdare=lamdatemp;
end