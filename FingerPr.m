

%%%%%%%%%%%%The main code contains 4 stages of search (Hirarchical Search)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all
warning off;
'running'
tic
%%%% Getting the images and Pre-processing
TempFile= '104_4.tif';
TemplateROI= (imread(['C:\MyResearch\fingerprint\DB3_B_ROItemplates\' TempFile]));
TemplateROI= double(TemplateROI);
AcceptedRatio= 4;
dirname = dir('C:\MyResearch\fingerprint\DB3_B_ROIinputs\*.tif');
Count= length(dirname);
for k= 1:Count
        '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        InputFile= dirname(k).name
        Input1ROI=   (imread(['C:\MyResearch\fingerprint\DB3_B_ROIinputs\' InputFile]));
        Input1ROI= double(Input1ROI);
        [NN MM]= size(TemplateROI);
        [SS1 SS2]= size(Input1ROI);
        SS= min(SS1,SS2);
 
        num= 1
 
        tetas0= 4;
        xf= 2+(MM/AcceptedRatio);
        xe= MM+SS-2-(MM/AcceptedRatio);
        yf= 2+(NN/AcceptedRatio);
        ye= NN+SS-2-(NN/AcceptedRatio);
        xs0= (xe-xf)/10;
        ys0= (ye-yf)/10;
        
%         step=[-6 tetas0 6; -100 xs0 100;-100 ys0 100];% [teta1 tetaStep teta2; x1 xStep x2; y1 yStep y2]
        step=[-6 tetas0 6; xf xs0 xe;yf ys0 ye];% [teta1 tetaStep teta2; x1 xStep x2; y1 yStep y2]
        WdRatio= 1;
        NumBin= 18;
%         ResizeCoef=.8;
        ResizeCoef=1;
        [NMI x y teta duration V]= FindTransform(num,TemplateROI,Input1ROI,step,WdRatio,NumBin,ResizeCoef)
        % if V==0
        %     'These Parameters Do Not Work!'
        % end
        result(num).tem= TempFile;
        result(num).inp= InputFile;
        result(num).step= step;
        result(num).WdRatio= WdRatio;
        result(num).NumBin= NumBin;
        result(num).ResizeCoef= ResizeCoef;
        result(num).results= [NMI x y teta duration V];
 
        num= 2
        xs1= 8;
        ys1= 8;
        tetas1= 3;
        step=[max(teta-(tetas0/2),-6) tetas1 min(teta+(tetas0/2),6); max(x-(xs0/2),xf) xs1 min(x+(xs0/2),xe);max(y-(ys0/2),yf) ys1 min(y+(ys0/2),ye)];% [teta1 tetaStep teta2; x1 xStep x2; y1 yStep y2]
        % WdRatio= 1;
        % NumBin= 18;
        ResizeCoef=1;
        [NMI x y teta duration V]= FindTransform(num,TemplateROI,Input1ROI,step,WdRatio,NumBin,ResizeCoef)
        % if V==0
        %     'These Parameters Do Not Work!'
        % end
        result(num).tem= TempFile;
        result(num).inp= InputFile;
        result(num).step= step;
        result(num).WdRatio= WdRatio;
        result(num).NumBin= NumBin;
        result(num).ResizeCoef= ResizeCoef;
        result(num).results= [NMI x y teta duration V];
 
        % pause;
        num= 3
        xs2= 3;%floor(xs1/2);
        ys2= 3;%floor(ys1/2);
        tetas2= 2;%floor(tetas1/2);
        step=[max(teta-tetas1,-6) tetas2 min(teta+tetas1,6); max(x-(xs1/2),xf) xs2 min(x+(xs1/2),xe);max(y-(ys1/2),yf) ys2 min(y+(ys1/2),ye)];% [teta1 tetaStep teta2; x1 xStep x2; y1 yStep y2]
        [NMI x y teta duration V]= FindTransform(num,TemplateROI,Input1ROI,step,WdRatio,NumBin,ResizeCoef)
        result(num).tem= TempFile;
        result(num).inp= InputFile;
        result(num).step= step;
        result(num).WdRatio= WdRatio;
        result(num).NumBin= NumBin;
        result(num).ResizeCoef= ResizeCoef;
        result(num).results= [NMI x y teta duration V];
 
        % pause;
        num= 4
        xs3= 1;%floor(xs1/2);
        ys3= 1;%floor(ys1/2);
        tetas3= 1;%floor(tetas1/2);
        step=[max(teta-tetas2,-6) tetas3 min(teta+tetas2,6); max(x-(xs2),xf) xs3 min(x+(xs2),xe);max(y-(ys2),yf) ys3 min(y+(ys2),ye)];% [teta1 tetaStep teta2; x1 xStep x2; y1 yStep y2]
        [NMI x y teta duration V]= FindTransform(num,TemplateROI,Input1ROI,step,WdRatio,NumBin,ResizeCoef)
        result(num).tem= TempFile;
        result(num).inp= InputFile;
        result(num).step= step;
        result(num).WdRatio= WdRatio;
        result(num).NumBin= NumBin;
        result(num).ResizeCoef= ResizeCoef;
        result(num).results= [NMI x y teta duration V];
        
        save(['C:\mtemp\' date 'FingerPrint' TempFile '_to_' InputFile 'Shab.mat'],'result')
end
time= toc
 
 
 


%%%%%%%%%%% a function which finds the best alignment (transformation) for
%%%%%%%%%%% a pair of images
function [NMImax x0 y0 teta0 myTime validity]= FindTransform(Num,TemplateROI,input1ROI,step,WdRatio,NumBin,ResizeCoef)
%%%%%  Validity=1: valid    Validity=0: more than one NMI=2     Validity=2:
%%%%%  more than one maximum
 
%%%%%%%%%%Finding the Gradient image for Template
        tic
        step =step';
        TemplateROI= imresize(TemplateROI,ResizeCoef);
        input1ROI= imresize(input1ROI,ResizeCoef);
        [N1 M1]= size(TemplateROI);
        Kx= [1 0 -1; 2 0 -2; 1 0 -1];
        Gx= conv2(Kx,TemplateROI);
        Gy= conv2(Kx',TemplateROI);
        Gx1= Gx(2:N1+1,2:M1+1);
        Gy1= Gy(2:N1+1,2:M1+1);
        W= 16; %The size of the block to find the Orientation Field
        %%%%%%%%%%%%% Finding Orientation Field of Template
        Sx= zeros(N1+W,M1+W);
        Sx(1+W/2:N1+W/2,1+W/2:M1+W/2)= Gx1;
        Sy= zeros(N1+W,M1+W);
        Sy(1+W/2:N1+W/2,1+W/2:M1+W/2)= Gy1;
 
        for i= (1+W/2):(N1+W/2)
            for j= (1+W/2):(M1+W/2)
                Teta1(i-W/2,j-W/2)=  .5*atan2((sum(sum(2*Sx(i-W/2:i+W/2,j-W/2:j+W/2).*Sy(i-W/2:i+W/2,j-W/2:j+W/2)))),(sum(sum(Sx(i-W/2:i+W/2,j-W/2:j+W/2).^2-Sy(i-W/2:i+W/2,j-W/2:j+W/2).^2))));
            end
        end
        %%%%%%%%%%% filtering the Orientation Field image
        Fx= 4;
        Fy= 4;
        TetaD1= Teta1*180/pi;
        TetaDF1= imfilter(TetaD1,ones(Fx,Fy))/(Fx*Fy);
        %%%%%%%%%%%%%%5 Finding the Direction Feature image
        Wd=floor(8*WdRatio);%%% the size of the block of averaging
        kn= floor(N1/Wd);
        km= floor(M1/Wd);
        for i= 1:kn
            for j= 1:km
                DirectionFeature1(i,j)= mean(mean((TetaDF1((i-1)*Wd+1:(i)*Wd,(j-1)*Wd+1:(j)*Wd))));
            end
        end
        Com1=0;
        Com2=2;
        se = strel('disk',2);
 
        NMIArray= [];
                    %%%%%%%%%%%%%%%%%%% Defining the Transformation
        TheFlag=1;
        k0=0;
        for teta= step(1):step(2):step(3)
            k1=0;
            k0= k0+1;
            for x= step(4):step(5):step(6)
                k2=0;
                k1= k1+1;
            for y= step(7):step(8):step(9)
                %%%%%%%%%%%%%%%
                if x==200 && y==200 && teta==-5
                    'ali'
                end
                %%%%%%%%%%%%
                    k2= k2+1;
                    imOver= myOverlap(input1ROI,x,y,teta,N1,M1);
                    if (size(imOver)==[0 0])
                        continue;
                    end
                    for i= 1:N1
                        for j= 1:M1
                            if imOver(i,j)>=0
                                imOverMask(i,j)=1;
                            else
                                imOverMask(i,j)=0;
                            end
 
                        end
                    end
 
                                %%%%%%%%%%%%% Finding the Direction Feature mask
                    for i= 1:kn
                        for j= 1:km
                            DirectionMask2(i,j)= sum(sum((imOverMask((i-1)*Wd+1:(i)*Wd,(j-1)*Wd+1:(j)*Wd))));
                        end
                    end
 
                    DirectionMaskBinary2= zeros(kn,km);
                    for i= 1:kn
                        for j= 1:km
                            if DirectionMask2(i,j)>=(Wd^2/1)
                                DirectionMaskBinary2(i,j)= 1;
                            end
                        end
                    end
 
 
                    DirectionMaskBinary2= imerode(DirectionMaskBinary2,se);
                    if sum(sum(DirectionMaskBinary2))<(kn*km/4)
                        FF(k0,k1,k2)= 0;
                        continue;
                    end
                    %%%%%%%%%%%% Now two same size images are ready to find the features
                    %%%%% Finding the Orientation Field
 
 
                    Gx= conv2(Kx,imOver);
                    Gy= conv2(Kx',imOver);
                    Gx2= Gx(2:N1+1,2:M1+1);
                    Gy2= Gy(2:N1+1,2:M1+1);
 
 
                    Gx2= Gx2.*imOverMask;
                    Gy2= Gy2.*imOverMask;
 
 
 
                    %%%%%%%%%%%%% Finding Orientation Field of Image2
                    Sx= zeros(N1+W,M1+W);
                    Sx(1+W/2:N1+W/2,1+W/2:M1+W/2)= Gx2;
                    Sy= zeros(N1+W,M1+W);
                    Sy(1+W/2:N1+W/2,1+W/2:M1+W/2)= Gy2;
                    %%%%%%%%%%%5 Finding the Orientation field of Input
                    for i= (1+W/2):(N1+W/2)
                        for j= (1+W/2):(M1+W/2)
                            Teta2(i-W/2,j-W/2)= .5*atan2((sum(sum(2*Sx(i-W/2:i+W/2,j-W/2:j+W/2).*Sy(i-W/2:i+W/2,j-W/2:j+W/2)))),(sum(sum(Sx(i-W/2:i+W/2,j-W/2:j+W/2).^2-Sy(i-W/2:i+W/2,j-W/2:j+W/2).^2))));
                        end
                    end
                    Teta2NaN= isnan(Teta2);
                    [t1,t2]= size(Teta2);
                    for i= 1:t1
                        for j= 1:t2
                            if   (Teta2NaN(i,j)==1)
                                Teta2(i,j)=0;
                            end
                        end
                    end
                    %%%%%%5 Filtering the Orientation Field image
                    TetaD2= Teta2*180/pi;
                    TetaDF2= imfilter(TetaD2,ones(Fx,Fy))/(Fx*Fy);
                    TetaDF2= TetaDF2.*imOverMask;
                    %%%%%%%%%%%%%%%% Finding the Direction Feature
 
                    for i= 1:kn
                        for j= 1:km
                            DirectionFeature2(i,j)= mean(mean((TetaDF2((i-1)*Wd+1:(i)*Wd,(j-1)*Wd+1:(j)*Wd))));
                        end
                    end
 
                    %%%%%%%%%%%%5 Calculating the NMI in the overlapped area
 
                    NMI= MIcalc(DirectionFeature1+90,DirectionFeature2+90,DirectionMaskBinary2,NumBin);
 
                    NMIArray= [NMIArray [NMI;teta;x;y;k0;k1;k2]];
                    FF(k0,k1,k2)= NMI;
                    if (NMI==2)
                        if TheFlag==0
                            validity=0;
                            x0= 0;
                            y0= 0;
                            teta0= 0;
                            NMImax=0;
                            myTime=0;
                            return;
                        else
                            TheFlag=0;
                        end
                    end
                end
            end
        end
        if prod(size(NMIArray))==0
            validity=2;
            x0= 0;
            y0= 0;
            teta0= 0;
            NMImax=0;
            myTime=0;
            return; 
        end;
        indMax= find(max(NMIArray(1,:))==NMIArray(1,:));
        if size(indMax)>1
            validity=2;
            x0= 0;
            y0= 0;
            teta0= 0;
            NMImax=0;
            myTime=0;
            return;
        else
            validity=1;
        end
        indMaxH= indMax(1);
        clear indMax;
        indMax= indMaxH;
        
        imOverMax= myOverlap(input1ROI,NMIArray(3,indMax),NMIArray(4,indMax),NMIArray(2,indMax),N1,M1);
%         figure
%         plot(NMIArray(1,:));
%         title(num2str(Num))
        
        figure
        subplot(1,2,1)
        imshow(uint8(immerge(imOverMax,uint8(TemplateROI),.5)));
        title(num2str(Num))
        h(:,:)= FF(NMIArray(5,indMax),:,:);
        
        subplot(1,2,2)
        surf((1:k2),(1:k1),h), hold on
        plot3(NMIArray(7,indMax),NMIArray(6,indMax),NMIArray(1,indMax),'*r','MarkerSize',12);
        zlabel('NMI');
        xlabel('Y');
        ylabel('X');
        title(num2str(Num))
        x0= NMIArray(3,indMax);
        y0= NMIArray(4,indMax);
        teta0= NMIArray(2,indMax);
        NMImax= NMIArray(1,indMax);
        myTime= toc/60;
end







%%% a function which constructs the transformed image and finds the
%%% overlapp and the corresponding Mask
function imOver= myOverlap(im,x,y,teta,N,M)
%     a= min(min(im));
%     if a==0
%         im= im+1;
%     end
    im= im+1;
%     R= double(imrotate(im,teta));
    R= (imrotate(im,teta));
    R=R-1;
    [s1,s2]= size(R);
    IM= (-1*ones(s1+2*N,s2+2*M));
    IM(N+1:N+s1,M+1:M+s2)= R;
%     imOver= IM(s1+1+y:s1+N+y,s2+1+x:s2+M+x);
%     yo= (s1+N)/2;
%     xo= (s2+M)/2;
    yo= 0;
    xo= 0;
 
    if ((yo+1+y)>0 && (xo+1+x)>0&& (yo+N+y)<=(s1+2*N)&&(xo+M+x)<=(s2+2*M))
        imOver= IM(yo+1+y:yo+N+y,xo+1+x:xo+M+x);
    else
        imOver= [];
    end
end


%%%%%%%%%%% This function finds the NMI over the overlapped area
function NMI= MIcalc(Im1,Im2,Mask2,n)
% assuming that element values of IM1 & Im2 are from 0 to 180
% Im1 & Im2 have same sizes.
    [s1 s2]= size(Im1);
%     n= 60;%%% the number of bins
    r= 180/(n-1);
    %%%%%%%%%%%%Finding I(A,B)
    FS= zeros(n,n);
    for i= 1:s1
        for j= 1:s2
            if Mask2(i,j)==1
                IX= floor(Im1(i,j)/r)+1;
                IY= floor(Im2(i,j)/r)+1;
                FS(IX,IY)= FS(IX,IY)+1;  
            end
        end
    end
    T= sum(sum(FS));
    FSp= FS/T;
    G= 0;
    for i= 1:n
        for j= 1:n
            if FSp(i,j)~=0
                G= G - FSp(i,j)*log(FSp(i,j));
            end
        end
    end
    %%%%%%%%%%%%% Finding I(A)
    FS= zeros(n);
    for i= 1:s1
        for j= 1:s2
            if Mask2(i,j)==1
                IX= floor(Im1(i,j)/r)+1;
                FS(IX)= FS(IX)+1; 
            end
        end
    end
    T= sum(sum(FS));
    FSp= FS/T;
    A= 0;
    for j= 1:n
        if FSp(j)~=0
            A= A - FSp(j)*log(FSp(j));
        end
    end
    %%%%%%%%%%%%% Finding I(B)
    FS= zeros(n);
    for i= 1:s1
        for j= 1:s2
            if Mask2(i,j)==1
                IX= floor(Im2(i,j)/r)+1;
                FS(IX)= FS(IX)+1;  
            end
        end
    end
    T= sum(sum(FS));
    FSp= FS/T;
    B= 0;
    for j= 1:n
        if FSp(j)~=0
            B= B - FSp(j)*log(FSp(j));
        end
    end 
    
    NMI= (A+B)/G;  %%%%NMI
end

