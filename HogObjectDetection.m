function [Objects,probObj]= HogObjectDetection(Mat,Images,modelo)
% Mat=0;
Objects=zeros(20,4); 
probObj=zeros(20,1); %probabilidad de cada obj+bandera para ver sino se intersecto con ninguno
n=0; 
% orden=[1 0.9 0.8 0.7 0.6 0.5 0.4 0.3]; %orden original
orden=[0.19 0.2 0.6 0.3 0.4 0.8 1.2];
% orden= [0.6 0.3 0.5 0.4 0.8 1.2]; % verificar esto!!
tam=7;
Pic = cell(1,1);
Matt=cell(1,1);
%Guardamos a todas la piramide de imagenes
Picture = imresize(Images,orden(1));
[w,h]=size(Picture);
if w<55 || h<55
   ij=2;
else
    ij=1;
end
 parfor i=ij:tam  
     Picture = imresize(Images,orden(i)); 
     [w,h]=size(Picture);
     if w>420 || h>420
        orden(i)=orden(i)-0.1;
        Picture = imresize(Images,orden(i));    
    end
    if Mat==0        
    else
        Mat1=imresize(Mat,orden(i));
        Matt{i}=logical(uint8(Mat1));        
    end
    Pic{i}=Picture; 
 end

%-------Extraemos las subventanas de cada imagen-----------------
for i=ij:tam
    mat_array = cell(1,1);
    [w,h]=size(Pic{i});  
    va11= w-50;
    va12=h-50;
    Ex=[];
    Ey=[];    
    [X,Y]= meshgrid(1:5:va11,1:6:va12);
    X=X(:);
    Y=Y(:);     
    if Mat==0  
        %primero cogemos los indices de los cuadrados        
        d=size(X,1);        
        parfor id=1:d  
                subimg= Pic{i}(X(id):X(id)+49, Y(id):Y(id)+49);
                hog = HoG(subimg);        
                y1=reshape(hog',900,1)'; 
                mat_array{id}=y1;
        end
        Ej=[X,Y];
        X1=vertcat(mat_array{:});        
    else     
        d=size(X,1);   
        area=sum(sum(Matt{i}==1));%Area con respecto al tamaño de la imagen
        parfor id=1:d            
                submat= Matt{i}(X(id):X(id)+49, Y(id):Y(id)+49);
                sumI=sum(sum(submat==1));
                umb1=sumI/area;
                umb2=sumI*0.04; %area con respecto al tamaño de la ventana 50x50                                       
                    if sumI>1193 && umb1>0.033 && umb2>50
                          subimg= Pic{i}(X(id):X(id)+49, Y(id):Y(id)+49);
                          hog = HoG(subimg);        
                          y1=reshape(hog',900,1)'; 
                          mat_array{id}=y1; 
                          Ex=[Ex,X(id)];
                          Ey=[Ey,Y(id)];
                    end
        end                   
        Ej=[Ex',Ey'];
%         mat_array(~cellfun('isempty',mat_array));  
        X1=vertcat(mat_array{:});
    end    
    Y1=zeros(size(X1,1),1);            
    [~, ~, dec_values] =svmpredict(Y1,X1,modelo.model);      
    VR=find(dec_values>0.4);    
    proR=dec_values(VR(:));
    
    if ~isempty(VR)  
       B=ones(size(VR,1),1)*50;
       Result=floor([Ej(VR,2) Ej(VR,1) B(:) B(:)]/(orden(i)));  
       for k=1:size(VR,1)
        %  subI= Pic{i}(Result(2):Result(2)+49, Result(1):Result(1)+49);
        %  [BBOX,ORIENTATION, SCORE] = step(Pic{i},Result(k,:));
           
           [band,probObj,probR2,Objects]=Overlap2(Result(k,:),Objects,n,proR(k,1),probObj);            
           if band==1  
              n=n+1;
              Objects(n,:)=Result(k,:);  
              probObj(n,1)= probR2; %probabilidad   probR          
           end
       end
    end
end
% Crop the initial array with Objects detected
Objects=Objects(1:n,:);
