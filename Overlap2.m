function [band,probObj,probR,Objects]=Overlap2(Result,Objects,tam,probR,probObj)
band2=0;
band=1;
xm=0;
porInter=10;
if tam>0
  for m=1:tam    
                %area de interseccion de ambos rectangulos
                area3 = rectint(Objects(m,:),Result);                
                area1= Objects(m,3)*Objects(m,4);%Objects
                area2= Result(1,3)*Result(1,4);%Result                
                por1= area3*100/area2;
                por2=area3*100/area1;                                
                if por2 > porInter || por1 > porInter
                    %significa q si se intersecta mas del porcentage se descarta
                    %a la ventana
                    band=0;                  
                    %HAY interseccion entre ambos rectg
                    band2= band2+1;                    
                     %verificar la intersccion junto con la probabilidad
                    if probObj(m,1) < probR(1,1) %si este objeto ya se intersecto con mas de 3 ents el obj actual se considera falso positivo                                                               
                        xm=m;                    
                    end            
                end
  end
  
  %Hacemos el intercambio
  if band2==1 && band==0 && xm>0
   Objects(xm,1)= Result(1,1);% y
   Objects(xm,2)= Result(1,2);% x
   Objects(xm,3)= Result(1,3);
   Objects(xm,4)= Result(1,4);   
   probObj(xm,1)= probR(1,1); 
  end  
end
