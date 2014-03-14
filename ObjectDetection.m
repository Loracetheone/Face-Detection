function [Objects, Picture,probObj,Ratio] = ObjectDetection(Picture,modelo)
    if (size(Picture,2) > size(Picture,1))
        Ratio = size(Picture,2) /350;
    else
        Ratio = size(Picture,1) /350;
    end
    Picture = imresize(Picture, [size(Picture,1) size(Picture,2) ]/ Ratio);

[fil,col,k]=size(Picture);
if k>1
    im = double(Picture);
    skinprob = computeSkinProbability(im);
    ImPiel = (skinprob > 0)+0;
    B = [0 1 0;1 1 1;0 1 0];
    IM2=ImPiel;%imdilate(ImPiel,B);  
    Mat=imdilate(imdilate(imdilate(imdilate(imdilate(imdilate(IM2,B),B),B),B),B),B);    
  %  figure;hold on; imshow(Mat,[])

    Picture=double(rgb2gray(Picture)); 
    [Objects,probObj] = HogObjectDetection(Mat,Picture,modelo);
else    
    Picture=im2double(Picture);
    Mat=0;
    [Objects,probObj] = HogObjectDetection(Mat, Picture,modelo);
end


