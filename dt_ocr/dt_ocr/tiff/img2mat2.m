function [X,Y] = img2mat2() 
X = [];Y=[];
for i  = 0:7
         foldstr =['D:\\data\\graph-DB\\nt4\\charImgs\\tiff\\charImgsep\\new2\\', int2str(i),'\\'];
         files = dir(foldstr);
         lens = length(files) - 2;
         for j = 1:lens
            imgstr = int2str(j);
            pathstr = [foldstr,imgstr,'.tiff'];
            tmpImg = imread(pathstr); tmpVec = reshape(tmpImg,[1,225]);
            tmpOut = (0:7) == i;
            X = [X; tmpVec];
            Y= [Y; tmpOut];
         end
end
