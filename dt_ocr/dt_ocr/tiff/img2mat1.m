
function [X, Y] = img2mat1()
X = [];Y=[];
for i  = 0:7
    for j = 0:19
         foldstr =['D:\\data\\graph-DB\\nt4\\charImgs\\tiff\\charImgsep\\new\\', int2str(i)];
         imgstr = int2str(j);
         pathstr = [foldstr, '\\',imgstr,'.tiff'];
         tmpImg = imread(pathstr); tmpVec = reshape(tmpImg,[1,225]);
         tmpOut = (0:7) == i;
         X = [X; tmpVec];
         Y= [Y; tmpOut];
    end
end