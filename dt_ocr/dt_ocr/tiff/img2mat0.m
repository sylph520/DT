%nntest

function outMat1 = img2mat0()
outMat = [];
for i  = 0:7
    for j = 0:19
         foldstr =['D:\\data\\graph-DB\\nt4\\charImgs\\tiff\\charImgsep\\new\\', int2str(i)];
         imgstr = int2str(j);
         pathstr = [foldstr, '\\',imgstr,'.tiff'];
         tmpImg = imread(pathstr); tmpVec = reshape(tmpImg,[1,225]); tmpVec = [tmpVec, i];
         outMat = [outMat; tmpVec];
    end
end
    