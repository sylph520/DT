
function Mat = img2mat(pathstr)
        tmpImg = imread(pathstr); 
        tmpVec = reshape(tmpImg,[1,225]); 
         Mat = tmpVec;
end
