convert8=1;
%Get Image folder
foldername = uigetdir;
%Iterate over all files in this folder
allfiles = dir(foldername);
for(i=1:numel(allfiles))
    %Check if file ends with .tif or .tiff
    ind = strfind([foldername '/' allfiles(i).name],'.ti');
    if(numel(ind) > 0)
        %Load image and save it as .tiff
        img = imread([foldername '/' allfiles(i).name]);
        if(convert8)
            %Convert image to 8 bit
            img = double(img)./4095;
         %   img = imadjust(img, [double(0);double(1)], [double(0);double(1)]);
            img = uint8(img.*255);
        end
        %Delete old file
        delete([foldername '/' allfiles(i).name]);
        %Replace .tif(f) with .png in filename
        allfiles(i).name = strrep(allfiles(i).name,'.tiff','.png');
        allfiles(i).name = strrep(allfiles(i).name,'.tif','.png');
        imwrite(img,[foldername '/' allfiles(i).name]);
    end
end