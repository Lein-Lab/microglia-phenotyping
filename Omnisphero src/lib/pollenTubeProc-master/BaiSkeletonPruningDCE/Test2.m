%I=double(imread('carriage-17.GIF'));
I = imread('4.tif');
%I = imresize(I,.75);
bw = im2bw(I);
%bw= 1-bw; %the shape must be black, i.e., values zero.

%[bw,I0,x,y,x1,y1,aa,bb]=div_skeleton_new(4,1,bw,4);
skel=div_skeleton_new(4,1,bw,8);
skel=parsiSkel(skel);
skel=imfill(skel,4,'holes');
%StartPoint = 1st endpoint

distImg=bwdist(bw);
maxBblNum=6;
%[subMatrix labelNum skelImg bubbles tips lbbImg lbbLen lbbSubs]=decomposeSkel(bw,[68 116],1,10,distImg,maxBblNum)

skel=div_skeleton_new(4,1,1-bw,4);
skel=parsiSkel(skel);
skel=imfill(skel,4,'holes');
distImg=bwdist(~bw);


imshow(bw);
hold on
plot(bb, aa, '.r');
plot(y1, x1, 'og');
plot(y, x, '.g');