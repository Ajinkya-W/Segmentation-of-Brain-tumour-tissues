%% reading Image of .T1ce, .FLAIR, seg.nii(gold standard mask) for patient 
%directory location for T1ce and Flair
%We need to change directory for applying algoon different patients
dirt1ce ='D:\BML\project\Data_BraTS2018\Brats18_2013_2_1\Brats18_2013_2_1_t1ce.nii';
dirgtmask='D:\BML\project\Data_BraTS2018\Brats18_2013_2_1\Brats18_2013_2_1_seg.nii';
dirfl='D:\BML\project\Data_BraTS2018\Brats18_2013_2_1\Brats18_2013_2_1_flair.nii';

%reading those .nii images as niftyread() is not supported in 15a of matlab
t1c_nii = load_nii(dirt1ce);
%t1c is t1 contrast enhanced image
t1c = t1c_nii.img;
t1c = imrotate(t1c,90);

fl_nii = load_nii(dirfl);
%fl is FLAIR image
fl = fl_nii.img;
fl = imrotate(fl,90);

gtmask_nii = load_nii(dirgtmask);
%gtmask is goldstandard segmented tumor image
gtmask = gtmask_nii.img;
gtmask = imrotate(gtmask,90);

%showing one slice of each
figure, imagesc(t1c(:,:,106)),colormap(gray),axis image,axis off,title('106th Slice of t1ce');

figure, imagesc(fl(:,:,106)),colormap(gray),axis image,axis off,title('106th Slice of Flair');

figure, imagesc(gtmask(:,:,106)),colormap(gray),axis image,axis off,title('106th Slice of Ground truth mask');
%% Region Growing
%UNCOMMENT AFTER RUNNING NORMAL ALGORITHM
% U_nii =  load_nii(t1c_nii);
% U = U_nii.img;
% U = imrotate(U,90);
% %figure, imagesc(U(:,:,105)),colormap(jet),axis image,axis off,title('slice mask');
% [y, x] = ginput();
% % img1 = double(Img(:,:,1));
%    x1 = 109;%96 88 4% 
%    y1 = 92;
% x1=round(x);
% y1=round(y);
% thres = 0.06*(max(max(V(:,:,106)))- min(min(V(:,:,106))));
% [polygon, maskk] = rg(V,[x1,y1,106],thres, 130, false, true, false);
% % index 
% figure,imagesc(mask(:,:,106)),colormap(gray),axis image,axis off,title('Result of Region Growing');
% dice=getDiceCoeff(U,maskk);
% nii=make_nii(maskk);
% save_nii(nii,'orig.nii');
%% Algorithm Begins    
dand=0;dplus=0;
    
for i = 1:155
    T1 = t1c(:,:,i);
    Fl = fl(:,:,i);
    gImg =gtmask(:,:,i);
    gImgf =gtmask(:,:,i);
    gImg1=gtmask(:,:,i);
    gImge=gtmask(:,:,i);
    % Applying thresholding of non-zero no. of pixels in ground truth and also associated max. value  of ground truth
    if max(gImg)<4 | nnz(gImg(:))<1000
        continue;
    end
    
    %% Normalizing the image
    T1N = uint8((double(T1-min(T1(:)))./double(max(T1(:))-min(T1(:))))*255);
    FlN = uint8((double(Fl-min(Fl(:)))./double(max(Fl(:))-min(Fl(:))))*255);
   
    gImg1(gImg(:)==1|gImg(:)==2)=0;
    gImgf(gImg(:)==1|gImg(:)==4)=0;
    gImge(gImg(:)==2|gImg(:)==4)=0;
    W=240;H=240;
    Img_maskt1 = zeros(W,H);
    Img_maskfl = zeros(W,H);
    %% Applying k means
    k=4;k1=5;
    %k1 for getting different segmentation in F1N image’s tumour tissues
    
    [mut,maskt1]=kmeans1(T1N,k);
    [muf,maskfl]=kmeans1(FlN,k1);
    Img_maskt1(maskt1(:)==4)=1;
    
    %getting a third portion also by filling with region 4 and then subtracting it from original
    Img_mask0=imfill(Img_maskt1,'holes');
    Img_mask0=Img_maskt1-Img_mask0;
    Img_maskfl(maskfl(:)==5)=1;

    %% morphological process on mask
    se = strel('disk', 1); 
    sef = strel('square', 1); 
    maskt1 = imopen(Img_maskt1, sef);
    maskfl = imopen(Img_maskfl, se);
    maske1 = imopen(Img_mask0, se);
       
    %% dice coefficeint computation
    %adding up all tumor tissues of a brain
    mask=maske1+maskfl+maskt1;
    
    dand=dand+nnz(gImg & mask);
    dplus=dplus+nnz(gImg) + nnz(mask);    
end
 dice = 2*dand/dplus;
 fprintf('dice coefficient : %f\t',dice);
  
 %% getting segmented imagesof tumors
 segImgt1=zeros(240,240,155);
 segImgf=zeros(240,240,155);
 segImge=zeros(240,240,155);
 segImg=zeros(240,240,155);
 for i=1:155
    maskt2=maskt1;
    maskt2(maskt1(:)>0)=1;maskt2(maskt1(:)==0)=0;
    segImgt1(:,:,i)=T1N.*uint8(maskt2);
    
    
    maskf2=maskfl;
    maskf2(maskfl(:)>0)=1;maskf2(maskfl(:)==0)=0;
    segImgf(:,:,i)=T1N.*uint8(maskf2);
    
    maske2=maske1;
    maske2(maske1(:)>0)=1;maske2(maske1(:)==0)=0;
    segImge(:,:,i)=T1N.*uint8(maske2);
    
    mask1=mask;
    mask1(mask(:)>0)=1;mask1(mask(:)==0)=0;
    segImg(:,:,i)=T1N.*uint8(mask1);
 end
 %% histograms of all tumor tissues
 figure,hist(segImg(segImg>0),100),title('Full mask');
 figure,hist(segImgt1(segImgt1>0),100),title('t1ce mask');
 figure,hist(segImgf(segImgf(:)>0),100),title('flair');
 figure,hist(segImge(segImge>0),100),title('filled one'); 

 %% mean
 
 meanIntensityt1=mean(nonzeros(segImgt1(:)));
 meanIntensityf=mean(nonzeros(segImgf(:)));
 meanIntensitye1=mean(nonzeros(segImge(:)));
 meanIntensity=mean(nonzeros(segImg(:)));
 
 %% sd
 
 sdt1=std(segImgt1(segImgt1(:)>0));
 sdf=std(segImgf(segImgf(:)>0));
 sd=std(segImg(segImg(:)>0));
 sde=std(segImge(segImge(:)>0));
 
 %% volume
 
 volt1=nnz(segImgt1(:));
 volf=nnz(segImgf(:));
 vol=nnz(segImg(:));
 vole=nnz(segImgf(:));