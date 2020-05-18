%% 
I=imread('cell.png'); %读取图片
I_light=double(I/255);  %转化为double类型并将矩阵值归一化

%RGB三路值
I_r = I_light(:,:,1);
I_g = I_light(:,:,2);
I_b = I_light(:,:,3);
[row,col]=size(I_r);%图片大小
%% 查看RGB分量图便于分析从哪个分量进行分析
% figure
% subplot(221),imshow(I_r),title('R')
% subplot(222),imshow(I_g),title('G')
% subplot(223),imshow(I_b),title('B')
% subplot(224),imshow(I),title('原图')
%% 阈值法分离红色细胞
I_r_t=I_r;

 for i=1:row
     for j=1:col
         if I_r(i,j)<0.445
             I_r_t(i,j)=0;
         else
             I_r_t(i,j)=255;
         end
     end
 end
 I_r_t=uint8(I_r_t);
%  imshow(I_r_t),title('阈值处理后红色细胞区域');
%开运算(去除小点,开运算不影响面积，先三次腐蚀再三次膨胀）
se=strel('disk',10);
I_ero=I_r_t;
i1=0;i2=0;
while i1<3
    I_ero=imerode(I_ero,se);
    i1=i1+1;
end
[L_1,NUM_1]=bwlabel(I_ero,8); %多次腐蚀后进行标记计数
I_ero_dilate=I_ero;
while i2<3
    I_ero_dilate=imdilate(I_ero_dilate,se);
    i2=i2+1;
end
imshow(I_ero_dilate)
%计算红色细胞面积
count_1=0;
for i=1:row
     for j=1:col
         if I_ero_dilate(i,j)==255
             count_1=count_1+1;
         end
     end
 end

%% 第二种细胞
%先进行阈值分离
I_g_t=I_g;
 for i=1:row
     for j=1:col
         if I_g(i,j)<0.55
             I_g_t(i,j)=0;
         else
             I_g_t(i,j)=255;
         end
     end
 end
 figure, imshow(I_g_t)
I_g_t_use=I_g_t;
 %因为粘连严重，使用最终腐蚀和条件腐蚀
m=row;n=col;
imgn=zeros(m,n);
preimg=imgn;

se= strel('square',3);
while sum(sum(preimg-I_g_t_use))~=0 
    preimg=I_g_t_use;
    I_g_t_use=I_g_t_use>0;
   
    [I_g_t_use, NUM_2]=bwlabel(I_g_t_use,8);      %标记不同区域，label是区域个数
    L_2=I_g_t_use;
    imgn=imerode(I_g_t_use,se);
    
    %腐蚀之后是否有哪个被标记的区域消失了
    Hist=zeros(1,NUM_2);            
    for i=1:m
        for j=1:n
            if imgn(i,j)~=0
                Hist(imgn(i,j))=imgn(i,j);  
            end
        end
    end

    %统计消失区域的标号
    H=[];
    for i=1:NUM_2
        if Hist(i)==0
            H=[H i];       
        end
    end
    
    %如果这个区域消失了，那么再把这个区域恢复过来
    if ~isempty(H)
        l=length(H);
        for i=1:m
            for j=1:n
                for k=1:l
                    if I_g_t_use(i,j)==H(k)   
                        imgn(i,j)=I_g_t_use(i,j);
                    end
                end
            end
        end   
    end
           
    I_g_t_use=imgn;
end

figure;
imshow(imdilate(I_g_t_use,se));    %再膨胀一下好看

%计算面积
count_2=0;
for i=1:row
     for j=1:col
         if I_g_t(i,j)==255
             count_2=count_2+1;
         end
     end
end
 
%% 原图除去已求细胞1和2的区域
area_sovled=(uint8(I_g_t)+I_ero_dilate);
for i=1:row
     for j=1:col
         if area_sovled(i,j)==255
             area_sovled(i,j)=0;
         else
             area_sovled(i,j)=1;
         end
     end
end
area_sovled=cat(3,area_sovled,area_sovled,area_sovled);
area_after=I.*area_sovled; 

imshow(area_after),title('除去已经处理的细胞1和2，剩下的区域')

%% 

area_adjust=area_after(:,:,1);
area_adjust1=area_after(:,:,2);
area_adjust2=area_after(:,:,3);
figure,imshow(area_adjust),title('adjust')

%% 获得4区域的mask
mask=zeros(m,n);
%利用hsv空间进行区域分割
[h,~,~]=rgb2hsv(I);
for i=1:row
     for j=1:col
         if h(i,j)<0.65&&h(i,j)>0.55
             mask(i,j)=1;
         else
             mask(i,j)=0;
         end
     end
end
%进行多次的膨胀腐蚀
mask=imdilate(mask,se);
mask=imdilate(mask,se);
mask=imdilate(mask,se);
mask=imdilate(mask,se);
mask=imdilate(mask,se);
mask=imdilate(mask,se);
mask=imerode(mask,se);
mask=imerode(mask,se);
mask=imerode(mask,se);
mask=imerode(mask,se);
mask=imerode(mask,se);
mask=imerode(mask,se);
mask=imerode(mask,se);
mask=imerode(mask,se);
mask=imerode(mask,se);
mask=imerode(mask,se);

imshow(mask),title('细胞4区域的mask')
%%
area_adjust_mask=area_adjust.*uint8(mask);
area_adjust_mask1=area_adjust1.*uint8(mask);
area_adjust_mask2=area_adjust2.*uint8(mask);
a=cat(3,area_adjust_mask,area_adjust_mask1,area_adjust_mask2);
imshow(area_adjust_mask);
area_adjust_mask=imadjust(area_adjust_mask);
figure,imshow(area_adjust_mask)
% a=gamma_light(a);
% figure,imshow(a),title('细胞4割分区域在原图上的展示')
% bk=ordfilt2(area_adjust_mask,1,ones(3,3),'symmetric');
% h=ones(3,3)/(9);
% bk=imfilter(bk,h,'replicate');
% area_adjust_mask=imsubtract(double(area_adjust_mask),double(bk));
% figure,imshow(area_adjust_mask),title('亮度调整')
% %% maybe nouse
% se1= strel('disk',2);
% area_adjust1=imerode(area_adjust_mask,se1);
% area_adjust1=imerode(area_adjust1,se1);
% area_adjust1=imdilate(area_adjust1,se1);
% 
% imshow(area_adjust1)
% bk=ordfilt2(area_adjust1,1,ones(3,3),'symmetric');
% h=ones(3,3)/(9);
% bk=imfilter(bk,h,'replicate');
% area_adjust2=imsubtract(double(area_adjust1),double(bk));
% for i=1:row
%      for j=1:col
%          if mask(i,j)==0
%              area_adjust2(i,j)=255;
%          end
%      end
% end
% figure,imshow(area_adjust2),title('亮度调整')
% 
% %% maybe nouse
% % for i=1:row
% %      for j=1:col
% %          if area_adjust_mask(i,j)>0&&area_adjust_mask(i,j)<0.01
% %              area_adjust_mask(i,j)=1;
% %          else
% %              area_adjust_mask(i,j)=0;
% %          end
% %      end
% % end
% % area_mask=area_adjust_mask;
% % imshow(area_mask)
% for i=1:row
%      for j=1:col
%          if area_adjust2(i,j)==0
%              area_adjust2(i,j)=255;
%          else
%              area_adjust2(i,j)=0;
%          end
%      end
% end
% imshow(area_adjust2)
% area_adjust3=area_adjust2;
% %%
% 
% se1= strel('disk',2);
% 
% area_adjust3=imerode(area_adjust2,se1);
% 
% % area_adjust3=imdilate(area_adjust3,se1);
% % area_adjust3=imdilate(area_adjust3,se1);
% 
% % [p,r,~]=imfindcircles(area_adjust,[1,40]);
% 
% % [centers_r, radii_r,~]=imfindcircles(area_adjust3,[10,35]);
% 
% imshow(area_adjust3);
%%  分水岭法的部分操作
I_gray=area_adjust_mask;
% gmag = imgradient(I_gray);
% imshow(gmag,[]),title('Gradient Magnitude')
se=strel('disk',13);
Io = imopen(I_gray,se);
imshow(Io),title('Opening')

Ie = imerode(I_gray,se);
Iobr = imreconstruct(Ie,I_gray);
imshow(Iobr)
title('Opening-by-Reconstruction')


Ioc = imclose(Io,se);
imshow(Ioc)
title('Opening-Closing')

Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
imshow(Iobrcbr)
title('Opening-Closing by Reconstruction')

fgm = imregionalmax(Iobrcbr);
% imshow(fgm)
% title('Regional Maxima of Opening-Closing by Reconstruction')

% I2 = labeloverlay(I_gray,fgm);
% imshow(I2)
% title('Regional Maxima Superimposed on Original Image')

se2 = strel(ones(5,5));
fgm2 = imclose(fgm,se2);
fgm3 = imerode(fgm2,se2);
figure,imshow(fgm3),title('需检测连通数量')
fgm4 = bwareaopen(fgm3,20);
I3 = labeloverlay(I_gray,fgm4);
% figure,imshow(I3)
% title('Modified Regional Maxima Superimposed on Original Image')

[L_4,NUM_4]=bwlabel(fgm3,4);
%计算部分4面积
count_4=0;
for i=1:row
     for j=1:col
         if mask(i,j)==1
             count_4=count_4+1;
         end
     end
end

%% 处理第三部分细胞（照样用到分水岭的部分操作）
area_12=area_sovled(:,:,1);
mask_123=imcomplement(mask);%二值图反像素
area_124=double(area_12).*mask_123; %area_124表示124区的细胞区域为0
% imshow(area_124)
area_124rgb=cat(3,area_124,area_124,area_124);
area_3=double(I)/255.*area_124rgb;  %原图乘以矩阵,将3区分离出来
figure,imshow(area_3)   %第三区细胞区域

area_3_gray=rgb2gray(area_3);

gmag = imgradient(area_3_gray);
% imshow(gmag,[]),title('Gradient Magnitude')
se=strel('disk',4);
Io = imopen(area_3_gray,se);
% imshow(Io),title('Opening')

Ie = imerode(area_3_gray,se);
Iobr = imreconstruct(Ie,area_3_gray);
% imshow(Iobr)
% title('Opening-by-Reconstruction')


Ioc = imclose(Io,se);
% imshow(Ioc)
% title('Opening-Closing')

Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
% imshow(Iobrcbr)
% title('Opening-Closing by Reconstruction')

fgm = imregionalmax(Iobrcbr);
figure,imshow(fgm);
title('Regional Maxima of Opening-Closing by Reconstruction')

I2 = labeloverlay(area_3_gray,fgm);
imshow(I2)
title('Regional Maxima Superimposed on Original Image')

se123=strel('disk',2);
fgm1=imdilate(fgm,se123);
[L_3,NUM_3]=bwlabel(fgm1,8);
figure,imshow(fgm)
%计算部分4面积
count_3=0;
for i=1:row
     for j=1:col
         if area_124(i,j)==1
             count_3=count_3+1;
         end
     end
end

%% 计算面积及输出
area_the_picture=row*col;
area_cell_1=count_1/(area_the_picture);
area_cell_2=count_2/(area_the_picture);
area_cell_3=count_3/(area_the_picture);
area_cell_4=count_4/(area_the_picture);
cell_1=strcat(['细胞1的数量为：',num2str(NUM_1),', 面积占比:',num2str(area_cell_1)]);
cell_2=strcat(['细胞2的数量为：',num2str(NUM_2),', 面积占比:',num2str(area_cell_2)]);
cell_3=strcat(['细胞3的数量为：',num2str(NUM_3),', 面积占比:',num2str(area_cell_3)]);
cell_4=strcat(['细胞4的数量为：',num2str(NUM_4),', 面积占比:',num2str(area_cell_4)]);
disp(cell_1);
disp(cell_2);
disp(cell_3);
disp(cell_4);
