% 
% cd 'C:/Users/20203/ECO/sequences'
% obj = VideoReader('C:\\Users\\20203\\Documents\\Tencent Files\\2020354644\\FileRecv\\lab480p.mp4');%������Ƶλ��
% numFrames = obj.NumberOfFrames;% ֡������
% % weishu = floor(log10(numFrames))+1
% %  for k = 1 : numFrames
% %      frame = read(obj,k);%��ȡ�ڼ�֡
% %     % imshow(frame);%��ʾ֡
% %       imwrite(frame,strcat('lab/img/',num2str(k,'%04d'),'.jpg'),'jpg');% ����֡
% %  end
%  
% %  frame = read(obj, 1);
% %  imshow(frame);
% %  h = imrect
% M= repmat([198 21.0000000000001 57 216], numFrames, 1);
% dlmwrite('lab/groundtruth_rect.txt',M,'delimiter',' ')
% cd 'C:/Users/20203/ECO/sequences'
% img_files = num2str((1:3789)',  ['lab/img/','%04i.jpg']);
% for i=1:3789
%     big = imread(img_files(i,:));
%     small=imresize(big,[270,360]);
%     imwrite(small, strcat('lab_small/img/',num2str(i,'%04d'),'.jpg'),'jpg');
% end

% cd 'C:/Users/20203/ECO/sequences'
% a = imread(strcat('lab_small/img/',num2str(1,'%04d'),'.jpg')); 
% imshow(a);
% h=imrect;

vid = videoinput('winvideo');
set(vid,'ReturnedColorSpace','rgb');
vidRes=get(vid,'VideoResolution');
width=vidRes(1);
height=vidRes(2);
nBands=get(vid,'NumberOfBands');
figure('Name', 'Matlab��������ͷ By Lyqmath', 'NumberTitle', 'Off', 'ToolBar', 'None', 'MenuBar', 'None');
hImage=image(zeros(vidRes(2),vidRes(1),nBands));
preview(vid,hImage);