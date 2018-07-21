clear all;close all;clc
[filename,filepath]=uigetfile('.jpg','输入一个图像');% read a picture
file=strcat(filepath,filename); 
IMG=imread(file);
IMG = double(rgb2gray(IMG));% Convert to gray image
figure('name','Gray image'),imshow(mat2gray(IMG));
m_num = 10;
n_num = 30;
I = imresize(IMG,[m_num,n_num]);% resize to be easy to verify. Also can be [30 10]
[m,n] = size(I);

 % X*Q make every column of X move to the right circularly 
 Q = [zeros(n-1,1),eye(n-1);1,zeros(1,n-1)];
 
 % P*X make every row of X move down 
 P = [zeros(1,m-1),1;eye(m-1),zeros(m-1,1)];
 
 % Initialize a four-dimensional tensor. Concretely,
 % the first and second dimension show the picture,
 % the third dimension repensents the picture moved to the right ordinal
 % times,
 % and the fourth dimension repensents the picture moved down ordinal times
 tensor = zeros(m,n,n,m);
for i = 1:n
    for j = 1:m
        tempQ = Q^mod(i-1,n);
        tempP = P^mod(j-1,m);
        tensor(:,:,i,j) = tempP*I*tempQ;
    end
end

 % Convert tensor to 2-D form that can be observed.
tt = zeros(m*n,n*m);
for i = 1:n
    for j = 1:m
        tt((i-1)*m+1:i*m,(j-1)*n+1:j*n) = tensor(:,:,i,j);
    end
end
figure('name','Cyclic matrix')
imshow(mat2gray(tt));
 % 'tt' is the partitioned matrix of the n row m column, in which
 % each bolck is m*n. So 'tt' is square matrix.

 % Obtain the four-dimensional Fourier transform of the tensor
T = fftn(tensor);

 % 'TT' is the observed form of 'T'
TT = zeros(m*n,n*m);
for i = 1:n
    for j = 1:m
        TT((i-1)*m+1:i*m,(j-1)*n+1:j*n) = T(:,:,i,j);
    end
end
figure('name','Fourier transform'),imshow(log(abs(TT)));

 % The observed form of the result of Fourier transform
 % is no diagonalization matrix.
 % Fortunately,the position of the nonzero elements is full-rank, and the 
 % num of the elements is same with the rank. 
 % If we relpace the nonzero elements with '1' and others '0', the matrix
 % is orthogonal matrix(正交矩阵). 
 % The product of it and its transpose is diagonalization matrix.
[r,c] = find(log(abs(TT))>2);
waoo = zeros(m*n,n*m);% waoo is just an onomatopoeia word(哇哦..)
figure('name','Green * is the position of the nonzeros elements');
for i = 1:size(r,1)
plot(r(i),c(i),'g*','markersize',12);hold on;
waoo(r(i),c(i)) = 1;
end
figure('name','Product of TT and the transpose of TT')
imshow(abs(waoo*waoo'))
DT = TT*waoo';


% The follows is the process in which only 2-D Fourier transform is used
% The range of the 2-D Fourier transform is the m*n block.
I = imresize(IMG,[m_num,n_num]);
[m,n] = size(I);
Q = [zeros(n-1,1),eye(n-1);1,zeros(1,n-1)];
P = [zeros(1,m-1),1;eye(m-1),zeros(m-1,1)];
B = zeros(m*n,n*m);
 for i = 1:n
     for j = 1:m
         tempQ = Q^mod(i-1,n);
         tempP = P^mod(j-1,m);
         B(m*(i-1)+1:m*i,n*(j-1)+1:n*j) = fft2(tempP*I*tempQ);
         XI(m*(i-1)+1:m*i,n*(j-1)+1:n*j) = tempP*I*tempQ;
     end
 end
 
 % The matrix C is same(or like) the matrix TT.
 C = zeros(m*n,n*m);
 for i = 1:m
     for j = 1:n
         Ctemp = zeros(n,m);
         for ii = 1:n
             for jj = 1:m
                 Ctemp(ii,jj) = B((ii-1)*m+i,(jj-1)*n+j);
             end
         end
         C((i-1)*n+1:i*n,(j-1)*m+1:j*m) = fft2(Ctemp);
         CI((i-1)*n+1:i*n,(j-1)*m+1:j*m) = Ctemp;
     end
 end
 figure
 subplot(1,2,1)
 imshow(mat2gray(log(abs(B))+10*ones(size(B))));subplot(1,2,2)
 D = log(abs(C));
 % remove the '-inf'
 inf_ind = isinf(D);
 [r,c] = find(inf_ind == 1);
 for i = 1:size(r,1)
 D(r(i),c(i)) = -20;
 end
D = D + 30*ones(size(D));
imshow(mat2gray(D));
[r c] = find(abs(C)>2);
waoo = zeros(m*n,n*m);
figure;
for i = 1:size(r,1)
plot(r(i),c(i),'g*','markersize',20);hold on;
waoo(r(i),c(i)) = 1;
end
figure
imshow(abs(waoo*waoo'))
DC = C*waoo';

 % Operate the original image with 2-D Fourier transform.
 % Find if the value of 2-D transform is same with diagonal element of 
 % the 4-D Fourier transform of the tensor.
 If = m*n*fft2(I);
 samenum = 0;
 d2 = diag(DC);
 d1 = diag(DT);
 flag=0;
for i = 1:size(If,1)
    for j = 1:size(If,2)
    flag = find(abs(d2-If(i,j))<abs(0.0001*If(i,j)));% d2 can be replaced with d1
    if(flag)
        samenum = samenum + 1;% the number of the same elements
        If(i,j) = -inf;d2(flag(1)) = inf;% d2 or d1 should be corresponding to line 138
    end
    end
end
 
% 1、思路的转变
% 转化一下实验的目的，不直接考虑4-D情形（实在是难以想象，也没有实在的分析工具）
% 既然，对循环矩阵做4-D傅里叶变换确实可以得到与2-D变换一一对应的值。
% 不如从回归公式（1-D原始信号的情形）出发，
% 尝试得到对角矩阵，或者非零元素两两不同行不同列的矩阵。
% 
% 2、原始信号一维形式 及 一维形式的2-D傅里叶变换
% 2.1
% 将图像转化为一维,即将矩阵改为行排式。具体地,
% Img = [1 2 3;
%        4 5 6];
% [m,n] = size(Img);
% vector1 = reshape(Img',1,6) = [1 2 3 4 5 6]
% 2.2
% 能不能对行排式进行矩阵乘法等价于对原矩阵进行2-D傅里叶变换呢？
% 简单可得，
%     reshape(F1*Img*F2') = vector1*kron(F1.',eye(n))*kron(eye(m),F2');
% 可以通过Fouriertest.m 验证
% 
% 3、回归公式
% 观察回归公式，我们还需要将行排式再排成矩阵，不妨，先排右移n次的n个行排式，
% 再排下移一次，右移n次的n个行排式，得到m*n个行排式组成的矩阵，也是方阵。
% 而且，每一列其实是4-D张量另外两个维度的行排式，而且同样是先排序号为1:n的那一维。
% vector = [1 2 3 4 5 6;
%           3 1 2 6 4 5;
%           2 3 1 5 6 4;
%           4 5 6 1 2 3;
%           6 4 5 3 1 2;
%           5 6 4 2 3 1];
% 再说一遍，vector的行向量是图像（即横纵）两个维度的行排式
%           vector的列向量是下移m次右移n次两个维度的行排式
% 即vector的第一列的元素是平移后的图像的（1,1）处的像素值
% [1 4 3 6 2 5]'分别是 
% 不平移 右移1次 右移2次 下移1次 右移1次下移1次 右移2次下移1次
%
% 通过一维实验，得到F*C(x)*F'(共轭转置)是对角的，但是这并不是2-D Fourier transform
% 2-D is F*C(x)*F.'(转置)
% 结合二维实验，确实有
%     原图像2-D傅里叶变换的结果与循环矩阵4-D傅里叶变换的结果的值是一一对应的。
% 所以，得到一个猜想：不做4-D傅里叶变换而是将F*C(x)*F'推广到4-D情形。
% 加之对行排式做乘法可以得到2-D傅里叶变换（这里使用共轭转置形式也可以...
% 以下两个公式是实验的核心公式：
% TF = vector*FF;将每一行当作一个行排式，做2-D傅里叶变换
% TTF = (TF'*FF)' = FF'*vector*FF;将结果的每一列当作一个行排式，做2-D傅里叶对角化
%                    |
%               此处必须为共轭转置
% 有趣的是，这个猜想是正确的。
% 可以将回归公式中所用的大矩阵X对角化。
% 而且，对角矩阵的对角元素是原图像二维傅里叶变换结果的行排式，可以看
%  'mmd' 和 'If' 来证明。
% 具体实现如下
% 如果你想方面地修改验证实验的矩阵规模,去掉下边的注释符号'%'
% The experiment is implemented as follows.
% If you want to modify the scale of experimental matrix, 
% remove the the symbol '%'. 
% The English version of the above discusion is following the experiment.
% m_num = 30;n_num = 10;
I = imresize(IMG,[m_num,n_num]);
[m,n] = size(I);

% X*Q make every column of X move to the right circularly 
Q = [zeros(n-1,1),eye(n-1);1,zeros(1,n-1)];
 
% P*X make every row of X move down 
P = [zeros(1,m-1),1;eye(m-1),zeros(m-1,1)];
for i = 1:n
     for j = 1:m
         tempQ = Q^mod(i-1,n);
         tempP = P^mod(j-1,m);
         XI(m*(i-1)+1:m*i,n*(j-1)+1:n*j) = tempP*I*tempQ;
     end
 end
vector = zeros(m*n,m*n);
for i = 1:n
    for j = 1:m
        vector(i+(j-1)*n,:) = reshape(XI(m*(i-1)+1:m*i,n*(j-1)+1:n*j).',1,m*n);
    end
end
fly = n;
 Fn = zeros(fly,fly);
 for i = 1:fly
     for j = 1:fly
         Fn(i,j) = exp(-2*pi*1i*(i-1)*(j-1)/n);
     end
 end
 fly = m;
 Fm = zeros(fly,fly);
 for i = 1:fly
     for j = 1:fly
         Fm(i,j) = exp(-2*pi*1i*(i-1)*(j-1)/m);
     end
 end
FF = kron(Fm.',eye(n))*kron(eye(m),Fn.');
% FF2 = kron(Fn.',eye(m))*kron(eye(n),Fm');
TF = vector*FF;
TTF = FF'*vector*FF;
figure('name','TF'),imshow(mat2gray(log(abs(TF))))
figure('name','The row of TF is row type of this matrix'),imshow(mat2gray(log(abs(CI))))
figure('name','TTF'),imshow(mat2gray(((log(abs(TTF)))>0).*(log(abs(TTF))<200)))
mmd = diag(TTF);
If = m*n*fft2(I);
fprintf('\nThe maximum of the difference of the modulus of corresponding elements is')
max(abs(mmd - reshape(If.',m*n,1)))

% ***********************************************************************
% *******************     1.Change of Thought      **********************
% ***********************************************************************
% The purpose of the experiment was directly considering the 4-D case.
% (it is hard to imagine and to analysis without Mathematical Tools)
%
% Since the result of the 4-D Fourier transform of the circular matrix
% (seriously,4-D tensor) corresponds to 2-D transform one by one,
% The purpose now:
% From the regression formula, try to get diagonal matrix 
% or the matrix, of which the nonzero elements are neither the same row nor the same column. 
% (From the case of 1-D original signal).
%
% ***********************************************************************
% **  2.1-D original signal and 1-D form of 2-D Fourier transform.  *****
% ***********************************************************************
% 2.1
% Converting the image into one dimension is 
% changing the matrix into row type. Specifically,
% Img = [1 2 3;
%        4 5 6];
% [m,n] = size(Img);
% vector1 = reshape(Img',1,m*n) = [1 2 3 4 5 6]
% 2.2
% Is there a matrix multiplication for the row type,
% equivalent to the 2-D Fourier transform of the original matrix ?
% Easy to get,
%     reshape(F1*Img*F2') = vector1*kron(F1.',eye(n))*kron(eye(m),F2');
% where reshape is MATLAB command
% It can be verified by Fouriertest.m
%
% ***********************************************************************
% *****************     3. regression formula     ***********************
% ***********************************************************************
% Observing the regression formula, we also need arrange row type into a matrix. 
% First, the n row type with 1:n shift to the right is arranged.
% Second, the n row type with 1:n shift to the right and 1 downwards shift 
% and so on. 
% Then, a mn-by-mn matrix is obtained.
% Moreover, each column is actually a row type of two other dimensions of 
% 4-D tensor, and the dimension with 1:n ordinal first sorted for each column.
% vector = [1 2 3 4 5 6;
%           3 1 2 6 4 5;
%           2 3 1 5 6 4;
%           4 5 6 1 2 3;
%           6 4 5 3 1 2;
%           5 6 4 2 3 1];
% Again, the row vectors of 'vector' are row type of two dimensions of image 
% (i.e. horizontal and vertical coordinates).
% The column vectors of 'vector' are row type of two dimensions.
% (i.e. moved down 1:m times, right shifts to 1:n times,)
% The first element of vector is the pixel value at the (1,1) of the 
% translational image. 
% The first colunm [1 4 3 6 2 5]' are respectively
% the (1,1) of the translational image with no translation, 1 right, 2 right, 
% 1 downwards, 1 downwards & 1 right, and 1 downwards & 2 right translation.
% 
% By 1-D experiments, F*C(x)*F'(conjugate transposition) is 
% diagonal, but this is not 2-D Fourier transform:
%                           2-D is F*C (x) *F.'(transpose)
% Combined with the 2-D experiment, it is true that the result of 2-D 
% Fourier transform of the original image is one-to-one corresponding
% with the diagonal element value of 4-D Fourier transform of the tensor.
% Therefore, we get a conjecture: 
% If instead of doing 4-D Fourier transform, we extend F*C(x)*F'to 4-D case,
% the matrix 'vector' can be diagonalized.
%
% Considering the multiplication of row type is equivalent to 2-D Fourier 
% transform, the following is the key formula of the experiment.
% TF = vector*FF; Take every row as a row type, do 2-D Fourier transform.
% TTF = (TF'*FF) = FF'*vector*FF; Every column of the result is regarded 
%                   |
%                 must be conjugated transposed here
% as a row type is operated with 2-D Fourier transform.
% Interestingly, this conjecture is true.
% And the diagonal elements is the row type of the 2-D transform of the
% orginal image (a m-by-n matrix). 'mmd' and 'If' can verify this.


 