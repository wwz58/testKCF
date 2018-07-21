clear all;close all;clc
[filename,filepath]=uigetfile('.jpg','输入一个图像');% read a picture
file=strcat(filepath,filename); 
IMG=imread(file);
IMG = double(rgb2gray(IMG));% Convert to gray image
figure('name','Gray image'),imshow(mat2gray(IMG));
I = imresize(IMG,[30,10]);% resize to be easy to verify. Also can be [30 10]
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


% The follows is the process in which only 2-D Fourier transform is used
% The range of the 2-D Fourier transform is the m*n block.
I = imresize(IMG,[10,30]);
[m,n] = size(I);
Q = [zeros(n-1,1),eye(n-1);1,zeros(1,n-1)];
P = [zeros(1,m-1),1;eye(m-1),zeros(m-1,1)];
B = zeros(m*n,n*m);
 for i = 1:n
     for j = 1:m
         tempQ = Q^mod(i-1,n);
         tempP = P^mod(j-1,m);
         B(m*(i-1)+1:m*i,n*(j-1)+1:n*j) = fft2(tempP*I*tempQ);
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