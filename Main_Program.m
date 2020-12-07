function varargout = Main_Program(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Main_Program_OpeningFcn, ...
                   'gui_OutputFcn',  @Main_Program_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end



function Main_Program_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;
set(handles.axes1,'Visible','off');
set(handles.axes2,'Visible','off');

guidata(hObject, handles);



function varargout = Main_Program_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;


% Pushbutton1 to appy Wiener Restored Image.
function pushbutton1_Callback(hObject, eventdata, handles)

% User Select image from folder 
[baseFileName, folder] = uigetfile({'*.jpg';'*.png';'*.bmp'}); 
fullFileName = fullfile(folder, baseFileName);
InputImage = imread(fullFileName);
% Write the image
imwrite(InputImage,'InputImage.png');
axes(handles.axes1);
imshow(InputImage);
title('Original Image');
I = double(InputImage);
OIm = I;
window_size = 3;
means = imfilter(I, fspecial('average', window_size), 'replicate');
sigmas = sqrt((I-means).^2/window_size^2);
sigmas = imfilter(sigmas, fspecial('average', window_size), 'replicate');

ENLs = (means./sigmas).^2;
sx2s = ((ENLs.*(sigmas).^2) - means.^2)./(ENLs + 1);
fbar = means + (sx2s.*(I-means)./(sx2s + (means.^2 ./ENLs)));
OIm(means~=0) = fbar(means~=0);

axes(handles.axes2);
imshow(uint8(OIm));
title('Wiener Filtered Image')
% Write the image
imwrite(OIm,'OutPutImage.png');
%% READ IMAGE FOR ANALYSIS
img= imread('InputImage.png');
cover_object1= imread('OutPutImage.png');

% PSNR, MSE, BER, SSIM, Entropy, and EME for Wiener applyed here.
psnr1=psnr((img),(cover_object1));
set(handles.edit1,'string',psnr1);
ENtrop = entropy(img);
set(handles.edit2,'string',ENtrop);
MSEE = immse(cover_object1,img);
set(handles.edit3,'string',MSEE);
ssimval = ssim(cover_object1,img);
set(handles.edit4,'string',ssimval);
BERR = Biter(cover_object1,img);
set(handles.edit5,'string',BERR);
err = msecode(img,cover_object1);
set(handles.edit6,'string',err);

% Pushbutton2 for Adapt Blind Deconvolution Restored Image.
function pushbutton2_Callback(hObject, eventdata, handles)

% User Select image from folder 
[baseFileName, folder] = uigetfile({'*.jpg';'*.png';'*.bmp'}); 
fullFileName = fullfile(folder, baseFileName);
I = imread(fullFileName);
axes(handles.axes1);
imshow(I)
title('Original Image')
% Write the image
imwrite(I,'InputImage.png');
PSF = fspecial('motion',2,3);
Blurred = imfilter(I,PSF,'circ','conv');
INITPSF = ones(size(PSF));
[J P] = deconvblind(Blurred,INITPSF,30);
axes(handles.axes2);
imshow(J)
title('Adapt Blind Deconvolution Restored Image')
% Write the image
imwrite(J,'OutPutImage.png');
%% READ IMAGE FOR ANALYSIS
img= imread('InputImage.png');
cover_object1= imread('OutPutImage.png');

% PSNR, MSE, BER, SSIM, Entropy, and EME for Adapt Blind Deconvolution Restored applyed here.
psnr1=psnr(img,cover_object1);
set(handles.edit1,'string',psnr1);
ENtrop = entropy(img);
set(handles.edit2,'string',ENtrop);
MSEE = immse(cover_object1,img);
set(handles.edit3,'string',MSEE);
ssimval = ssim(cover_object1,img);
set(handles.edit4,'string',ssimval);
BERR = Biter(cover_object1,img);
set(handles.edit5,'string',BERR);
err = msecode(img,cover_object1);
set(handles.edit6,'string',err)

% pushbutton3 is for Lucy-Richardson Restored Image.
function pushbutton3_Callback(hObject, eventdata, handles)

% User Select image from folder 
[baseFileName, folder] = uigetfile({'*.jpg';'*.png';'*.bmp'}); 
fullFileName = fullfile(folder, baseFileName);
I = imread(fullFileName);

axes(handles.axes1);
imshow(I)
title('Original Image')
imwrite(I,'InputImage.png');
PSF = fspecial('gaussian',5,5);
V = 0.002;
luc1 = deconvlucy(uint8(I),PSF,5);
axes(handles.axes2);
imshow(luc1)
title(' Lucy-Richardson Restored Image')

imwrite(luc1,'OutPutImage.png');
%% READ IMAGE FOR ANALYSIS
img= imread('InputImage.png');
cover_object1= imread('OutPutImage.png');

% PSNR, MSE, BER, SSIM, Entropy, and EME for Lucy-Richardson Restored  applyed here.
psnr1=psnr((img),(cover_object1));
set(handles.edit1,'string',psnr1);
ENtrop = entropy(img);
set(handles.edit2,'string',ENtrop);
MSEE = immse(cover_object1,img);
set(handles.edit3,'string',MSEE);
ssimval = ssim(cover_object1,img);
set(handles.edit4,'string',ssimval);
BERR = Biter(cover_object1,img);
set(handles.edit5,'string',BERR);
err = msecode(img,cover_object1);
set(handles.edit6,'string',err)


% Pushbutton4 is the Proposed Algorithm.
function pushbutton4_Callback(hObject, eventdata, handles)


%% % User  load image from folder image
[baseFileName, folder] = uigetfile({'*.jpg';'*.png';'*.bmp'}); 
fullFileName = fullfile(folder, baseFileName);
inputImage = imread(fullFileName);

axes(handles.axes1);
imshow(inputImage);title('Input image')
imwrite(inputImage,'InputImage.png');
pause(0.01)
%% deconvolution
inputImage=im2double(inputImage);
ground_truth = [40 40]; 
[u h] = ProposedAlgorithm(inputImage, ground_truth); 

%% display result
axes(handles.axes2);
imshow(u);title('Proposed Algorithm Output image')

imwrite(u,'OutPutImage.png');

%%

envCfg = coder.gpuEnvConfig('host');
envCfg.DeepLibTarget = 'cudnn';
envCfg.DeepCodegen = 1;
envCfg.Quiet = 1;

%% |mobilenetv2_predict|  is applyed here
net = mobilenetv2();


%%
type('mobilenetv2_predict.m')

%% Run MEX Code Generation
% the MobileNet-v2 network.
cfg = coder.gpuConfig('mex');
cfg.TargetLang = 'C++';
cfg.DeepLearningConfig = coder.DeepLearningConfig('cudnn');
cover_object1= imread('OutPutImage.png');
cover_object=imresize(cover_object1,[224 224]);
        predict_scores = mobilenetv2_predict(double(cover_object));
%% 
% Classify the data
[scores,indx] = sort(predict_scores, 'descend');
classNames = net.Layers(end).ClassNames;
classNamesTop = classNames(indx(1:5));
figure
barh(scores(5:-1:1),'magenta')
xlabel('Probability')
yticklabels(classNamesTop(5:-1:1))
title('Top Five Predictions By Using MobileNet-v2')



clear mex;

%%
%% READ IMAGE FOR ANALYSIS
img= imread('InputImage.png');

% PSNR, MSE, BER, SSIM, Entropy, and EME for proposed Algorithm applyed here.
psnr1=PSNR_RGB1(double(img),double(cover_object1));
set(handles.edit1,'string',psnr1);
ENtrop = entropy(img);
set(handles.edit2,'string',ENtrop);
MSEE = immse(cover_object1,img);
set(handles.edit3,'string',MSEE);
ssimval = ssim(cover_object1,img);
set(handles.edit4,'string',ssimval);
BERR = Biter(cover_object1,img);
set(handles.edit5,'string',BERR);
err = msecode(img,cover_object1);
set(handles.edit6,'string',err);


% Pushbutton5 is for Correlation Cofficient.
function pushbutton5_Callback(hObject, eventdata, handles)

% Read input image for corelation cofficient
img=imread('InputImage.png');

figure
subplot(221)
imshow('InputImage.png');title('Original image');
[M N]=size(img(:,:,3));
im=double(img(:,:,3)); % to see the cc of othetr channel
im=mod(im,256);
im1=im(1:M-1,1:N-1);
[row col]=size(im1);
X=[];
Y=[];
%% (X,Y) coordinate
for i=1:1000
    x=randi(row);
    X=[X;x];
    
    y=randi(col);
    Y=[Y;y];
end
A=[X Y];
B=[X Y+1];
pixelXY=[]; pixelXY1=[];
for i=1: length(X)
pixelXY=[pixelXY ;im(X(i),Y(i))];
%%(X,Y+1) coordinate


pixelXY1=[pixelXY1;im(X(i),Y(i)+1)];
end
subplot(222)
plot(pixelXY,pixelXY1,'r.');
title('Correlation of original adjacent pixels in Horizontal direction');
xlabel('Pixel location (x,y)');
ylabel('pixel location (x+1,y)');


L=1000;
s=0;

for i=1:L
    
     s=s+pixelXY(i);

end
Ex=s/L;
tic
sm=0;
for i=1:L
    
    sm=sm+(pixelXY(i)-Ex)^2;
end
Dx=sm/L;
toc
%for any adjacent pixels


s1=0;

for i=1:L
    
     s1=s1+pixelXY1(i);

end
Ey=s/L;
sm1=0;
for i=1:L
    
    sm1=sm1+(pixelXY1(i)-Ey)^2;
end
Dy=sm1/L;
convXY=[];
for i=1:L
    
    
xx=(pixelXY(i)-Ex)*(pixelXY1(i)-Ey);
convXY=[convXY;xx];
end
convXY=1/L*sum(convXY);
rxy=convXY/(sqrt(Dx)*sqrt(Dy));

%% Read Out put image for corellation cofficient
img= imread('OutPutImage.png');
subplot(223)
imshow(img);title('Restored Image');
[M N]=size(img(:,:,3));
im=double(img(:,:,3)); % to see the cc of othetr channel
im=mod(im,256);
im1=im(1:M-1,1:N-1);
[row col]=size(im1);
X=[];
Y=[];
%% (X,Y) coordinate
for i=1:1000
    x=randi(row);
    X=[X;x];
    
    y=randi(col);
    Y=[Y;y];
end
A=[X Y];
B=[X Y+1];
pixelXY=[]; pixelXY1=[];
for i=1: length(X)
pixelXY=[pixelXY ;im(X(i),Y(i))];
%%(X,Y+1) coordinate


pixelXY1=[pixelXY1;im(X(i),Y(i)+1)];
end
subplot(224)
plot(pixelXY,pixelXY1,'b.');
title('Correlation of Restored Image adjacent pixels in Horizontal direction');
xlabel('Pixel location (x,y)');
ylabel('pixel location (x+1,y)');


L=1000;
s=0;

for i=1:L
    
     s=s+pixelXY(i);

end
Ex=s/L;
tic
sm=0;
for i=1:L
    
    sm=sm+(pixelXY(i)-Ex)^2;
end
Dx=sm/L;
toc
%for any adjacent pixels


s1=0;

for i=1:L
    
     s1=s1+pixelXY1(i);

end
Ey=s/L;
sm1=0;
for i=1:L
    
    sm1=sm1+(pixelXY1(i)-Ey)^2;
end
Dy=sm1/L;
convXY=[];
for i=1:L
    
    
xx=(pixelXY(i)-Ex)*(pixelXY1(i)-Ey);
convXY=[convXY;xx];
end
convXY=1/L*sum(convXY);
rxy=convXY/(sqrt(Dx)*sqrt(Dy));

% pushbutton6 is  for Histrogram Analysis.
function pushbutton6_Callback(hObject, eventdata, handles)

img= imread('InputImage.png');
cover_object1= imread('OutPutImage.png');

figure
subplot(221)
imshow('InputImage.jpg');title('Original image');
subplot(222)
imhist(img);title('Original image Histrogram');

subplot(223)
imshow(cover_object1);title('Restored Image');
subplot(224)
imhist(cover_object1);title('Restored Image Histrogram');


% Pushbutton7 to restart the program.
function pushbutton7_Callback(hObject, eventdata, handles)

close all;
clc;
clear all;
Main_Program


%  pushbutton8 is to close the program.
function pushbutton8_Callback(hObject, eventdata, handles)

clc;
close all;
clear all;


function edit1_Callback(hObject, eventdata, handles)



% Edit1 is used to display the result
function edit1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)


% Edit2 is used to display the result
function edit2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)


% Edit3 is used to display the result
function edit3_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)


% Edit4 is used to display the result
function edit4_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)


% Edit5 is used to display the result
function edit5_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)


% Edit6 is used to display the result
function edit6_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
