function varargout = SelectROI(varargin)
% SELECTROI MATLAB code for SelectROI.fig
%      SELECTROI, by itself, creates a new SELECTROI or raises the existing
%      singleton*.
%
%      H = SELECTROI returns the handle to a new SELECTROI or the handle to
%      the existing singleton*.
%
%      SELECTROI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECTROI.M with the given input arguments.
%
%      SELECTROI('Property','Value',...) creates a new SELECTROI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SelectROI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SelectROI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SelectROI

% Last Modified by GUIDE v2.5 19-Jan-2017 14:43:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SelectROI_OpeningFcn, ...
                   'gui_OutputFcn',  @SelectROI_OutputFcn, ...
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
% End initialization code - DO NOT EDIT


% --- Executes just before SelectROI is made visible.
function SelectROI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SelectROI (see VARARGIN)

% Choose default command line output for SelectROI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SelectROI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

set(handles.axes1,'visible','off');
set(handles.axes2,'visible','off');
set(handles.axes3,'visible','off');

% --- Outputs from this function are returned to the command line.
function varargout = SelectROI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% ================= Part 1: Read DCM file ====================
% Open Select DCM file
pathName = uigetdir('','Please pick a Directory');
folderName = dir(pathName);
folderNumber = length(folderName);
totalHdrNumber = 0;
global hdrPathName;
for i = 1:folderNumber
    hdrName = dir([pathName,'\',folderName(i).name,'\','*.hdr']);
    hdrNumber = length(hdrName);
    for j = 1:hdrNumber
        totalHdrNumber = totalHdrNumber + 1;
        hdrPathName(totalHdrNumber).pathName = [pathName,'\',folderName(i).name,'\'];
        hdrPathName(totalHdrNumber).fileName = hdrName(j).name;
    end
end

firstHdr = spm_vol([hdrPathName(1).pathName,hdrPathName(1).fileName]); 
firstImg = spm_read_vols(firstHdr);
img = firstImg(:,:,1);

% img = imread('I1.bmp');
axes(handles.axes1);
set(handles.axes1,'visible','on');
imshow(img,'DisplayRange',[]);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% ================= Part 1: SELECT ROI ====================

global hdrPathName;

totalHdrNumber = length(hdrPathName);
for hdrIndex = 1:totalHdrNumber
    hdrData = spm_vol([hdrPathName(hdrIndex).pathName,hdrPathName(hdrIndex).fileName]); 
    hdrImgData = spm_read_vols(hdrData);
    
    [height,width,imageNumber]=size(hdrImgData);
    for imageIndex = 1:imageNumber
        img = hdrImgData(:,:,imageIndex);
        axes(handles.axes1);
        set(handles.axes1,'visible','on');
        imshow(img,'DisplayRange',[]);
        img = double(img);
        
        [dx,dy] = gradient(img);
        imgGradient = sqrt(dx.*dx+dy.*dy);
        
        sizeOfData = 0;
        [m n] = size(img);
        y_max = 0;x_max = 0;
        y_min = n+1;x_min = m+1;
        for i=1:m
          for j=1:n
              if(img(i,j) ~=0)
                  sizeOfData = sizeOfData+1;
                  if(i<x_min) x_min = i;end
                  if(i>x_max) x_max = i;end
                  if(j<y_min) y_min = j;end
                  if(j>y_max) y_max = j;end
                  imgData(1,sizeOfData) = img(i,j);
                  imgGradientData(1,sizeOfData) = imgGradient(i,j);
              end
           end
        end
        if(sizeOfData == 0) continue;end;
        imgROI = img(x_min:x_max,y_min:y_max);
        
        axes(handles.axes2);
        set(handles.axes2,'visible','on');
        imshow(imgROI,'DisplayRange',[]);
        
        % ����ֱ��ͼ
        [m,n] = size(imgData);
        his = zeros(1,256);
        for k = 0:255   his(k+1)=length(find(imgData == k))/(m*n);end
        axes(handles.axes3);
        set(handles.axes3,'visible','on');
        bar(0:255,his,'g');

       %% ================= Part 2: RETRIVE FEATURES ====================
        % ================= HISTOGRAM ================
        % Total number of histogram based features: 9
        HIS_FEAT = histogramfeature(imgData);

        % ================= GRADIENT ================
        % Total number of absolute gradient based features: 5
        GRA_FEAT = gradientfeature(imgGradientData);

        % ================= RUN LENGTH MATRIX ================
        % Features are computed for 4 (2D images) or 13 (3D images) various directions.
        % Total number of run length matrix based features: 44 (2D) or 143 (3D)
        [GLRLMS,SI] = grayrlmatrix(img,'NumLevels',64,'G',[]);
        RL_STATS = grayrlprops(GLRLMS,4);

        % ================= COOCURRENCE MATRIX ================
        % Features are computed for 5 between-pixels distances (1, 2, 3, 4, 5) and for 4 (2D images) or 13 (3D images) various directions.
        % Total number of co-occurrence matrix based features: 220 (2D) or 715 (3D)
        for D = 1:5
            GLCMS{4*D-3} = graycomatrix(img,'Offset',[0,D],'NumLevels',64,'G',[]);
            GLCMS{4*D-2} = graycomatrix(img,'Offset',[-D,0],'NumLevels',64,'G',[]);
            GLCMS{4*D-1} = graycomatrix(img,'Offset',[-D,D],'NumLevels',64,'G',[]);
            GLCMS{4*D} = graycomatrix(img,'Offset',[D,D],'NumLevels',64,'G',[]);
        end
        CM_STATS = graycomyprops(GLCMS,20);

        % ================= AUTOREGRESSIVE MODEL  ================
        % Total number of autoregressive model based features: 5 
        % AUTOREG_FEAT = regfeature(bmp);

        % ================= HAAR WAVELET  ================
        % Feature is computed at 5 scales within four frequency bands LL, LH, HL and HH. 
        % Total number of Haar wavelet based features: 20 
        HAAR_FEAT = zeros(5,4);
        [c s] = wavedec2(img,5,'haar');

        for i = 1:5
            LL=appcoef2(c,s,'haar',i);
            LH=detcoef2('h',c,s,i);
            HL=detcoef2('v',c,s,i);
            HH=detcoef2('d',c,s,i); 
            HAAR_FEAT(i,:) = haarfeature(LL,LH,HL,HH);
        end
        
        %% ================= Part 3: OUTPUT FEATURES ====================
        xls_path = [hdrPathName(hdrIndex).pathName,'\output_L',num2str(imageNumber),'.xls'];

        Table = {'Histogram based features'};
        [status message] = xlswrite(xls_path,Table,'Sheet1','A1');
        TableTitle = {'Mean' 'Variance' 'Skewness' 'Kurtosis' 'Perc.01%' 'Perc.10%' 'Perc.50%' 'Perc.90%' 'Perc.99%'};
        [status message] = xlswrite(xls_path,TableTitle,'Sheet1','A2');
        TableData = HIS_FEAT;
        xlswrite(xls_path,TableData,'Sheet1','A3');

        Table = {'Gradient based features'};
        xlswrite(xls_path,Table,'Sheet1','A5');
        TableTitle = {'GrMean' 'GrVariance' 'GrSkewness' 'GrKurtosis' 'GrNonZeros'};
        xlswrite(xls_path,TableTitle,'Sheet1','A6');
        TableData = GRA_FEAT;
        xlswrite(xls_path,TableData,'Sheet1','A7');

        Table = {'Run Length Matrix based features'};
        xlswrite(xls_path,Table,'Sheet1','A9');
        % degree = 0
        TableTitle = {'dgr0_ShrtREmp' 'dgr0_LngREmph' 'dgr0_GLevNonU' 'dgr0_RLNonUni' 'dgr0_RunPerc' 'dgr0_LowGREmph' 'dgr0_HighGREmph' 'dgr0_ShrtRLowGREmph' 'dgr0_ShrtRHighGREmph' 'dgr0_LngRLowGREmph' 'dgr0_LngRHighGREmph'};
        xlswrite(xls_path,TableTitle,'Sheet1','A10');
        TableData = RL_STATS(1,:);
        xlswrite(xls_path,TableData,'Sheet1','A11');
        % degree = 90
        TableTitle = {'dgr90_ShrtREmp' 'dgr90_LngREmph' 'dgr90_GLevNonU' 'dgr90_RLNonUni' 'dgr90_RunPerc' 'dgr90_LowGREmph' 'dgr90_HighGREmph' 'dgr90_ShrtRLowGREmph' 'dgr90_ShrtRHighGREmph' 'dgr90_LngRLowGREmph' 'dgr90_LngRHighGREmph'};
        xlswrite(xls_path,TableTitle,'Sheet1','A12');
        TableData = RL_STATS(2,:);
        xlswrite(xls_path,TableData,'Sheet1','A13');
        % degree = 45
        TableTitle = {'dgr45_ShrtREmp' 'dgr45_LngREmph' 'dgr45_GLevNonU' 'dgr45_RLNonUni' 'dgr45_RunPerc' 'dgr45_LowGREmph' 'dgr45_HighGREmph' 'dgr45_ShrtRLowGREmph' 'dgr45_ShrtRHighGREmph' 'dgr45_LngRLowGREmph' 'dgr45_LngRHighGREmph'};
        xlswrite(xls_path,TableTitle,'Sheet1','A14');
        TableData = RL_STATS(3,:);
        xlswrite(xls_path,TableData,'Sheet1','A15');
        % degree = 135
        TableTitle = {'dgr135_ShrtREmp' 'dgr135_LngREmph' 'dgr135_GLevNonU' 'dgr135_RLNonUni' 'dgr135_RunPerc' 'dgr135_LowGREmph' 'dgr135_HighGREmph' 'dgr135_ShrtRLowGREmph' 'dgr135_ShrtRHighGREmph' 'dgr135_LngRLowGREmph' 'dgr135_LngRHighGREmph'};
        xlswrite(xls_path,TableTitle,'Sheet1','A16');
        TableData = RL_STATS(4,:);
        xlswrite(xls_path,TableData,'Sheet1','A17');

        Table = {'Co-occurrence Matrix based features'};
        xlswrite(xls_path,Table,'Sheet1','A19');

        line = 19;
        for i = 1:5
            str = num2str(i);
            TableTitle = {['S(',str,',0)AngScMom'] ['S(',str,',0)Contrast'] ['S(',str,',0)Correlat'] ['S(',str,',0)SumOfSqs'] ['S(',str,',0)InvDfMom'] ['S(',str,',0)SumAverg'] ['S(',str,',0)SumVarnc'] ['S(',str,',0)SumEntrp'] ['S(',str,',0)Entropy'] ['S(',str,',0)DifVarnc'] ['S(',str,',0)DifEntrp'] ['S(',str,',0)Homogeneity']};
            line = line + 1;
            line_num = strcat('A',num2str(line));
            xlswrite(xls_path,TableTitle,'Sheet1',line_num);
            line = line + 1;
            line_num = strcat('A',num2str(line));
            TableData = CM_STATS(4*i-3,:);
            xlswrite(xls_path,TableData,'Sheet1',line_num);

            TableTitle = {['S(0,',str,')AngScMom'] ['S(0,',str,')Contrast'] ['S(0,',str,')Correlat'] ['S(0,',str,')SumOfSqs'] ['S(0,',str,')InvDfMom'] ['S(0,',str,')SumAverg'] ['S(0,',str,')SumVarnc'] ['S(0,',str,')SumEntrp'] ['S(0,',str,')Entropy'] ['S(0,',str,')DifVarnc'] ['S(0,',str,')DifEntrp'] ['S(0,',str,')Homogeneity']};
            line = line + 1;
            line_num = strcat('A',num2str(line));
            xlswrite(xls_path,TableTitle,'Sheet1',line_num);
            line = line + 1;
            line_num = strcat('A',num2str(line));
            TableData = CM_STATS(4*i-2,:);
            xlswrite(xls_path,TableData,'Sheet1',line_num);

            TableTitle = {['S(',str,',',str,')AngScMom'] ['S(',str,',',str,')Contrast'] ['S(',str,',',str,')Correlat'] ['S(',str,',',str,')SumOfSqs'] ['S(',str,',',str,')InvDfMom'] ['S(',str,',',str,')SumAverg'] ['S(',str,',',str,')SumVarnc'] ['S(',str,',',str,')SumEntrp'] ['S(',str,',',str,')Entropy'] ['S(',str,',',str,')DifVarnc'] ['S(',str,',',str,')DifEntrp'] ['S(',str,',',str,')Homogeneity']};
            line = line + 1;
            line_num = strcat('A',num2str(line));
            xlswrite(xls_path,TableTitle,'Sheet1',line_num);
            line = line + 1;
            line_num = strcat('A',num2str(line));
            TableData = CM_STATS(4*i-1,:);
            xlswrite(xls_path,TableData,'Sheet1',line_num);

            TableTitle ={['S(',str,',-',str,')AngScMom'] ['S(',str,',-',str,')Contrast'] ['S(',str,',-',str,')Correlat'] ['S(',str,',-',str,')SumOfSqs'] ['S(',str,',-',str,')InvDfMom'] ['S(',str,',-',str,')SumAverg'] ['S(',str,',-',str,')SumVarnc'] ['S(',str,',-',str,')SumEntrp'] ['S(',str,',-',str,')Entropy'] ['S(',str,',-',str,')DifVarnc'] ['S(',str,',-',str,')DifEntrp'] ['S(',str,',-',str,')Homogeneity']}; 
            line = line + 1;
            line_num = strcat('A',num2str(line));
            xlswrite(xls_path,TableTitle,'Sheet1',line_num);
            line = line + 1;
            line_num = strcat('A',num2str(line));
            TableData = CM_STATS(4*i,:);
            xlswrite(xls_path,TableData,'Sheet1',line_num);
        end

        Table = {'Haar Wavelet based features'};
        xlswrite(xls_path,Table,'Sheet1','A61');
        line = 61;
        for i = 1:5
            str = num2str(i);
            TableTitle = {['WavEnLL_',str] ['WavEnLH',str] ['WavEnHL',str] ['WavEnHH',str]};
            line = line + 1;
            line_num = strcat('A',num2str(line));
            xlswrite(xls_path,TableTitle,'Sheet1',line_num);

            line = line + 1;
            line_num = strcat('A',num2str(line));
            TableData = CM_STATS(i,:);
            xlswrite(xls_path,TableData,'Sheet1',line_num);
        end   
    end    
end

