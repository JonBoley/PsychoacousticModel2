function varargout = MODEL2GUI(varargin)
% MODEL2GUI Application M-file for MODEL2GUI.fig
% Simulates MPEG Psychoacoustic Model 2
% (Note: Some small liberties were taken, and this code does not
% adhere strictly to the standard.  However, it should still work
% well as an instructional tool.)

% Original Code by Jonathan Boley and Vishweshwara Rao
% Music Engineering Department, University of Miami
% GUI adaptation by Jonathan Boley (10 Dec 2003)
% To report bugs or fixes, please send email to jdb@jboley.com

global Y
if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);

	if nargout > 0
		varargout{1} = fig;
	end
    t=1/44100:1/44100:1;
    Y = sin(2*pi*t*5000);
    Initialize(1);

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch
		disp(lasterr);
	end

end

% --------------------------------------------------------------------
function Initialize(l)

WriteTables(1);

% Initialization
global iblen Fs earlyblock prevblock r f iblen_index
iblen = 512;
Fs = 44100;
earlyblock = [];
prevblock = [];
r = zeros(3,1024);
f = zeros(3,1024);

for iblen_index = 0:1
    CalculateAll(1)
end % End Main Loop

% --------------------------------------------------------------------
function varargout = InputSelect_Callback(h, eventdata, handles, varargin)
% New Input Function
Input = get(h,'Value');

clear global;   % Clear all global variables
WriteTables(1);

% Initialization
global iblen Fs earlyblock prevblock r f Y
iblen = 512;
Fs = 44100;
earlyblock = [];
prevblock = [];
r = zeros(3,1024);
f = zeros(3,1024);

switch Input
case 1
    t=1/44100:1/44100:1;
    Y = sin(2*pi*t*5000);
    Initialize(1);
case 2
    t=1/Fs:1/Fs:1;                 % 1 sec @ 44.1kHz sample rate
    freqs = [50 300 725 1200 1800 ...
        4000 7000 12500];
    Y=zeros(1,Fs);
    for i=1:8, Y=Y+(sin(2*pi*t*freqs(i)))/8; end;
    Initialize(1);
case 3
    t=1/44100:1/44100:1;
    Y = rand(1,44100);
    [B,A]=ellip(2,5, 40,55/512);
    Y = filter(B,A,Y);
    Y = Y + sin(2*pi*t*6000);
    Initialize(1);
case 4
    name = uigetfile('*.wav', 'Please choose a WAV file');
    if name ~= 0
        X = wavread(name,[1 100000]);
        Y=X';
        Initialize(1);
    end
end

% Re-Plot
value = get(handles.PlotSelect,'Value');
Replot(value);

% --------------------------------------------------------------------
function varargout = PlotSelect_Callback(h, eventdata, handles, varargin)
% New plot
value = get(h,'Value');
Replot(value);

% --------------------------------------------------------------------
function varargout = NextFrame_Callback(h, eventdata, handles, varargin)
% Next Frame

global iblen iblen_index Y

if iblen_index <= floor(length(Y)/iblen)
    iblen_index = iblen_index + 1;
    CalculateAll(1);
    
    % Re-Plot
    value = get(handles.PlotSelect,'Value');
    Replot(value);
else
    error('There are no more samples left in this signal');
end % End Main If Statement

% --------------------------------------------------------------------
function CalculateAll(l)
% Calculate All Variables

global iblen Fs earlyblock prevblock newblock r f
global Y s r_hat f_hat cw e cb en cbb spreadplot epart npart
global tbb SNRb bcb nbb nbw thrw SMR THRw SMR2 Fs iblen_index

newblock = Y((iblen_index*iblen)+1:(iblen_index*iblen)+iblen);

% Step 1 - Reconstruct 1024 Samples of the Input Signal
[earlyblock, prevblock, newblock, s] = ...
    Reconstruct(earlyblock, prevblock, newblock, iblen);

if length(s) == 1024
    % Step 2 - Calculate the Complex Spectrum of the Input Signal
    [r,f] = Spectrum(s, r, f);

    % Step 3 - Calculate a Predicted r and f
    [r_hat, f_hat] = Predict(r, f);

    % Step 4 - Calculate the Unpredictability Measure cw
    cw = CalculateCw(r, f, r_hat, f_hat);

    % Step 5 - Calculate the Energy and Unpredictability in the Threshold Calculation Partitions
    [e, cb] = Energy_Unpredictability(r, cw);
        
    % Step 6 - Convolve the Partitioned Energy and Unpredictability with the Spreading Function
    load ('tables.mat'); % Tables and Spreading Functions
    [en, cbb] = Spread(e, cb);

    % Step 7 - Convert cbb to tbb, the Tonality Index
    tbb = CalcTonality(cbb);

    % Step 8 - Calculate the Required SNR in Each Partition
    SNRb = CalcSNR(tbb);

    % Step 9 - Calculate the Power Ratio
    bcb = CalcPwrRatio(SNRb);

    % Step 10 - Calculation of Actual Energy Threshold, nbb
    nbb = CalcNbb(en, bcb);

    % Step 11 - Spread the Threshold Energy over FFT Lines, Yielding nbw
    nbw = CalcNb(nbb);

    % Step 12 - Include Absolute Thresholds, Yielding the Final Energy Threshold of Audibility, thrw
    thrw = CalcThresh(nbw);

    % Step 13 - Pre-Echo Control
    % NOT USED FOR LAYERS 1,2

    % Step 14 - Calculate the Signal-to-Mask Ratios, SMRn
    [SMR(iblen_index+1,:),epart,npart] = CalcSMR(r, thrw);
end

% --------------------------------------------------------------------
function Replot(value)
% Re-plot
global Y s r f r_hat f_hat cw e cb en cbb spreadplot tbb SNRb bcb
global nbb nbw thrw SMR THRw SMR2 Fs epart npart iblen_index
load ('tables.mat'); % Tables and Spreading Functions
figure(2); clf;
switch value
case 1
    plot(s); title('Reconstructed 1024 Samples');
    xlabel('Sample Number'); ylabel('Amplitude');
case 2
    subplot(2,1,1), semilogx([0:Fs/1024:(Fs/2)],r(3,1:513)); title('Magnitude of FFT'); 
    subplot(2,1,2), semilogx([0:Fs/1024:(Fs/2)],f(3,1:513)); title('Phase of FFT');
    xlabel('Frequency (Hz)');
case 3
    subplot(2,1,1), semilogx([0:Fs/1024:(Fs/2)],r_hat(1:513)); title('Predicted Magnitude');
    subplot(2,1,2), semilogx([0:Fs/1024:(Fs/2)],f_hat(1:513)); title('Predicted Phase');
    xlabel('Frequency (Hz)');
case 4
    semilogx([0:Fs/1024:(Fs/2)],cw(1:513)); title('Unpredictability Measure');
    xlabel('Frequency (Hz)');
case 5
    for i=1:57, spreadplot(i,:) = sprdngf(i,:)*e(i); end
    subplot(2,1,1), plot(spreadplot', 'g:'); hold on; plot(e, 'b'); hold off;
    title('Energy (blue) and Spreading Functions (green)');
    subplot(2,1,2), plot(en); title('Energy w/ Spreading');
    xlabel('Sub-Band Number');
case 6
    for i=1:57, spreadplot(i,:) = sprdngf(i,:)*cb(i); end
    subplot(2,1,1), plot(spreadplot', 'g:'); hold on; plot(cb, 'b'); hold off;
    title('Unpredictability (blue) and Spreading Functions (green)');
    subplot(2,1,2), plot(cbb);
    title('Unpredictability w/ Spreading');
    xlabel('Sub-Band Number');
case 7
    plot(tbb); title('Tonality Index'); xlabel('Sub-Band Number');
case 8
    plot(SNRb); title('SNR'); xlabel('Sub-Band Number');
case 9
    plot(bcb); title('Power Ratio'); xlabel('Sub-Band Number');
case 10
    plot(nbb); title('Actual Energy Threshold'); xlabel('Sub-Band Number');
case 11
    semilogx([0:Fs/1024:(Fs/2)],abstable, 'b'); hold on;
    semilogx([0:Fs/1024:(Fs/2)],nbw, 'r');
    title('Threshold Energy (red) and Threshold in Quiet (blue)');
    xlabel('Frequency (Hz)');
    hold off;
case 12
    semilogx([0:Fs/1024:(Fs/2)],thrw); title('Absolute Threshold');
    xlabel('Frequency (Hz)');
case 13
    % Cut off everything < 0 (just for display)
    epart = 10*log10(epart);
    npart = 10*log10(npart);
    SMR2 = SMR(iblen_index+1,:);
    epart(find(epart<0))=0;
    npart(find(npart<0))=0;
    SMR2(find(SMR2<0))=0;

    subplot(2,1,1), semilogx(epart, 'r'); hold on; plot(real(npart), 'g'); hold off;
    title('Energy (Red) and Masking Threshold (Green)');
    subplot(2,1,2), semilogx(SMR2); title('SMR');
    xlabel('Sub-Band Number');
end

%-------------------------

function [earlyblock, prevblock, newblock, s] = ...
    Reconstruct(earlyblock, prevblock, newblock, iblen);
if iblen >= 512, block = [prevblock newblock];
else block = [earlyblock prevblock newblock]; end
earlyblock=prevblock;
prevblock=newblock;
if length(block) >= 1024
    s = block(end-1023:end);  % Newest 1024 samples
else
    s = [];
end

%-------------------------

function [r,f] = Spectrum(s, r, f)
sw = s .* (0.5 - 0.5*cos((2*pi*([1:1024]-0.5))/1024));  % Hann Window
r(1:2,:) = r(2:3,:); % Shift previous two magnitude values
f(1:2,:) = f(2:3,:); % Shift previous two phase values
r(3,:) = abs(fft(sw));
f(3,:) = angle(fft(sw));

%-------------------------

function [r_hat, f_hat] = Predict(r, f)
r_hat = 2.0 * r(2,:) - r(1,:);
f_hat = 2.0 * f(2,:) - f(1,:);

%-------------------------

function cw = CalculateCw(r, f, r_hat, f_hat)
for i=1:1024
    cw(i) = sqrt((r(3,i)*cos(f(3,i)) - r_hat(i)*cos(f_hat(i)))^2 + (r(3,i)*sin(f(3,i)) ...
        - r_hat(i)*sin(f_hat(i)))^2)/(r(3,i)+abs(r_hat(i)));
end

%-------------------------

function [e, cb] = Energy_Unpredictability(r, cw)
load ('tables.mat'); % Tables and Spreading Functions
w_lo(1:57) = Table3D3b([1:57],1);
w_hi(1:57) = Table3D3b([1:57],2);
for w=1:57
    e(w) = sum((r(3,w_lo(w):w_hi(w))).^2);
    cb(w)= sum(cw(w_lo(w):w_hi(w)).*(r(3,w_lo(w):w_hi(w))).^2);
end

%-------------------------

function [en, cbb] = Spread(e, cb)
load ('tables.mat'); % Tables and Spreading Functions
ecb=zeros(1,57); ct=zeros(1,57);
for i=1:57
    for j=1:57
        ecb(i) = ecb(i) + (e(j) * sprdngf(j,i));
        ct(i) = ct(i) + (cb(j) * sprdngf(j,i));
    end
end
cbb = ct;%./ecb;

rnormb = 1 ./ (sum(sprdngf,1));
en = ecb.*rnormb;

%-------------------------

function tbb = CalcTonality(cbb)
tbb = -0.299 - 0.43*log(cbb);
tbb=tbb/(max(tbb)-min(tbb));
tbb=1-(tbb-min(tbb));

%-------------------------

function SNRb = CalcSNR(tbb)
load ('tables.mat'); % Tables and Spreading Functions
NMTb = 5.5;  % Downshift for Noise Masking Tone (in dB)
SNRb = max(Table3D3b(:,4)', tbb .* Table3D3b(:,5)'+(1-tbb)*NMTb);

%-------------------------

function bcb = CalcPwrRatio(SNRb)
bcb = 10.^(-SNRb/10);

%-------------------------

function nbb = CalcNbb(en, bcb)
nbb = en .* bcb;

%-------------------------

function nb = CalcNb(nbb)
load ('tables.mat'); % Tables and Spreading Functions
w_lo(1:57) = Table3D3b([1:57],1);
w_hi(1:57) = Table3D3b([1:57],2);
for b=1:57
    for w=1:513
        if ((w>=w_lo(b))&(w<=w_hi(b)))
            nb(w) = nbb(b)/(w_hi(b)-w_lo(b)+1);
        end;
    end
end
% TRANSFORM FFT'S TO SPL VALUES
fftmax = 61559;  % max(abs(fft(1kHz tone)))... defined as 96dB
% sig=sin(2*pi*1000*[1/44100:1/44100:1024/44100]);
% win=(0.5 - 0.5*cos((2*pi*([1:1024]-0.5))/1024));
% fftmax = max(fft(sig.*win).*conj(fft(sig.*win))); % 61559
nb = 96 + 10*log10(abs(nb)/fftmax);

%-------------------------

function thrw = CalcThresh(nb)
load ('tables.mat'); % Tables and Spreading Functions
thrw = max(nb, abstable); % both nb and abstable are in dBSPL

%-------------------------

function [SMR,epart,npart] = CalcSMR(r, thrw)
w_low=1;
for i=1:31
    epart(i)=sum((r(3,w_low:w_low+16)).^2);
    if i < 13
        npart(i) = sum(thrw(w_low:w_low+16));
    else
        npart(i) = min(thrw(w_low:w_low+16)) * 17;
    end
    w_low=w_low+16;
end
SMR = 10*log10(epart)-npart;

%-------------------------

function WriteTables(a)
% Table 3-D.3b "Calculation Partition Table"
%             w(low) | w(high) | bval  |  minval  |  TMN
global Table3D3b abstable sprdngf
Table3D3b = [   1       1       0.00        0.0     24.5;...
                2       2       0.43        0.0     24.5;...
                3       3       0.86        0.0     24.5;...
                4       4       1.29        20.0    24.5;...
                5       5       1.72        20.0    24.5;...
                6       6       2.15        20.0    24.5;...
                7       7       2.58        20.0    24.5;...
                8       8       3.01        20.0    24.5;...
                9       9       3.45        20.0    24.5;...
                10      10      3.88        20.0    24.5;...
                11      11      4.28        20.0    24.5;...
                12      12      4.67        20.0    24.5;...
                13      13      5.06        20.0    24.5;...
                14      14      5.42        20.0    24.5;...
                15      15      5.77        20.0    24.5;...
                16      16      6.11        17.0    24.5;...
                17      19      6.73        17.0    24.5;...
                20      22      7.61        15.0    24.5;...
                23      25      8.44        10.0    24.5;...
                26      28      9.21        7.0     24.5;...
                29      31      9.88        7.0     24.5;...
                32      34      10.51       4.4     25.0;...
                35      37      11.11       4.5     25.6;...
                38      40      11.65       4.5     26.2;...
                41      44      12.24       4.5     26.7;...
                45      48      12.85       4.5     27.4;...
                49      52      13.41       4.5     27.9;...
                53      56      13.94       4.5     28.4;...
                57      60      14.42       4.5     28.9;...
                61      64      14.86       4.5     29.4;...
                65      69      15.32       4.5     29.8;...
                70      74      15.79       4.5     30.3;...
                75      80      16.26       4.5     30.8;...
                81      86      16.73       4.5     31.2;...
                87      93      17.19       4.5     31.7;...
                94      100     17.62       4.5     32.1;...
                101     108     18.05       4.5     32.5;...
                109     116     18.45       4.5     32.9;...
                117     124     18.83       4.5     33.3;...
                125     134     19.21       4.5     33.7;...
                135     144     19.60       4.5     34.1;...
                145     155     20.00       4.5     34.5;...
                156     166     20.38       4.5     34.9;...
                167     177     20.74       4.5     35.2;...
                178     192     21.12       4.5     35.6;...
                193     207     21.48       4.5     36.0;...
                208     222     21.84       4.5     36.3;...
                223     243     22.20       4.5     36.7;...
                244     264     22.56       4.5     37.1;...
                265     286     22.91       4.5     37.4;...
                287     314     23.26       4.5     37.8;...
                315     342     23.60       4.5     38.1;...
                343     371     23.95       4.5     38.4;...
                372     401     24.30       4.5     38.8;...
                402     431     24.65       4.5     39.1;...
                432     469     25.00       4.5     39.5;...
                470     513     25.33       3.5     39.8 ];
        
abstable=[  45.05           25.87           18.70           14.85           12.41           10.72 ...
            9.47            8.50            7.73            7.1             6.56            6.11 ...
            5.72            5.37            5.07            4.79            4.55            4.32 ...
            4.11            3.92            3.74            3.57            3.4             3.25 ...
            3.1             2.95            2.81            2.67            2.53            2.39 ...
            2.25            2.11            1.97            1.83            1.68            1.53 ...
            1.38            1.23            1.07            .9              .74             .56  ...
            .39             .21             .02             -.17            -.36            -.56 ...
            -.96            -.96            -1.37           -1.37           -1.79           -1.79 ...
            -2.21           -2.21           -2.63           -2.63           -3.03           -3.03 ...
            -3.41           -3.41           -3.77           -3.77           -4.09           -4.09 ...
            -4.37           -4.37           -4.6            -4.6            -4.78           -4.78 ...
            -4.91           -4.91           -4.97           -4.97           -4.98           -4.98 ...
            -4.92           -4.92           -4.81           -4.81           -4.65           -4.65 ...
            -4.43           -4.43           -4.17           -4.17           -3.87           -3.87 ...
            -3.54           -3.54           -3.19           -3.19           -2.82           -2.82 ...
            -2.06*ones(1,4) -1.33*ones(1,4) -.64*ones(1,4)  -.04*ones(1,4)  0.47*ones(1,4)  .89*ones(1,4) ...
            1.23*ones(1,4)  1.51*ones(1,4)  1.74*ones(1,4)  1.93*ones(1,4)  2.11*ones(1,4)  2.28*ones(1,4) ...
            2.45*ones(1,4)  2.63*ones(1,4)  2.82*ones(1,4)  3.03*ones(1,4)  3.25*ones(1,4)  3.49*ones(1,4) ...
            3.74*ones(1,4)  4.02*ones(1,4)  4.32*ones(1,4)  4.64*ones(1,4)  4.98*ones(1,4)  5.35*ones(1,4) ...
            6.15*ones(1,8)  7.07*ones(1,8)  8.10*ones(1,8)  9.25*ones(1,8)  10.54*ones(1,8) 11.97*ones(1,8) ...
            13.56*ones(1,8) 15.3*ones(1,8)  17.23*ones(1,8) 19.33*ones(1,8) 21.64*ones(1,8) 24.15*ones(1,8) ...
            26.88*ones(1,8) 29.84*ones(1,8) 33.04*ones(1,8) 36.51*ones(1,8) 40.24*ones(1,8) 44.26*ones(1,8) ...
            48.58*ones(1,8) 53.21*ones(1,8) 58.17*ones(1,8) 63.48*ones(1,8) 69.13*ones(1,96)];
 abstable(465:513)=69.13;
 
 % The Spreading Function
% i = bark value of the signal being spread
% j = bark value of the band being spread into
bval = Table3D3b(:,3);
for i=1:57,
    for j=1:57,
        tmpx(i,j) = 1.05*(bval(j)-bval(i));
        x(i,j) = 8 * min((tmpx(i,j)-0.5)^2 - 2*(tmpx(i,j)-0.5),0);
        tmpy(i,j) = 15.811389 + 7.5*(tmpx(i,j)+0.474) - 17.5*sqrt(1.0+(tmpx(i,j)+0.474)^2);
        if tmpy(i,j) < -100
            sprdngf(i,j)=0;
        else
            sprdngf(i,j)=10^((x(i,j)+tmpy(i,j))/10);
        end
    end
end
 
 save('tables.mat', 'Table3D3b', 'abstable', 'sprdngf');
 