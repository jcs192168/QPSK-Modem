%% ECE 161C Cascaded Midterm Project
%Jean-Claude Sleiman
%PID: A12934981
clear all;
close all;

%Intial Setup 8 Samples per Symbol
n_num = 10000;
sym = (floor(2*rand(1,n_num))-0.5)/0.5+j*(floor(2*rand(1,n_num))-0.5)/0.5;
Delay = 6;
sps = 8;
input_signal = upsample(sym, sps);
h0 = rcosine(1, sps, 'sqrt', 0.5, Delay); %pulse shaping
h0 = h0/max(h0); % normalize h0
h1=conv(h0,h0)/(h0*h0');  % convolution of h0 with itself --> h1
pulse_shaped_signal = filter(h1, 1, input_signal);
pulse_shaped_signal = downsample(pulse_shaped_signal,2); %8-->4 sps
pulse_shaped_signal = pulse_shaped_signal';

%Noise Control
noise_var = 0;
noisy_signal = 0.2*noise_var*(randn(1,length(pulse_shaped_signal))+1j*randn(1,length(pulse_shaped_signal)));
received_signal = pulse_shaped_signal + noisy_signal';

%Channel & Frequency Offset
ch = [1 0 0 0 0 0.2 0 0 0 0 0 0 1j*.05];
received_signal = filter(ch,1,received_signal);
received_signal = downsample(received_signal(2:end),2); %4-->2 sps
received_signal = received_signal .* exp(1j*2*pi*(1:length(received_signal))*0.1)';
%% FREQUENCY LOCK LOOP (FLL) (Howard)
h = rcosine(1,4,'sqrt',0.5,6);
h = h/max(h);
hm= h(1:2:end);

h2 = hm/(hm*hm');
hilb=j*[-1/5 0 -1/13 0 -1/11 0 -1/9 0 -1/7 0 -1/5 0 -1/3 0 -1 0 1 0 1/3 0 1/5 0 1/7 0 1/9 0 1/11 0 1/13 0 1/15]*1/pi;
hilb(16)=1/2;
hilb=hilb.*kaiser(31,5.5)';
gg1=j*h2.*(-3:1/4:3).*kaiser(25,4)';
gg1a=gg1;
gg1b=conv(gg1a,hilb);
g1=gg1b(16:55-15);
% link to Cid
x9 = received_signal;

fn=6;
fs=300;
eta=sqrt(2)/2;
theta=(pi*fn/fs);
k_i_f=4*theta*theta/(1+2*eta*theta+theta*theta);
k_p_f=4*eta*theta/(1+2*eta*theta+theta*theta);
reg=zeros(1,25);
accum=0;
int=0;
x7 = zeros(1,length(x9));
for nn=1:length(x7)
    x7(nn)=x9(nn)*exp(-j*2*pi*accum);
    reg=[x7(nn) reg(1:24)];
    be_1(nn)=reg*g1';
    be_2(nn)=reg*conj(g1)';
    f_err=be_2(nn)*conj(be_2(nn))-be_1(nn)*conj(be_1(nn));
    %f_err=abs(be_1(nn))-abs(be_2(nn));
    int=int+k_i_f*f_err;
    lpfltr(nn)=int+k_p_f*f_err;
    accum_sv(nn)=accum;
    accum=accum+lpfltr(nn);
end
%% TIMING LOOP (TL) (JC)
hh=rcosine(1,4,'sqrt',0.5,6);  %49 taps
hh=hh/max(hh);
n_dat=n_num;
% Howard Link
x3 = x7;
% 32 Polyphase matched filters operating at 2-smpls/sym
reg=zeros(1,21);
g=rcosine(1,64,'sqrt',0.5,6);  %769 taps
g=g/(hh(1:2:49)*g(1:32:769)');
dgx=conv(g,[1 0 -1]*64/2);
dg=dgx(2:770);
h2=reshape([0 0 0 hh(1:49)],4,13); 
g2=reshape([zeros(1,15) g(1:769) zeros(1,16)],32,25); % [hh_t]Polyphase interpolator filter, 25 taps per filter
dg2=reshape([zeros(1,15) dg(1:769) zeros(1,16)],32,25); %polyphase derivative filter

theta_0= 2*pi/200; %change BW
eta=sqrt(2)/2;
eta=4*eta;

k_i_t= (4*theta_0*theta_0)/(1+2*eta*theta_0+theta_0*theta_0);
k_p_t= (4*eta*theta_0)/(1+2*eta*theta_0+theta_0*theta_0);

reg_t=zeros(1,25); %interoplator register
int_t=0.0;

ndx_strt=22;   % Can Change starting index

accum_t=ndx_strt; %initial index stored in accumulator
accum_t_sv=zeros(1,n_dat); %save accumulator

pp=+1;                          % pp = +1 or -1 switches starting intervals
mm=1;                           % output clock at 1-sample per symbol
for nn=20:2:length(x3)-2
    reg_t=[x3(nn) reg_t(1:24)];     % new sample in matched filter register
    pntr=floor(accum_t);            % point to a coefficient set
    y_t1 =reg_t*g2(pntr,:)';        % MF output time sample
    dy_t1=reg_t*dg2(pntr,:)';       % dMF output time sample
    x6(nn)=y_t1;                    % save MF output sample
    
    reg_t=[x3(nn+1) reg_t(1:24)];   % new sample in matched filter register
    y_t2=reg_t*g2(pntr,:)';      % [y2] point to a coefficient set
    y2_sv(nn) = y_t2;
    x6(nn+1)=y_t2;                   % MF output time sample
    dy_t2=reg_t*dg2(pntr,:)';    % dMF output time sample
    
    det_t=pp*real(y_t2)*real(dy_t2);  % y*y_dot product (timing error)
    det_t_sv(mm)=det_t;             % save timing error
    int_t=int_t+k_i_t*det_t;        % Loop filter integrator
    sm_t =int_t+k_p_t*det_t;        % loop filter output (smoothing to estimate DC level)
    
    accum_t_sv(mm)=accum_t;         % save timing accumulator content
    mm=mm+1;                        % increment symbol clock
    accum_t=accum_t+sm_t;           % update accumulator (save DC, accum starts rising)
    if accum_t>33                   % test for accumulator overflow
        accum_t=accum_t-32;
        pp=-pp;
    end
    
    if accum_t<1                    % test for accumulator underflow
        accum_t=accum_t+32;
        pp=-pp;
    end
end
%% PHASE LOCK LOOP (PLL) /EQUALIZER (VIDYA/CID)
% JC Link
received_signal = x6(1:end-1);
n_num = n_dat;
% LMS Equalizer
% PLL section
% parameters for phase PLL
theta_0= 2*pi/200;
eta = sqrt(2)/2;
denom = (1+2*eta*theta_0+theta_0*theta_0);
k_i = (4*theta_0*theta_0)/denom;
k_p = (4*eta*theta_0)/denom;

% Initialize Internal registers
accum=0;
int=0;
hld=1;
reg=zeros(1,20);
wts=zeros(1,20);
wts(1)=1;

% Arrays to store results for figures
accum_sv=zeros(1,2*n_num);
err_sv=zeros(1,n_num);
lp_sv=zeros(1,n_num);

% output equalizer and Phase locked
y=zeros(1,2*n_num);

m=1;
mu=0.01;
for n=1:2:length(received_signal)
    prd=received_signal(n)*exp(-j*2*pi*accum);    % input down convert
    reg=[prd reg(1:19)];                      % insert in equalizer filter
    yy=reg*wts';                                  % compute equalizer filter output
    y(n)=yy;                                        % save output
    if abs(yy)>0.1;                            % test for signal strength
        y_det=sign(real(yy))+j*sign(imag(yy));  % slice (detect) yy
    else
        y_det=0+j*0;
    end
    yy_sv(m)=y_det;
    det_prd=conj(yy)*y_det;            % form conjugate product
    det=-angle(det_prd)/(2*pi);        % measure angle error
    err_sv(m)=det;                           % save angle error
    int=int+k_i*det;                           % loop filter integrator
    lp=int+k_p*det;                          % loop filter output
    lp_sv(m)=lp;                              % save loop filter output
    eq_err=yy-y_det;                       % equalizer error
    eq_err_sv(m)=eq_err;               % save equalizer error
    wts=wts-mu*conj(eq_err)*reg;  % update weights
    m=m+1;                                    % increment symbol clock
    hld=1*lp/2;                                % scale lp and store in hld
    accum_sv(n)=accum;               % save accumulator
    accum=accum+hld;                  % Increment accumulator
    
    prd=received_signal(n+1)*exp(-j*2*pi*accum);  % second input sample, not symbol sample
    reg=[prd reg(1:19)];                         % loop filter (sample rate)
    yy=reg*wts';
    y(n+1)=yy;
    accum_sv(n+1)=accum;
    accum=accum+hld;
end

%Final Plots
scatterplot(received_signal(1:2:end)) %input signal
scatterplot(y(1:2:end)) %output signal

