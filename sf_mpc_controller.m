function sf_mpc_controller(s)
setup(s);

function setup(s)
s.NumDialogPrms = 4;

s.NumInputPorts  = 5;
s.NumOutputPorts = 4;

s.SetPreCompInpPortInfoToDynamic;
s.SetPreCompOutPortInfoToDynamic;

s.InputPort(1).Dimensions = 1;
s.InputPort(1).DirectFeedthrough = true;
s.InputPort(2).Dimensions = 1;
s.InputPort(2).DirectFeedthrough = true;
s.InputPort(3).Dimensions = 1;
s.InputPort(3).DirectFeedthrough = true;
s.InputPort(4).Dimensions = 1;
s.InputPort(4).DirectFeedthrough = true;
s.InputPort(5).Dimensions = 1;
s.InputPort(5).DirectFeedthrough = true;

s.OutputPort(1).Dimensions = 1;
s.OutputPort(2).Dimensions = 1;
s.OutputPort(3).Dimensions = 1;
s.OutputPort(4).Dimensions = 1;

s.SampleTimes = [s.DialogPrm(1).Data 0];
s.SimStateCompliance = 'DefaultSimState';
s.RegBlockMethod('SetInputPortSamplingMode', @SetInpPortFrameData);
s.RegBlockMethod('Outputs', @Outputs);

function SetInpPortFrameData(block, idx, fd)
  block.InputPort(idx).SamplingMode = fd;
  block.OutputPort(1).SamplingMode  = fd;
  block.OutputPort(2).SamplingMode  = fd;
  block.OutputPort(3).SamplingMode  = fd;
  block.OutputPort(4).SamplingMode  = fd;

function Outputs(s)
Te = s.DialogPrm(1).Data;
Np = s.DialogPrm(2).Data; Nu = s.DialogPrm(3).Data;
L = s.DialogPrm(4).Data;

vR = s.InputPort(1).Data; omegaR = s.InputPort(2).Data;
e0 = [s.InputPort(3).Data;s.InputPort(4).Data;s.InputPort(5).Data];
Ad = ...
[  cos(Te*omegaR), sin(Te*omegaR), (vR - vR*cos(Te*omegaR))/omegaR;
  -sin(Te*omegaR), cos(Te*omegaR),     (vR*sin(Te*omegaR))/omegaR;
              0,            0,                          1];
Bd = ...
[          -sin(Te*omegaR)/omegaR, (vR*(sin(Te*omegaR) - Te*omegaR))/omegaR^2;
  (2*sin((Te*omegaR)/2)^2)/omegaR,    -(2*vR*sin((Te*omegaR)/2)^2)/omegaR^2;
                             0,                                   -Te];

Q = {100*eye(3),1*eye(3),diag([100 100 1]),eye(3)};
R = {1*eye(2),100*eye(2),diag([1 10]),eye(2)};

i_cas = 1;
% Utiliser Q{i_cas}, R{i_cas}

A_COMPLETER = nan;

S{i_cas} = A_COMPLETER;

% ...
 
% Critère 1/2 xi^T H xi où xi = [E^T U^T]^T
% avec E = [e[1];...;e[N]] et U = [uB[1];...;uB[N-1]]
H = A_COMPLETER;
f = A_COMPLETER;

%
% Contraintes égalités
Ae = A_COMPLETER;
Be = A_COMPLETER;

%
% Contraintes inégalités
uMIN = [-10;-10]; uMAX = [10;10];
zMIN = -pi; zMAX = pi;

Ai = A_COMPLETER;
Bi = A_COMPLETER;

%
% ENSUITE DECOMMENTER : [xiEU,fval,exitflag] = quadprog(H,f,Ai,Bi,Ae,Be);
% if (exitflag ~= 1), exitflag, error('Revoir programme optim'); end

uB = [vR;omegaR]; %MODIFIER BIEN SUR...
uF = [0;0]; %A_COMPLETER;

s.OutputPort(1).Data = uB(1)+uF(1);
s.OutputPort(2).Data = uB(2)+uF(2);
s.OutputPort(3).Data = uB(1);
s.OutputPort(4).Data = uB(2);
