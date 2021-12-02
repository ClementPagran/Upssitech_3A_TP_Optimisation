function sf_nonholonomic_robot(s)
setup(s);

function setup(s)
s.NumDialogPrms = 1;

s.NumInputPorts  = 2;
s.NumOutputPorts = 3;

s.SetPreCompInpPortInfoToDynamic;
s.SetPreCompOutPortInfoToDynamic;

s.InputPort(1).Dimensions = 1;
s.InputPort(1).DirectFeedthrough = false;
s.InputPort(2).Dimensions = 1;
s.InputPort(2).DirectFeedthrough = false;

s.OutputPort(1).Dimensions = 1;
s.OutputPort(2).Dimensions = 1;
s.OutputPort(3).Dimensions = 1;

s.SampleTimes = [0 0];

s.NumContStates = 3;

s.SimStateCompliance = 'DefaultSimState';

s.RegBlockMethod('SetInputPortSamplingMode', @SetInpPortFrameData);
  
s.RegBlockMethod('InitializeConditions', @InitializeConditions);
s.RegBlockMethod('Outputs', @Outputs);     % Required
s.RegBlockMethod('Derivatives', @Derivatives);

function SetInpPortFrameData(block, idx, fd)
  block.InputPort(idx).SamplingMode = fd;
  block.OutputPort(1).SamplingMode  = fd;
  block.OutputPort(2).SamplingMode  = fd;
  block.OutputPort(3).SamplingMode  = fd;

function InitializeConditions(s)
s.ContStates.Data = s.DialogPrm(1).Data;

function Outputs(s)
s.OutputPort(1).Data = s.ContStates.Data(1);
s.OutputPort(2).Data = s.ContStates.Data(2);
s.OutputPort(3).Data = s.ContStates.Data(3);

function Derivatives(s)
v = s.InputPort(1).Data;
omega = s.InputPort(2).Data;
x = s.ContStates.Data(1);
y = s.ContStates.Data(2);
theta = s.ContStates.Data(3);
s.Derivatives.Data = [v*cos(theta);v*sin(theta);omega];
