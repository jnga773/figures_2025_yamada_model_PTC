function varargout=F_bcs_PO(action,varargin)
%% Automatically generated with matlabFunction
%#ok<*DEFNU,*INUSD,*INUSL>

switch action
  case 'nargs'
   varargout{1}=1;
   return
  case 'nout'
   varargout{1}=4;
   return
  case 'argrange'
   varargout{1}=struct('u',1:10);
   return
  case 'argsize'
   varargout{1}=struct('u',10);
   return
  case 'vector'
   varargout{1}=struct('u',1);
   return
  case 'extension'
   varargout{1}='rhs';
   return
  case 'maxorder'
   varargout{1}=2;
   return
end
nout=4;
order=varargin{1};
f=str2func(sprintf('F_bcs_PO_%s_%d',action,order));
varargout=cell(nout,1);
[varargout{:}]=f(varargin{2:end});
end



function [out1,out2,out3,out4] = F_bcs_PO_rhs_0(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12,in13,in14,in15,in16,in17,in18,in19,in20)
%F_bcs_PO_rhs_0
%    [OUT1,OUT2,OUT3,OUT4] = F_bcs_PO_rhs_0(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12,IN13,IN14,IN15,IN16,IN17,IN18,IN19,IN20)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    16-Nov-2024 18:10:21

out1 = in1-in4;
if nargout > 1
    out2 = in2-in5;
end
if nargout > 2
    out3 = in3-in6;
end
if nargout > 3
    out4 = -in7.*(in1-in8+in1.*in3);
end
end


function [out1,out2,out3,out4] = F_bcs_PO_rhs_1(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12,in13,in14,in15,in16,in17,in18,in19,in20)
%F_bcs_PO_rhs_1
%    [OUT1,OUT2,OUT3,OUT4] = F_bcs_PO_rhs_1(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12,IN13,IN14,IN15,IN16,IN17,IN18,IN19,IN20)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    16-Nov-2024 18:10:22

out1 = in11-in14;
if nargout > 1
    out2 = in12-in15;
end
if nargout > 2
    out3 = in13-in16;
end
if nargout > 3
    out4 = -in17.*(in1-in8+in1.*in3)-in7.*(in11-in18+in1.*in13+in3.*in11);
end
end


function [out1,out2,out3,out4] = F_bcs_PO_rhs_2(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12,in13,in14,in15,in16,in17,in18,in19,in20)
%F_bcs_PO_rhs_2
%    [OUT1,OUT2,OUT3,OUT4] = F_bcs_PO_rhs_2(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12,IN13,IN14,IN15,IN16,IN17,IN18,IN19,IN20)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    16-Nov-2024 18:10:22

out1 = 0.0;
if nargout > 1
    out2 = 0.0;
end
if nargout > 2
    out3 = 0.0;
end
if nargout > 3
    out4 = in17.*(in11-in18+in1.*in13+in3.*in11).*-2.0-in7.*in11.*in13.*2.0;
end
end

