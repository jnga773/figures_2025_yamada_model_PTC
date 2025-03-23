function varargout=F_bcs_VAR(action,varargin)
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
   varargout{1}=struct('u',1:9);
   return
  case 'argsize'
   varargout{1}=struct('u',9);
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
f=str2func(sprintf('F_bcs_VAR_%s_%d',action,order));
varargout=cell(nout,1);
[varargout{:}]=f(varargin{2:end});
end



function [out1,out2,out3,out4] = F_bcs_VAR_rhs_0(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12,in13,in14,in15,in16,in17,in18)
%F_bcs_VAR_rhs_0
%    [OUT1,OUT2,OUT3,OUT4] = F_bcs_VAR_rhs_0(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12,IN13,IN14,IN15,IN16,IN17,IN18)

%    This function was generated by the Symbolic Math Toolbox version 24.2.
%    24-Mar-2025 12:23:19

out1 = in4-in1.*in7;
if nargout > 1
    out2 = in5-in2.*in7;
end
if nargout > 2
    out3 = in6-in3.*in7;
end
if nargout > 3
    out4 = -in8+in1.*conj(in1)+in2.*conj(in2)+in3.*conj(in3);
end
end


function [out1,out2,out3,out4] = F_bcs_VAR_rhs_1(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12,in13,in14,in15,in16,in17,in18)
%F_bcs_VAR_rhs_1
%    [OUT1,OUT2,OUT3,OUT4] = F_bcs_VAR_rhs_1(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12,IN13,IN14,IN15,IN16,IN17,IN18)

%    This function was generated by the Symbolic Math Toolbox version 24.2.
%    24-Mar-2025 12:23:19

out1 = in13-in1.*in16-in7.*in10;
if nargout > 1
    out2 = in14-in2.*in16-in7.*in11;
end
if nargout > 2
    out3 = in15-in3.*in16-in7.*in12;
end
if nargout > 3
    out4 = -in17+in1.*conj(in10)+in10.*conj(in1)+in2.*conj(in11)+in11.*conj(in2)+in3.*conj(in12)+in12.*conj(in3);
end
end


function [out1,out2,out3,out4] = F_bcs_VAR_rhs_2(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12,in13,in14,in15,in16,in17,in18)
%F_bcs_VAR_rhs_2
%    [OUT1,OUT2,OUT3,OUT4] = F_bcs_VAR_rhs_2(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12,IN13,IN14,IN15,IN16,IN17,IN18)

%    This function was generated by the Symbolic Math Toolbox version 24.2.
%    24-Mar-2025 12:23:19

out1 = in10.*in16.*-2.0;
if nargout > 1
    out2 = in11.*in16.*-2.0;
end
if nargout > 2
    out3 = in12.*in16.*-2.0;
end
if nargout > 3
    out4 = in10.*conj(in10).*2.0+in11.*conj(in11).*2.0+in12.*conj(in12).*2.0;
end
end

