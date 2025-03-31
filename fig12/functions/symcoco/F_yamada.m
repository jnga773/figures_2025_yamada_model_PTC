function varargout=F_yamada(action,varargin)
%% Automatically generated with matlabFunction
%#ok<*DEFNU,*INUSD,*INUSL>

switch action
  case 'nargs'
   varargout{1}=2;
   return
  case 'nout'
   varargout{1}=3;
   return
  case 'argrange'
   varargout{1}=struct('x',1:3,'p',4:7);
   return
  case 'argsize'
   varargout{1}=struct('x',3,'p',4);
   return
  case 'vector'
   varargout{1}=struct('x',1,'p',1);
   return
  case 'extension'
   varargout{1}='rhs';
   return
  case 'maxorder'
   varargout{1}=2;
   return
end
nout=3;
order=varargin{1};
f=str2func(sprintf('F_yamada_%s_%d',action,order));
varargout=cell(nout,1);
[varargout{:}]=f(varargin{2:end});
end



function [out1,out2,out3] = F_yamada_rhs_0(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12,in13,in14)
%F_yamada_rhs_0
%    [OUT1,OUT2,OUT3] = F_yamada_rhs_0(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12,IN13,IN14)

%    This function was generated by the Symbolic Math Toolbox version 25.1.
%    31-Mar-2025 15:55:59

out1 = -in4.*(in1-in5+in1.*in3);
if nargout > 1
    out2 = -in4.*(in2-in6+in2.*in3.*in7);
end
if nargout > 2
    out3 = -in3.*(-in1+in2+1.0);
end
end


function [out1,out2,out3] = F_yamada_rhs_1(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12,in13,in14)
%F_yamada_rhs_1
%    [OUT1,OUT2,OUT3] = F_yamada_rhs_1(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12,IN13,IN14)

%    This function was generated by the Symbolic Math Toolbox version 25.1.
%    31-Mar-2025 15:55:59

out1 = -in11.*(in1-in5+in1.*in3)-in4.*(in8-in12+in1.*in10+in3.*in8);
if nargout > 1
    out2 = -in4.*(in9-in13+in2.*in3.*in14+in2.*in7.*in10+in3.*in7.*in9)-in11.*(in2-in6+in2.*in3.*in7);
end
if nargout > 2
    out3 = -in10.*(-in1+in2+1.0)+in3.*(in8-in9);
end
end


function [out1,out2,out3] = F_yamada_rhs_2(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12,in13,in14)
%F_yamada_rhs_2
%    [OUT1,OUT2,OUT3] = F_yamada_rhs_2(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12,IN13,IN14)

%    This function was generated by the Symbolic Math Toolbox version 25.1.
%    31-Mar-2025 15:55:59

out1 = in11.*(in8-in12+in1.*in10+in3.*in8).*-2.0-in4.*in8.*in10.*2.0;
if nargout > 1
    out2 = in11.*(in9-in13+in2.*in3.*in14+in2.*in7.*in10+in3.*in7.*in9).*-2.0-in4.*(in2.*in10.*in14.*2.0+in3.*in9.*in14.*2.0+in7.*in9.*in10.*2.0);
end
if nargout > 2
    out3 = in10.*(in8-in9).*2.0;
end
end

