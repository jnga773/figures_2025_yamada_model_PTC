function varargout=F_seg4(action,varargin)
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
   varargout{1}=struct('x',1:3,'p',4:16);
   return
  case 'argsize'
   varargout{1}=struct('x',3,'p',13);
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
f=str2func(sprintf('F_seg4_%s_%d',action,order));
varargout=cell(nout,1);
[varargout{:}]=f(varargin{2:end});
end



function [out1,out2,out3] = F_seg4_rhs_0(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12,in13,in14,in15,in16,in17,in18,in19,in20,in21,in22,in23,in24,in25,in26,in27,in28,in29,in30,in31,in32)
%F_seg4_rhs_0
%    [OUT1,OUT2,OUT3] = F_seg4_rhs_0(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12,IN13,IN14,IN15,IN16,IN17,IN18,IN19,IN20,IN21,IN22,IN23,IN24,IN25,IN26,IN27,IN28,IN29,IN30,IN31,IN32)

%    This function was generated by the Symbolic Math Toolbox version 25.1.
%    24-Jun-2025 15:16:30

out1 = -in4.*in8.*in9.*(in1-in5+in1.*in3);
if nargout > 1
    out2 = -in4.*in8.*in9.*(in2-in6+in2.*in3.*in7);
end
if nargout > 2
    out3 = -in3.*in8.*in9.*(-in1+in2+1.0);
end
end


function [out1,out2,out3] = F_seg4_rhs_1(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12,in13,in14,in15,in16,in17,in18,in19,in20,in21,in22,in23,in24,in25,in26,in27,in28,in29,in30,in31,in32)
%F_seg4_rhs_1
%    [OUT1,OUT2,OUT3] = F_seg4_rhs_1(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12,IN13,IN14,IN15,IN16,IN17,IN18,IN19,IN20,IN21,IN22,IN23,IN24,IN25,IN26,IN27,IN28,IN29,IN30,IN31,IN32)

%    This function was generated by the Symbolic Math Toolbox version 25.1.
%    24-Jun-2025 15:16:31

t2 = in1.*in3;
t3 = in2.*in3.*in7;
t4 = -in1;
t5 = -in5;
t6 = -in6;
t7 = in2+t4+1.0;
t8 = in1+t2+t5;
out1 = -in4.*in8.*in25.*t8-in4.*in9.*in24.*t8-in8.*in9.*in20.*t8-in4.*in8.*in9.*(in17-in21+in1.*in19+in3.*in17);
if nargout > 1
    t9 = in2+t3+t6;
    out2 = -in4.*in8.*in9.*(in18-in22+in2.*in3.*in23+in2.*in7.*in19+in3.*in7.*in18)-in4.*in8.*in25.*t9-in4.*in9.*in24.*t9-in8.*in9.*in20.*t9;
end
if nargout > 2
    out3 = -in3.*in8.*in25.*t7-in3.*in9.*in24.*t7-in8.*in9.*in19.*t7+in3.*in8.*in9.*(in17-in18);
end
end


function [out1,out2,out3] = F_seg4_rhs_2(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12,in13,in14,in15,in16,in17,in18,in19,in20,in21,in22,in23,in24,in25,in26,in27,in28,in29,in30,in31,in32)
%F_seg4_rhs_2
%    [OUT1,OUT2,OUT3] = F_seg4_rhs_2(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12,IN13,IN14,IN15,IN16,IN17,IN18,IN19,IN20,IN21,IN22,IN23,IN24,IN25,IN26,IN27,IN28,IN29,IN30,IN31,IN32)

%    This function was generated by the Symbolic Math Toolbox version 25.1.
%    24-Jun-2025 15:16:31

t2 = in1.*in3;
t3 = in1.*in19;
t4 = in3.*in17;
t5 = in2.*in3.*in7;
t6 = in2.*in3.*in23;
t7 = in2.*in7.*in19;
t8 = in3.*in7.*in18;
t9 = -in1;
t10 = -in5;
t11 = -in6;
t12 = -in18;
t13 = -in21;
t14 = -in22;
t15 = in17+t12;
t16 = in2+t9+1.0;
t17 = in1+t2+t10;
t18 = in2+t5+t11;
t19 = in17+t3+t4+t13;
out1 = in4.*in8.*in25.*t19.*-2.0-in4.*in9.*in24.*t19.*2.0-in8.*in9.*in20.*t19.*2.0-in4.*in24.*in25.*t17.*2.0-in8.*in20.*in25.*t17.*2.0-in9.*in20.*in24.*t17.*2.0-in4.*in8.*in9.*in17.*in19.*2.0;
if nargout > 1
    t20 = in18+t6+t7+t8+t14;
    out2 = in4.*in8.*in25.*t20.*-2.0-in4.*in9.*in24.*t20.*2.0-in8.*in9.*in20.*t20.*2.0-in4.*in24.*in25.*t18.*2.0-in8.*in20.*in25.*t18.*2.0-in9.*in20.*in24.*t18.*2.0-in4.*in8.*in9.*(in2.*in19.*in23.*2.0+in3.*in18.*in23.*2.0+in7.*in18.*in19.*2.0);
end
if nargout > 2
    out3 = in3.*in8.*in25.*t15.*2.0+in3.*in9.*in24.*t15.*2.0+in8.*in9.*in19.*t15.*2.0-in3.*in24.*in25.*t16.*2.0-in8.*in19.*in25.*t16.*2.0-in9.*in19.*in24.*t16.*2.0;
end
end

