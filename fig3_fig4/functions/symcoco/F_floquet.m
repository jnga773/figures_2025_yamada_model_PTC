function varargout=F_floquet(action,varargin)
%% Automatically generated with matlabFunction
%#ok<*DEFNU,*INUSD,*INUSL>

switch action
  case 'nargs'
   varargout{1}=2;
   return
  case 'nout'
   varargout{1}=6;
   return
  case 'argrange'
   varargout{1}=struct('x',1:6,'p',7:13);
   return
  case 'argsize'
   varargout{1}=struct('x',6,'p',7);
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
nout=6;
order=varargin{1};
f=str2func(sprintf('F_floquet_%s_%d',action,order));
varargout=cell(nout,1);
[varargout{:}]=f(varargin{2:end});
end



function [out1,out2,out3,out4,out5,out6] = F_floquet_rhs_0(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12,in13,in14,in15,in16,in17,in18,in19,in20,in21,in22,in23,in24,in25,in26)
%F_floquet_rhs_0
%    [OUT1,OUT2,OUT3,OUT4,OUT5,OUT6] = F_floquet_rhs_0(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12,IN13,IN14,IN15,IN16,IN17,IN18,IN19,IN20,IN21,IN22,IN23,IN24,IN25,IN26)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    12-Feb-2025 17:14:03

out1 = -in7.*in13.*(in1-in8+in1.*in3);
if nargout > 1
    out2 = -in7.*in13.*(in2-in9+in2.*in3.*in10);
end
if nargout > 2
    t2 = in3.*in6.*in13;
    t3 = -in1;
    t4 = in2+t3+1.0;
    out3 = -in3.*in13.*t4;
end
if nargout > 3
    out4 = -t2+in4.*in7.*in13.*(in3+1.0);
end
if nargout > 4
    out5 = t2+in5.*in7.*in13.*(in3.*in10+1.0);
end
if nargout > 5
    out6 = in6.*in13.*t4+in1.*in4.*in7.*in13+in2.*in5.*in7.*in10.*in13;
end
end


function [out1,out2,out3,out4,out5,out6] = F_floquet_rhs_1(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12,in13,in14,in15,in16,in17,in18,in19,in20,in21,in22,in23,in24,in25,in26)
%F_floquet_rhs_1
%    [OUT1,OUT2,OUT3,OUT4,OUT5,OUT6] = F_floquet_rhs_1(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12,IN13,IN14,IN15,IN16,IN17,IN18,IN19,IN20,IN21,IN22,IN23,IN24,IN25,IN26)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    12-Feb-2025 17:14:03

t2 = in1.*in3;
t3 = in3.*in10;
t4 = in3+1.0;
t6 = in3.*in6.*in26;
t7 = in3.*in13.*in19;
t8 = in6.*in13.*in16;
t9 = -in1;
t10 = -in8;
t11 = -in9;
t12 = -in15;
t5 = in2.*t3;
t13 = t3+1.0;
t14 = in14+t12;
t15 = in2+t9+1.0;
t16 = in1+t2+t10;
out1 = -in7.*in13.*(in14-in21+in1.*in16+in3.*in14)-in7.*in26.*t16-in13.*in20.*t16;
if nargout > 1
    t17 = in2+t5+t11;
    out2 = -in7.*in13.*(in15-in22+in15.*t3+in2.*in3.*in23+in2.*in10.*in16)-in7.*in26.*t17-in13.*in20.*t17;
end
if nargout > 2
    out3 = in3.*in13.*t14-in3.*in26.*t15-in13.*in16.*t15;
end
if nargout > 3
    out4 = -t6-t7-t8+in4.*in7.*in13.*in16+in4.*in7.*in26.*t4+in4.*in13.*in20.*t4+in7.*in13.*in17.*t4;
end
if nargout > 4
    out5 = t6+t7+t8+in5.*in7.*in26.*t13+in5.*in13.*in20.*t13+in7.*in13.*in18.*t13+in5.*in7.*in13.*(in3.*in23+in10.*in16);
end
if nargout > 5
    out6 = -in6.*in13.*t14+in6.*in26.*t15+in13.*in19.*t15+in1.*in4.*in7.*in26+in1.*in4.*in13.*in20+in1.*in7.*in13.*in17+in4.*in7.*in13.*in14+in2.*in5.*in7.*in10.*in26+in2.*in5.*in7.*in13.*in23+in2.*in5.*in10.*in13.*in20+in2.*in7.*in10.*in13.*in18+in5.*in7.*in10.*in13.*in15;
end
end


function [out1,out2,out3,out4,out5,out6] = F_floquet_rhs_2(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12,in13,in14,in15,in16,in17,in18,in19,in20,in21,in22,in23,in24,in25,in26)
%F_floquet_rhs_2
%    [OUT1,OUT2,OUT3,OUT4,OUT5,OUT6] = F_floquet_rhs_2(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12,IN13,IN14,IN15,IN16,IN17,IN18,IN19,IN20,IN21,IN22,IN23,IN24,IN25,IN26)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    12-Feb-2025 17:14:03

t2 = in3.*in10;
t3 = in1.*in16;
t4 = in3.*in14;
t5 = in3.*in23;
t6 = in10.*in16;
t7 = in3+1.0;
t11 = -in1;
t12 = -in15;
t13 = -in21;
t14 = -in22;
t16 = in3.*in19.*in26.*2.0;
t17 = in6.*in16.*in26.*2.0;
t18 = in13.*in16.*in19.*2.0;
t8 = in2.*t5;
t9 = in2.*t6;
t10 = in15.*t2;
t15 = t2+1.0;
t19 = in14+t12;
t20 = in2+t11+1.0;
t21 = t5+t6;
t22 = in14+t3+t4+t13;
out1 = in7.*in26.*t22.*-2.0-in13.*in20.*t22.*2.0-in20.*in26.*(in1-in8+in1.*in3).*2.0-in7.*in13.*in14.*in16.*2.0;
if nargout > 1
    t23 = in15+t8+t9+t10+t14;
    out2 = in20.*in26.*(in2-in9+in2.*t2).*-2.0-in7.*in26.*t23.*2.0-in13.*in20.*t23.*2.0-in7.*in13.*(in15.*t5.*2.0+in15.*t6.*2.0+in2.*in16.*in23.*2.0);
end
if nargout > 2
    out3 = in3.*in26.*t19.*2.0+in13.*in16.*t19.*2.0-in16.*in26.*t20.*2.0;
end
if nargout > 3
    out4 = -t16-t17-t18+in4.*in7.*in16.*in26.*2.0+in4.*in13.*in16.*in20.*2.0+in7.*in13.*in16.*in17.*2.0+in4.*in20.*in26.*t7.*2.0+in7.*in17.*in26.*t7.*2.0+in13.*in17.*in20.*t7.*2.0;
end
if nargout > 4
    out5 = t16+t17+t18+in5.*in7.*in26.*t21.*2.0+in5.*in13.*in20.*t21.*2.0+in7.*in13.*in18.*t21.*2.0+in5.*in20.*in26.*t15.*2.0+in7.*in18.*in26.*t15.*2.0+in13.*in18.*in20.*t15.*2.0+in5.*in7.*in13.*in16.*in23.*2.0;
end
if nargout > 5
    out6 = in6.*in26.*t19.*-2.0-in13.*in19.*t19.*2.0+in19.*in26.*t20.*2.0+in1.*in4.*in20.*in26.*2.0+in1.*in7.*in17.*in26.*2.0+in1.*in13.*in17.*in20.*2.0+in4.*in7.*in14.*in26.*2.0+in4.*in13.*in14.*in20.*2.0+in7.*in13.*in14.*in17.*2.0+in2.*in5.*in7.*in23.*in26.*2.0+in2.*in5.*in10.*in20.*in26.*2.0+in2.*in5.*in13.*in20.*in23.*2.0+in2.*in7.*in10.*in18.*in26.*2.0+in2.*in7.*in13.*in18.*in23.*2.0+in2.*in10.*in13.*in18.*in20.*2.0+in5.*in7.*in10.*in15.*in26.*2.0+in5.*in7.*in13.*in15.*in23.*2.0+in5.*in10.*in13.*in15.*in20.*2.0+in7.*in10.*in13.*in15.*in18.*2.0;
end
end

