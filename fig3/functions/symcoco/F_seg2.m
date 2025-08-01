function varargout=F_seg2(action,varargin)
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
   varargout{1}=struct('x',1:6,'p',7:18);
   return
  case 'argsize'
   varargout{1}=struct('x',6,'p',12);
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
f=str2func(sprintf('F_seg2_%s_%d',action,order));
varargout=cell(nout,1);
[varargout{:}]=f(varargin{2:end});
end



function [out1,out2,out3,out4,out5,out6] = F_seg2_rhs_0(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12,in13,in14,in15,in16,in17,in18,in19,in20,in21,in22,in23,in24,in25,in26,in27,in28,in29,in30,in31,in32,in33,in34,in35,in36)
%F_seg2_rhs_0
%    [OUT1,OUT2,OUT3,OUT4,OUT5,OUT6] = F_seg2_rhs_0(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12,IN13,IN14,IN15,IN16,IN17,IN18,IN19,IN20,IN21,IN22,IN23,IN24,IN25,IN26,IN27,IN28,IN29,IN30,IN31,IN32,IN33,IN34,IN35,IN36)

%    This function was generated by the Symbolic Math Toolbox version 25.1.
%    26-Jul-2025 15:36:32

t2 = -in1;
t3 = in13-1.0;
out1 = in7.*t3.*(in1-in8+in1.*in3);
if nargout > 1
    out2 = in7.*t3.*(in2-in9+in2.*in3.*in10);
end
if nargout > 2
    t4 = in2+t2+1.0;
    out3 = in3.*t3.*t4;
end
if nargout > 3
    t5 = in3.*in6.*t3;
    out4 = t5-in4.*in7.*t3.*(in3+1.0);
end
if nargout > 4
    out5 = -t5-in5.*in7.*t3.*(in3.*in10+1.0);
end
if nargout > 5
    out6 = -in6.*t3.*t4+in4.*in7.*t2.*t3-in2.*in5.*in7.*in10.*t3;
end
end


function [out1,out2,out3,out4,out5,out6] = F_seg2_rhs_1(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12,in13,in14,in15,in16,in17,in18,in19,in20,in21,in22,in23,in24,in25,in26,in27,in28,in29,in30,in31,in32,in33,in34,in35,in36)
%F_seg2_rhs_1
%    [OUT1,OUT2,OUT3,OUT4,OUT5,OUT6] = F_seg2_rhs_1(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12,IN13,IN14,IN15,IN16,IN17,IN18,IN19,IN20,IN21,IN22,IN23,IN24,IN25,IN26,IN27,IN28,IN29,IN30,IN31,IN32,IN33,IN34,IN35,IN36)

%    This function was generated by the Symbolic Math Toolbox version 25.1.
%    26-Jul-2025 15:36:32

t2 = in1.*in3;
t3 = in3.*in10;
t4 = in3+1.0;
t6 = in3.*in6.*in31;
t7 = -in1;
t8 = -in8;
t9 = -in9;
t10 = -in20;
t11 = in13-1.0;
t5 = in2.*t3;
t12 = t3+1.0;
t13 = in19+t10;
t14 = in2+t7+1.0;
t15 = in3.*in24.*t11;
t16 = in6.*in21.*t11;
t17 = in1+t2+t8;
out1 = in7.*t11.*(in19-in26+in1.*in21+in3.*in19)+in7.*in31.*t17+in25.*t11.*t17;
if nargout > 1
    t18 = in2+t5+t9;
    out2 = in7.*t11.*(in20-in27+in20.*t3+in2.*in3.*in28+in2.*in10.*in21)+in7.*in31.*t18+in25.*t11.*t18;
end
if nargout > 2
    out3 = in3.*in31.*t14-in3.*t11.*t13+in21.*t11.*t14;
end
if nargout > 3
    out4 = t6+t15+t16-in4.*in7.*in21.*t11-in4.*in7.*in31.*t4-in4.*in25.*t4.*t11-in7.*in22.*t4.*t11;
end
if nargout > 4
    out5 = -t6-t15-t16-in5.*in7.*in31.*t12-in5.*in25.*t11.*t12-in7.*in23.*t11.*t12-in5.*in7.*t11.*(in3.*in28+in10.*in21);
end
if nargout > 5
    out6 = -in6.*in31.*t14+in6.*t11.*t13-in24.*t11.*t14-in4.*in7.*in19.*t11+in4.*in7.*in31.*t7+in4.*in25.*t7.*t11+in7.*in22.*t7.*t11-in2.*in5.*in7.*in10.*in31-in2.*in5.*in7.*in28.*t11-in2.*in5.*in10.*in25.*t11-in2.*in7.*in10.*in23.*t11+in5.*in7.*in10.*t10.*t11;
end
end


function [out1,out2,out3,out4,out5,out6] = F_seg2_rhs_2(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12,in13,in14,in15,in16,in17,in18,in19,in20,in21,in22,in23,in24,in25,in26,in27,in28,in29,in30,in31,in32,in33,in34,in35,in36)
%F_seg2_rhs_2
%    [OUT1,OUT2,OUT3,OUT4,OUT5,OUT6] = F_seg2_rhs_2(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12,IN13,IN14,IN15,IN16,IN17,IN18,IN19,IN20,IN21,IN22,IN23,IN24,IN25,IN26,IN27,IN28,IN29,IN30,IN31,IN32,IN33,IN34,IN35,IN36)

%    This function was generated by the Symbolic Math Toolbox version 25.1.
%    26-Jul-2025 15:36:32

t2 = in3.*in10;
t3 = in1.*in21;
t4 = in3.*in19;
t5 = in3.*in28;
t6 = in10.*in21;
t7 = in3+1.0;
t11 = -in1;
t12 = -in20;
t13 = -in26;
t14 = -in27;
t15 = in13-1.0;
t17 = in3.*in24.*in31.*2.0;
t18 = in6.*in21.*in31.*2.0;
t8 = in2.*t5;
t9 = in2.*t6;
t10 = in20.*t2;
t16 = t2+1.0;
t19 = in19+t12;
t20 = in2+t11+1.0;
t21 = t5+t6;
t22 = in21.*in24.*t15.*2.0;
t23 = in19+t3+t4+t13;
out1 = in7.*in31.*t23.*2.0+in25.*t15.*t23.*2.0+in25.*in31.*(in1-in8+in1.*in3).*2.0+in7.*in19.*in21.*t15.*2.0;
if nargout > 1
    t24 = in20+t8+t9+t10+t14;
    out2 = in25.*in31.*(in2-in9+in2.*t2).*2.0+in7.*in31.*t24.*2.0+in25.*t15.*t24.*2.0+in7.*t15.*(in20.*t5.*2.0+in20.*t6.*2.0+in2.*in21.*in28.*2.0);
end
if nargout > 2
    out3 = in3.*in31.*t19.*-2.0+in21.*in31.*t20.*2.0-in21.*t15.*t19.*2.0;
end
if nargout > 3
    out4 = t17+t18+t22-in4.*in7.*in21.*in31.*2.0-in4.*in21.*in25.*t15.*2.0-in7.*in21.*in22.*t15.*2.0-in4.*in25.*in31.*t7.*2.0-in7.*in22.*in31.*t7.*2.0-in22.*in25.*t7.*t15.*2.0;
end
if nargout > 4
    out5 = -t17-t18-t22-in5.*in7.*in31.*t21.*2.0-in5.*in25.*in31.*t16.*2.0-in7.*in23.*in31.*t16.*2.0-in5.*in25.*t15.*t21.*2.0-in7.*in23.*t15.*t21.*2.0-in23.*in25.*t15.*t16.*2.0-in5.*in7.*in21.*in28.*t15.*2.0;
end
if nargout > 5
    out6 = in6.*in31.*t19.*2.0-in24.*in31.*t20.*2.0+in24.*t15.*t19.*2.0-in1.*in4.*in25.*in31.*2.0-in1.*in7.*in22.*in31.*2.0-in4.*in7.*in19.*in31.*2.0-in1.*in22.*in25.*t15.*2.0-in4.*in19.*in25.*t15.*2.0-in7.*in19.*in22.*t15.*2.0-in2.*in5.*in7.*in28.*in31.*2.0-in2.*in5.*in10.*in25.*in31.*2.0-in2.*in7.*in10.*in23.*in31.*2.0-in5.*in7.*in10.*in20.*in31.*2.0-in2.*in5.*in25.*in28.*t15.*2.0-in2.*in7.*in23.*in28.*t15.*2.0-in2.*in10.*in23.*in25.*t15.*2.0-in5.*in7.*in20.*in28.*t15.*2.0-in5.*in10.*in20.*in25.*t15.*2.0-in7.*in10.*in20.*in23.*t15.*2.0;
end
end

