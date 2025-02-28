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
   varargout{1}=struct('x',1:6,'p',7:19);
   return
  case 'argsize'
   varargout{1}=struct('x',6,'p',13);
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



function [out1,out2,out3,out4,out5,out6] = F_seg2_rhs_0(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12,in13,in14,in15,in16,in17,in18,in19,in20,in21,in22,in23,in24,in25,in26,in27,in28,in29,in30,in31,in32,in33,in34,in35,in36,in37,in38)
%F_seg2_rhs_0
%    [OUT1,OUT2,OUT3,OUT4,OUT5,OUT6] = F_seg2_rhs_0(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12,IN13,IN14,IN15,IN16,IN17,IN18,IN19,IN20,IN21,IN22,IN23,IN24,IN25,IN26,IN27,IN28,IN29,IN30,IN31,IN32,IN33,IN34,IN35,IN36,IN37,IN38)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    28-Feb-2025 11:28:35

t2 = -in1;
t3 = in14-1.0;
out1 = in7.*in11.*t3.*(in1-in8+in1.*in3);
if nargout > 1
    out2 = in7.*in11.*t3.*(in2-in9+in2.*in3.*in10);
end
if nargout > 2
    t4 = in2+t2+1.0;
    out3 = in3.*in11.*t3.*t4;
end
if nargout > 3
    t5 = in3.*in6.*in11.*t3;
    out4 = t5-in4.*in7.*in11.*t3.*(in3+1.0);
end
if nargout > 4
    out5 = -t5-in5.*in7.*in11.*t3.*(in3.*in10+1.0);
end
if nargout > 5
    out6 = -in6.*in11.*t3.*t4+in4.*in7.*in11.*t2.*t3-in2.*in5.*in7.*in10.*in11.*t3;
end
end


function [out1,out2,out3,out4,out5,out6] = F_seg2_rhs_1(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12,in13,in14,in15,in16,in17,in18,in19,in20,in21,in22,in23,in24,in25,in26,in27,in28,in29,in30,in31,in32,in33,in34,in35,in36,in37,in38)
%F_seg2_rhs_1
%    [OUT1,OUT2,OUT3,OUT4,OUT5,OUT6] = F_seg2_rhs_1(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12,IN13,IN14,IN15,IN16,IN17,IN18,IN19,IN20,IN21,IN22,IN23,IN24,IN25,IN26,IN27,IN28,IN29,IN30,IN31,IN32,IN33,IN34,IN35,IN36,IN37,IN38)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    28-Feb-2025 11:28:35

t2 = in1.*in3;
t3 = in3.*in10;
t4 = in3+1.0;
t6 = -in1;
t7 = -in8;
t8 = -in9;
t9 = -in21;
t10 = in14-1.0;
t11 = in3.*in6.*in11.*in33;
t5 = in2.*t3;
t12 = t3+1.0;
t13 = in20+t9;
t14 = in2+t6+1.0;
t15 = in3.*in6.*in30.*t10;
t16 = in3.*in11.*in25.*t10;
t17 = in6.*in11.*in22.*t10;
t18 = in1+t2+t7;
out1 = in7.*in11.*in33.*t18+in7.*in30.*t10.*t18+in11.*in26.*t10.*t18+in7.*in11.*t10.*(in20-in27+in1.*in22+in3.*in20);
if nargout > 1
    t19 = in2+t5+t8;
    out2 = in7.*in11.*in33.*t19+in7.*in30.*t10.*t19+in11.*in26.*t10.*t19+in7.*in11.*t10.*(in21-in28+in21.*t3+in2.*in3.*in29+in2.*in10.*in22);
end
if nargout > 2
    out3 = in3.*in11.*in33.*t14-in3.*in11.*t10.*t13+in3.*in30.*t10.*t14+in11.*in22.*t10.*t14;
end
if nargout > 3
    out4 = t11+t15+t16+t17-in4.*in7.*in11.*in22.*t10-in4.*in7.*in11.*in33.*t4-in4.*in7.*in30.*t4.*t10-in4.*in11.*in26.*t4.*t10-in7.*in11.*in23.*t4.*t10;
end
if nargout > 4
    out5 = -t11-t15-t16-t17-in5.*in7.*in11.*in33.*t12-in5.*in7.*in30.*t10.*t12-in5.*in11.*in26.*t10.*t12-in7.*in11.*in24.*t10.*t12-in5.*in7.*in11.*t10.*(in3.*in29+in10.*in22);
end
if nargout > 5
    out6 = -in6.*in11.*in33.*t14+in6.*in11.*t10.*t13-in6.*in30.*t10.*t14-in11.*in25.*t10.*t14-in4.*in7.*in11.*in20.*t10+in4.*in7.*in11.*in33.*t6+in4.*in7.*in30.*t6.*t10+in4.*in11.*in26.*t6.*t10+in7.*in11.*in23.*t6.*t10-in2.*in5.*in7.*in10.*in11.*in33-in2.*in5.*in7.*in10.*in30.*t10-in2.*in5.*in7.*in11.*in29.*t10-in2.*in5.*in10.*in11.*in26.*t10-in2.*in7.*in10.*in11.*in24.*t10+in5.*in7.*in10.*in11.*t9.*t10;
end
end


function [out1,out2,out3,out4,out5,out6] = F_seg2_rhs_2(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12,in13,in14,in15,in16,in17,in18,in19,in20,in21,in22,in23,in24,in25,in26,in27,in28,in29,in30,in31,in32,in33,in34,in35,in36,in37,in38)
%F_seg2_rhs_2
%    [OUT1,OUT2,OUT3,OUT4,OUT5,OUT6] = F_seg2_rhs_2(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12,IN13,IN14,IN15,IN16,IN17,IN18,IN19,IN20,IN21,IN22,IN23,IN24,IN25,IN26,IN27,IN28,IN29,IN30,IN31,IN32,IN33,IN34,IN35,IN36,IN37,IN38)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    28-Feb-2025 11:28:35

t2 = in1.*in3;
t3 = in3.*in10;
t4 = in1.*in22;
t5 = in3.*in20;
t6 = in3.*in29;
t7 = in10.*in22;
t8 = in3+1.0;
t13 = -in1;
t14 = -in8;
t15 = -in9;
t16 = -in21;
t17 = -in27;
t18 = -in28;
t19 = in14-1.0;
t21 = in3.*in6.*in30.*in33.*2.0;
t22 = in3.*in11.*in25.*in33.*2.0;
t23 = in6.*in11.*in22.*in33.*2.0;
t9 = in2.*t3;
t10 = in2.*t6;
t11 = in2.*t7;
t12 = in21.*t3;
t20 = t3+1.0;
t24 = in20+t16;
t25 = in2+t13+1.0;
t26 = t6+t7;
t27 = in3.*in25.*in30.*t19.*2.0;
t28 = in6.*in22.*in30.*t19.*2.0;
t29 = in11.*in22.*in25.*t19.*2.0;
t30 = in1+t2+t14;
t32 = in20+t4+t5+t17;
out1 = in7.*in11.*in33.*t32.*2.0+in7.*in30.*in33.*t30.*2.0+in11.*in26.*in33.*t30.*2.0+in7.*in30.*t19.*t32.*2.0+in11.*in26.*t19.*t32.*2.0+in26.*in30.*t19.*t30.*2.0+in7.*in11.*in20.*in22.*t19.*2.0;
if nargout > 1
    t31 = in2+t9+t15;
    t33 = in21+t10+t11+t12+t18;
    out2 = in7.*in11.*in33.*t33.*2.0+in7.*in30.*in33.*t31.*2.0+in11.*in26.*in33.*t31.*2.0+in7.*in30.*t19.*t33.*2.0+in11.*in26.*t19.*t33.*2.0+in26.*in30.*t19.*t31.*2.0+in7.*in11.*t19.*(in21.*t6.*2.0+in21.*t7.*2.0+in2.*in22.*in29.*2.0);
end
if nargout > 2
    out3 = in3.*in11.*in33.*t24.*-2.0+in3.*in30.*in33.*t25.*2.0+in11.*in22.*in33.*t25.*2.0-in3.*in30.*t19.*t24.*2.0-in11.*in22.*t19.*t24.*2.0+in22.*in30.*t19.*t25.*2.0;
end
if nargout > 3
    out4 = t21+t22+t23+t27+t28+t29-in4.*in7.*in11.*in22.*in33.*2.0-in4.*in7.*in22.*in30.*t19.*2.0-in4.*in7.*in30.*in33.*t8.*2.0-in4.*in11.*in22.*in26.*t19.*2.0-in4.*in11.*in26.*in33.*t8.*2.0-in7.*in11.*in22.*in23.*t19.*2.0-in7.*in11.*in23.*in33.*t8.*2.0-in4.*in26.*in30.*t8.*t19.*2.0-in7.*in23.*in30.*t8.*t19.*2.0-in11.*in23.*in26.*t8.*t19.*2.0;
end
if nargout > 4
    out5 = -t21-t22-t23-t27-t28-t29-in5.*in7.*in11.*in33.*t26.*2.0-in5.*in7.*in30.*in33.*t20.*2.0-in5.*in11.*in26.*in33.*t20.*2.0-in7.*in11.*in24.*in33.*t20.*2.0-in5.*in7.*in30.*t19.*t26.*2.0-in5.*in11.*in26.*t19.*t26.*2.0-in7.*in11.*in24.*t19.*t26.*2.0-in5.*in26.*in30.*t19.*t20.*2.0-in7.*in24.*in30.*t19.*t20.*2.0-in11.*in24.*in26.*t19.*t20.*2.0-in5.*in7.*in11.*in22.*in29.*t19.*2.0;
end
if nargout > 5
    et1 = in6.*in11.*in33.*t24.*2.0-in6.*in30.*in33.*t25.*2.0-in11.*in25.*in33.*t25.*2.0+in6.*in30.*t19.*t24.*2.0+in11.*in25.*t19.*t24.*2.0-in25.*in30.*t19.*t25.*2.0-in1.*in4.*in7.*in30.*in33.*2.0-in1.*in4.*in11.*in26.*in33.*2.0-in1.*in7.*in11.*in23.*in33.*2.0-in4.*in7.*in11.*in20.*in33.*2.0-in1.*in4.*in26.*in30.*t19.*2.0-in1.*in7.*in23.*in30.*t19.*2.0-in1.*in11.*in23.*in26.*t19.*2.0-in4.*in7.*in20.*in30.*t19.*2.0-in4.*in11.*in20.*in26.*t19.*2.0-in7.*in11.*in20.*in23.*t19.*2.0-in2.*in5.*in7.*in10.*in30.*in33.*2.0-in2.*in5.*in7.*in11.*in29.*in33.*2.0-in2.*in5.*in10.*in11.*in26.*in33.*2.0-in2.*in7.*in10.*in11.*in24.*in33.*2.0-in5.*in7.*in10.*in11.*in21.*in33.*2.0-in2.*in5.*in7.*in29.*in30.*t19.*2.0-in2.*in5.*in10.*in26.*in30.*t19.*2.0-in2.*in5.*in11.*in26.*in29.*t19.*2.0-in2.*in7.*in10.*in24.*in30.*t19.*2.0;
    et2 = in2.*in7.*in11.*in24.*in29.*t19.*-2.0-in2.*in10.*in11.*in24.*in26.*t19.*2.0-in5.*in7.*in10.*in21.*in30.*t19.*2.0-in5.*in7.*in11.*in21.*in29.*t19.*2.0-in5.*in10.*in11.*in21.*in26.*t19.*2.0-in7.*in10.*in11.*in21.*in24.*t19.*2.0;
    out6 = et1+et2;
end
end

