function varargout=F_bcs_PR_segs(action,varargin)
%% Automatically generated with matlabFunction
%#ok<*DEFNU,*INUSD,*INUSL>

switch action
  case 'nargs'
   varargout{1}=1;
   return
  case 'nout'
   varargout{1}=22;
   return
  case 'argrange'
   varargout{1}=struct('u',1:49);
   return
  case 'argsize'
   varargout{1}=struct('u',49);
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
nout=22;
order=varargin{1};
f=str2func(sprintf('F_bcs_PR_segs_%s_%d',action,order));
varargout=cell(nout,1);
[varargout{:}]=f(varargin{2:end});
end



function [out1,out2,out3,out4,out5,out6,out7,out8,out9,out10,out11,out12,out13,out14,out15,out16,out17,out18,out19,out20,out21,out22] = F_bcs_PR_segs_rhs_0(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12,in13,in14,in15,in16,in17,in18,in19,in20,in21,in22,in23,in24,in25,in26,in27,in28,in29,in30,in31,in32,in33,in34,in35,in36,in37,in38,in39,in40,in41,in42,in43,in44,in45,in46,in47,in48,in49,in50,in51,in52,in53,in54,in55,in56,in57,in58,in59,in60,in61,in62,in63,in64,in65,in66,in67,in68,in69,in70,in71,in72,in73,in74,in75,in76,in77,in78,in79,in80,in81,in82,in83,in84,in85,in86,in87,in88,in89,in90,in91,in92,in93,in94,in95,in96,in97,in98)
%F_bcs_PR_segs_rhs_0
%    [OUT1,OUT2,OUT3,OUT4,OUT5,OUT6,OUT7,OUT8,OUT9,OUT10,OUT11,OUT12,OUT13,OUT14,OUT15,OUT16,OUT17,OUT18,OUT19,OUT20,OUT21,OUT22] = F_bcs_PR_segs_rhs_0(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12,IN13,IN14,IN15,IN16,IN17,IN18,IN19,IN20,IN21,IN22,IN23,IN24,IN25,IN26,IN27,IN28,IN29,IN30,IN31,IN32,IN33,IN34,IN35,IN36,IN37,IN38,IN39,IN40,IN41,IN42,IN43,IN44,IN45,IN46,IN47,IN48,IN49,IN50,IN51,IN52,IN53,IN54,IN55,IN56,IN57,IN58,IN59,IN60,IN61,IN62,IN63,IN64,IN65,IN66,IN67,IN68,IN69,IN70,IN71,IN72,IN73,IN74,IN75,IN76,IN77,IN78,IN79,IN80,IN81,IN82,IN83,IN84,IN85,IN86,IN87,IN88,IN89,IN90,IN91,IN92,IN93,IN94,IN95,IN96,IN97,IN98)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    05-Mar-2025 22:17:18

out1 = in1-in25;
if nargout > 1
    out2 = in2-in26;
end
if nargout > 2
    out3 = in3-in27;
end
if nargout > 3
    out4 = in7-in19;
end
if nargout > 4
    out5 = in8-in20;
end
if nargout > 5
    out6 = in9-in21;
end
if nargout > 6
    out7 = -in37.*(in1-in38+in1.*in3);
end
if nargout > 7
    t2 = -in10;
    out8 = in22+t2;
end
if nargout > 8
    t3 = -in11;
    out9 = in23+t3;
end
if nargout > 9
    t4 = -in12;
    out10 = in24+t4;
end
if nargout > 10
    out11 = -in28+in4.*in45;
end
if nargout > 11
    out12 = -in29+in5.*in45;
end
if nargout > 12
    out13 = -in30+in6.*in45;
end
if nargout > 13
    out14 = sqrt(abs(in28).^2+abs(in29).^2+abs(in30).^2)-1.0;
end
if nargout > 14
    out15 = -in1+in31;
end
if nargout > 15
    out16 = -in2+in32;
end
if nargout > 16
    out17 = -in3+in33;
end
if nargout > 17
    out18 = -in13+in16-in47.*cos(in48);
end
if nargout > 18
    out19 = -in14+in17;
end
if nargout > 19
    out20 = -in15+in18-in47.*sin(in48);
end
if nargout > 20
    out21 = t2.*(conj(in7)-conj(in34))+t3.*(conj(in8)-conj(in35))+t4.*(conj(in9)-conj(in36));
end
if nargout > 21
    out22 = -in46+(in7-in34).^2+(in8-in35).^2+(in9-in36).^2;
end
end


function [out1,out2,out3,out4,out5,out6,out7,out8,out9,out10,out11,out12,out13,out14,out15,out16,out17,out18,out19,out20,out21,out22] = F_bcs_PR_segs_rhs_1(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12,in13,in14,in15,in16,in17,in18,in19,in20,in21,in22,in23,in24,in25,in26,in27,in28,in29,in30,in31,in32,in33,in34,in35,in36,in37,in38,in39,in40,in41,in42,in43,in44,in45,in46,in47,in48,in49,in50,in51,in52,in53,in54,in55,in56,in57,in58,in59,in60,in61,in62,in63,in64,in65,in66,in67,in68,in69,in70,in71,in72,in73,in74,in75,in76,in77,in78,in79,in80,in81,in82,in83,in84,in85,in86,in87,in88,in89,in90,in91,in92,in93,in94,in95,in96,in97,in98)
%F_bcs_PR_segs_rhs_1
%    [OUT1,OUT2,OUT3,OUT4,OUT5,OUT6,OUT7,OUT8,OUT9,OUT10,OUT11,OUT12,OUT13,OUT14,OUT15,OUT16,OUT17,OUT18,OUT19,OUT20,OUT21,OUT22] = F_bcs_PR_segs_rhs_1(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12,IN13,IN14,IN15,IN16,IN17,IN18,IN19,IN20,IN21,IN22,IN23,IN24,IN25,IN26,IN27,IN28,IN29,IN30,IN31,IN32,IN33,IN34,IN35,IN36,IN37,IN38,IN39,IN40,IN41,IN42,IN43,IN44,IN45,IN46,IN47,IN48,IN49,IN50,IN51,IN52,IN53,IN54,IN55,IN56,IN57,IN58,IN59,IN60,IN61,IN62,IN63,IN64,IN65,IN66,IN67,IN68,IN69,IN70,IN71,IN72,IN73,IN74,IN75,IN76,IN77,IN78,IN79,IN80,IN81,IN82,IN83,IN84,IN85,IN86,IN87,IN88,IN89,IN90,IN91,IN92,IN93,IN94,IN95,IN96,IN97,IN98)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    05-Mar-2025 22:17:19

out1 = in50-in74;
if nargout > 1
    out2 = in51-in75;
end
if nargout > 2
    out3 = in52-in76;
end
if nargout > 3
    out4 = in56-in68;
end
if nargout > 4
    out5 = in57-in69;
end
if nargout > 5
    out6 = in58-in70;
end
if nargout > 6
    out7 = -in86.*(in1-in38+in1.*in3)-in37.*(in50-in87+in1.*in52+in3.*in50);
end
if nargout > 7
    out8 = -in59+in71;
end
if nargout > 8
    out9 = -in60+in72;
end
if nargout > 9
    out10 = -in61+in73;
end
if nargout > 10
    out11 = -in77+in4.*in94+in45.*in53;
end
if nargout > 11
    out12 = -in78+in5.*in94+in45.*in54;
end
if nargout > 12
    out13 = -in79+in6.*in94+in45.*in55;
end
if nargout > 13
    t2 = abs(in28);
    t3 = abs(in29);
    t4 = abs(in30);
    t5 = conj(in28);
    t6 = conj(in29);
    t7 = conj(in30);
    out14 = ((t2.*(in77.*t5+in28.*conj(in77)).*1.0./sqrt(in28.*t5)+t3.*(in78.*t6+in29.*conj(in78)).*1.0./sqrt(in29.*t6)+t4.*(in79.*t7+in30.*conj(in79)).*1.0./sqrt(in30.*t7)).*1.0./sqrt(t2.^2+t3.^2+t4.^2))./2.0;
end
if nargout > 14
    out15 = -in50+in80;
end
if nargout > 15
    out16 = -in51+in81;
end
if nargout > 16
    out17 = -in52+in82;
end
if nargout > 17
    t8 = cos(in48);
    t9 = sin(in48);
    out18 = -in62+in65-in96.*t8+in47.*in97.*t9;
end
if nargout > 18
    out19 = -in63+in66;
end
if nargout > 19
    out20 = -in64+in67-in96.*t9-in47.*in97.*t8;
end
if nargout > 20
    out21 = -in59.*(conj(in7)-conj(in34))-in60.*(conj(in8)-conj(in35))-in61.*(conj(in9)-conj(in36))-in10.*(conj(in56)-conj(in83))-in11.*(conj(in57)-conj(in84))-in12.*(conj(in58)-conj(in85));
end
if nargout > 21
    out22 = -in95+(in7-in34).*(in56-in83).*2.0+(in8-in35).*(in57-in84).*2.0+(in9-in36).*(in58-in85).*2.0;
end
end


function [out1,out2,out3,out4,out5,out6,out7,out8,out9,out10,out11,out12,out13,out14,out15,out16,out17,out18,out19,out20,out21,out22] = F_bcs_PR_segs_rhs_2(in1,in2,in3,in4,in5,in6,in7,in8,in9,in10,in11,in12,in13,in14,in15,in16,in17,in18,in19,in20,in21,in22,in23,in24,in25,in26,in27,in28,in29,in30,in31,in32,in33,in34,in35,in36,in37,in38,in39,in40,in41,in42,in43,in44,in45,in46,in47,in48,in49,in50,in51,in52,in53,in54,in55,in56,in57,in58,in59,in60,in61,in62,in63,in64,in65,in66,in67,in68,in69,in70,in71,in72,in73,in74,in75,in76,in77,in78,in79,in80,in81,in82,in83,in84,in85,in86,in87,in88,in89,in90,in91,in92,in93,in94,in95,in96,in97,in98)
%F_bcs_PR_segs_rhs_2
%    [OUT1,OUT2,OUT3,OUT4,OUT5,OUT6,OUT7,OUT8,OUT9,OUT10,OUT11,OUT12,OUT13,OUT14,OUT15,OUT16,OUT17,OUT18,OUT19,OUT20,OUT21,OUT22] = F_bcs_PR_segs_rhs_2(IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11,IN12,IN13,IN14,IN15,IN16,IN17,IN18,IN19,IN20,IN21,IN22,IN23,IN24,IN25,IN26,IN27,IN28,IN29,IN30,IN31,IN32,IN33,IN34,IN35,IN36,IN37,IN38,IN39,IN40,IN41,IN42,IN43,IN44,IN45,IN46,IN47,IN48,IN49,IN50,IN51,IN52,IN53,IN54,IN55,IN56,IN57,IN58,IN59,IN60,IN61,IN62,IN63,IN64,IN65,IN66,IN67,IN68,IN69,IN70,IN71,IN72,IN73,IN74,IN75,IN76,IN77,IN78,IN79,IN80,IN81,IN82,IN83,IN84,IN85,IN86,IN87,IN88,IN89,IN90,IN91,IN92,IN93,IN94,IN95,IN96,IN97,IN98)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    05-Mar-2025 22:17:19

out1 = 0.0;
if nargout > 1
    out2 = 0.0;
end
if nargout > 2
    out3 = 0.0;
end
if nargout > 3
    out4 = 0.0;
end
if nargout > 4
    out5 = 0.0;
end
if nargout > 5
    out6 = 0.0;
end
if nargout > 6
    out7 = in86.*(in50-in87+in1.*in52+in3.*in50).*-2.0-in37.*in50.*in52.*2.0;
end
if nargout > 7
    out8 = 0.0;
end
if nargout > 8
    out9 = 0.0;
end
if nargout > 9
    out10 = 0.0;
end
if nargout > 10
    out11 = in53.*in94.*2.0;
end
if nargout > 11
    out12 = in54.*in94.*2.0;
end
if nargout > 12
    out13 = in55.*in94.*2.0;
end
if nargout > 13
    t2 = abs(in28);
    t3 = abs(in29);
    t4 = abs(in30);
    t5 = conj(in28);
    t6 = conj(in29);
    t7 = conj(in30);
    t8 = conj(in77);
    t9 = conj(in78);
    t10 = conj(in79);
    t11 = cos(in48);
    t12 = sin(in48);
    t13 = in97.^2;
    t14 = t2.^2;
    t15 = t3.^2;
    t16 = t4.^2;
    t17 = in28.*t5;
    t18 = in29.*t6;
    t19 = in30.*t7;
    t20 = in28.*t8;
    t21 = in77.*t5;
    t22 = in29.*t9;
    t23 = in78.*t6;
    t24 = in30.*t10;
    t25 = in79.*t7;
    t26 = t20+t21;
    t27 = t22+t23;
    t28 = t24+t25;
    t29 = 1.0./sqrt(t17);
    t30 = 1.0./sqrt(t18);
    t31 = 1.0./sqrt(t19);
    t35 = t14+t15+t16;
    t32 = t26.^2;
    t33 = t27.^2;
    t34 = t28.^2;
    out14 = 1.0./t35.^(3.0./2.0).*(t2.*t26.*t29+t3.*t27.*t30+t4.*t28.*t31).^2.*(-1.0./4.0)+(1.0./sqrt(t35).*((t29.^2.*t32)./2.0+(t30.^2.*t33)./2.0+(t31.^2.*t34)./2.0-(t2.*t29.^3.*t32)./2.0-(t3.*t30.^3.*t33)./2.0-(t4.*t31.^3.*t34)./2.0+in77.*t2.*t8.*t29.*2.0+in78.*t3.*t9.*t30.*2.0+in79.*t4.*t10.*t31.*2.0))./2.0;
end
if nargout > 14
    out15 = 0.0;
end
if nargout > 15
    out16 = 0.0;
end
if nargout > 16
    out17 = 0.0;
end
if nargout > 17
    out18 = in96.*in97.*t12.*2.0+in47.*t11.*t13;
end
if nargout > 18
    out19 = 0.0;
end
if nargout > 19
    out20 = in96.*in97.*t11.*-2.0+in47.*t12.*t13;
end
if nargout > 20
    out21 = in59.*(conj(in56)-conj(in83)).*-2.0-in60.*(conj(in57)-conj(in84)).*2.0-in61.*(conj(in58)-conj(in85)).*2.0;
end
if nargout > 21
    out22 = (in56-in83).^2.*2.0+(in57-in84).^2.*2.0+(in58-in85).^2.*2.0;
end
end

