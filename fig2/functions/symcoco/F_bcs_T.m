function varargout=F_bcs_T(action,varargin)
%% Automatically generated with matlabFunction
%#ok<*DEFNU,*INUSD,*INUSL>

switch action
  case 'nargs'
   varargout{1}=1;
   return
  case 'nout'
   varargout{1}=1;
   return
  case 'argrange'
   varargout{1}=struct('u',1:1);
   return
  case 'argsize'
   varargout{1}=struct('u',1);
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
nout=1;
order=varargin{1};
f=str2func(sprintf('F_bcs_T_%s_%d',action,order));
varargout=cell(nout,1);
[varargout{:}]=f(varargin{2:end});
end



function out1 = F_bcs_T_rhs_0(in1,in2)
%F_bcs_T_rhs_0
%    OUT1 = F_bcs_T_rhs_0(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    11-Mar-2025 14:12:41

out1 = in1-1.0;
end


function out1 = F_bcs_T_rhs_1(in1,in2)
%F_bcs_T_rhs_1
%    OUT1 = F_bcs_T_rhs_1(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    11-Mar-2025 14:12:41

out1 = in2;
end


function out1 = F_bcs_T_rhs_2(in1,in2)
%F_bcs_T_rhs_2
%    OUT1 = F_bcs_T_rhs_2(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    11-Mar-2025 14:12:41

out1 = 0.0;
end

