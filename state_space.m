function [A, B] = state_space(f_xdot, X0, U0)
%--------------------------------------------------------------------------
%
% FUNCTION:		state_space
%
% PURPOSE:		Calculate the state space representation of a system near a
%					given trim condition
%
% SYNTAX:		[A, B] = state_space(f_xdot, X0, U0)
%
% INPUTS:		f_xdot	- Handle to derivate function xdot = f_xdot(X, U)
%				X0		- Trimmed condition state vector
%				U0		- Trimmed condition control vector
%
% OUTPUTS:		A	- System to sytem propagation matrix
%				B	- Control to system propagation matrix
%				Xdot = A(X-X0) + B(X-X0)
%
% AUTHOR:		Nicholas Lawrance
%
% CREATED:		November 2007
%
% MODIFIED:     November 2007
%
% See also:		control_test, plot_states
%--------------------------------------------------------------------------

% --- INTIALISATION --- %

% Calculate initial state rate vector for equilibrium condition (should be
% ~0)

Xdot0 = f_xdot(X0, U0);		% Assign initial state rate vector

delta_x = 1e-6;				% State perturbation value
delta_u = 1e-6;				% Control perturbation value

A = zeros(length(X0));		% Initialise A matrix
B = zeros(length(Xdot0), length(U0));
							% Initialise B matrix

% --- CALCULATE A MATRIX --- %
U = U0;

for i = 1:length(Xdot0)		% For all state variables
    X = X0;					% Reset states to equilibrium values
    X(i) = X(i) + delta_x;	% Add perturbation to a single state variable

	% Evaluate state rate change due to perturbation (from non-linear sim)
    Xdot = f_xdot(X, U);

	% Assign column of A for each state variable
    A(:,i) = (Xdot - Xdot0)/(delta_x);
end


% --- CALCULATE B MATRIX --- %
X = X0;						% Reset state variables vector

for i = 1:length(U0)		% For all control variables
    U = U0;					% Reset control variable vector
    
    U(i) = U(i) + delta_u;	% Add perturbation to single control variable
    
	Xdot = f_xdot(X, U);

	% Assign column of B for each control variable
    B(:,i) = (Xdot - Xdot0)/(delta_u);
end