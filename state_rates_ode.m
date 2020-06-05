function X_dot = state_rates_ode(t, y, U, aero, param)

X_dot = state_rates(y, U, aero, param);