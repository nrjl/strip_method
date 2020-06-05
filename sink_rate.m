function edot = sink_rate(w_in, V, h)

% global airspeed
global w_wind
w_wind = w_in;

[X, U] = trim_controls(V, h);
Xdot = state_rates(X, U);

edot = Xdot(12);
fprintf(1, 'Sink rate = %0.5g m/s\n', edot);