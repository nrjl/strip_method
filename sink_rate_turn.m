function edot = sink_rate_turn(w_in, V, bank, h)

global w_wind

w_wind = w_in;

[X, U] = trim_bankedturn(V, bank, h);
Xdot = state_rates(X, U);

edot = Xdot(12);
fprintf(1, 'Sink rate = %0.5g m/s\n', edot);