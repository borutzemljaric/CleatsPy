import cProfile, pstats, io

pr = cProfile.Profile()
pr.enable()

# These should satisfy ydot = f(t, y) at t=0
zacpog_ode_dot = sistem_enacb(
    0,  # t=0
    zacpog_ode,  # Current state [positions; velocities]
    U, V, hh, ne_c, tipverige,
    e_zp, Cdi, Qe_gg,
    UV_indices,
    qc, qv, Ce, ni, nd, M,
    odvisneU, odvisneUdot,
    cab, le, le_full, Ip
)

res = make_system_dae(
    U, V, hh, ne_c, tipverige, e_zp, Cdi, Qe_gg, UV_indices,
    qc, qv, Ce, ni, nd, M, odvisneU, odvisneUdot,
    cab, le, le_full, Ip
)

# Returns: [velocities; accelerations]  which is exactly ydot!


pr.disable()

s = io.StringIO()
ps = pstats.Stats(pr, stream=s).sort_stats("cumulative")
ps.print_stats(20)
print(s.getvalue())
breakpoint()