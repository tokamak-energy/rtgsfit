# GSFit settings for `st40_magrun16`

GSFit's `st40_setup_for_rtgsfit` settings, made as close as possible to those of the post-shot
GSFit reconstructions (GSFIT RUN04, e.g. pulse 14816). RT-GSFit's `st40_magrun16` preshot,
`RTGSFIT` 99000230 `PRESHOT09`, is built from them, and `tests/test_st40_magrun16.py` runs GSFit
itself with them.

* `GSFIT_code_settings.json`: the `st40_mdsplus_with_digital_filter` reader with ELMAG RUN16, MAG
  RUN16 and PSU2COIL BEST (the run RUN04 used for each pulse: RUN06 for 14816, RUN02 for 11089),
  and RUN04's 81x161 grid. RT-GSFit's preshot uses its own `gsfit/grid`.
* `sensor_weights_flux_loops.json`, `sensor_weights_bp_probe.json`, `sensor_weights_rogowski_coils.json`:
  RUN04's archived `include`, `weight` and `expected_value`.
* `source_function_p_prime.json`, `source_function_ff_prime.json`: **`n_dof: 1`** rather than RUN04's 2,
  so that GSFit's p' and FF' have RT-GSFit's single degree of freedom each. This is the one intended
  difference from RUN04; with `n_dof: 2` GSFit reproduces RUN04.

`st40_mdsplus_with_digital_filter` does not divide the flux-loop weights by 2 pi, unlike GSFit's
other readers, so its flux-loop weights are effectively 2 pi larger. That is what RUN04 did.
