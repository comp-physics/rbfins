function stability_diagnostics()
  %STABILITY_DIAGNOSTICS Comprehensive stability analysis and debugging functions
  %
  % This file contains all debugging and stability analysis functions for the
  % RBF-FD Navier-Stokes simulation. Functions are designed to provide detailed
  % root cause analysis when simulations become unstable.
  %
  % Main functions:
  %   - log_step_stats: Print diagnostic statistics for current time step
  %   - analyze_stability: Comprehensive stability analysis and warnings
  %   - analyze_instability_cause: Detailed analysis of what caused instability
  %   - provide_stability_recommendations: Give actionable advice for stability issues
  %   - classify_node_idx: Determine location type of a node index
  %   - save_debug_snapshot: Save simulation state for debugging
  %
  % Usage: These functions are called automatically from simulate.m when
  %        debug logging is enabled via config.logging parameters.
  
  error('This is a function library. Call individual functions directly.');
end

function log_step_stats(j, Wprev, Wcur, dt, Dx, Dy, D0_12_x, D0_12_y, xy1, h_min, h_med, cfg) %#ok<DEFNU>
  %LOG_STEP_STATS Print diagnostic statistics for current time step

  Uprev = Wprev(1:end / 2);
  Vprev = Wprev(end / 2 + 1:end);
  U = Wcur(1:end / 2);
  V = Wcur(end / 2 + 1:end);

  % Basic statistics
  ke = 0.5 * (sum(U.^2) + sum(V.^2));
  dWdt = norm(Wcur - Wprev, 2) / max(dt, eps);
  div = D0_12_x * U + D0_12_y * V;
  div_l2 = norm(div, 2);
  div_max = max(abs(div));
  umax = max(abs(U));
  vmax = max(abs(V));
  umax_i = find(abs(U) == umax, 1);
  vmax_i = find(abs(V) == vmax, 1);
  pt_u = xy1(umax_i, :);
  pt_v = xy1(vmax_i, :);

  % CFL analysis
  Vmag_max = max(sqrt(max(U.^2 + V.^2, 0)));
  CFL_min = NaN;
  CFL_med = NaN;
  if isfinite(h_min) && h_min > 0
    CFL_min = Vmag_max * dt / h_min;
  end
  if isfinite(h_med) && h_med > 0
    CFL_med = Vmag_max * dt / h_med;
  end

  % Basic output
  fprintf('[step %5d] KE=%.3e |dW/dt|=%.3e div(L2)=%.3e div(max)=%.3e |u|max=%.3e@[(%.3g,%.3g)] |v|max=%.3e@[(%.3g,%.3g)] CFL[min/med]=[%.3g/%.3g]\n', ...
          j, ke, dWdt, div_l2, div_max, umax, pt_u(1), pt_u(2), vmax, pt_v(1), pt_v(2), CFL_min, CFL_med);

  % Enhanced stability analysis
  if cfg.logging.stability_analysis
    analyze_stability(j, U, V, Uprev, Vprev, div, div_l2, div_max, CFL_min, CFL_med, dt, xy1, cfg);
  end
end

function analyze_stability(j, U, V, Uprev, Vprev, div, div_l2, div_max, CFL_min, CFL_med, dt, xy1, cfg)
  %ANALYZE_STABILITY Comprehensive stability analysis and warnings

  warnings = {};
  critical_issues = {};

  % 1. CFL Condition Analysis
  if cfg.logging.cfl_monitoring
    if isfinite(CFL_min) && CFL_min > cfg.logging.cfl_critical_threshold
      critical_issues{end + 1} = sprintf('CRITICAL CFL violation: %.3f > %.3f (time step too large)', CFL_min, cfg.logging.cfl_critical_threshold);
    elseif isfinite(CFL_min) && CFL_min > cfg.logging.cfl_warning_threshold
      warnings{end + 1} = sprintf('CFL warning: %.3f > %.3f (approaching instability)', CFL_min, cfg.logging.cfl_warning_threshold);
    end
  end

  % 2. Mass Conservation Analysis
  if cfg.logging.mass_conservation
    if div_max > cfg.logging.divergence_critical_threshold
      critical_issues{end + 1} = sprintf('CRITICAL mass conservation failure: max|div|=%.3e > %.3e', div_max, cfg.logging.divergence_critical_threshold);
    elseif div_l2 > cfg.logging.divergence_warning_threshold
      warnings{end + 1} = sprintf('Mass conservation degrading: L2(div)=%.3e > %.3e', div_l2, cfg.logging.divergence_warning_threshold);
    end
  end

  % 3. Velocity Explosion Detection
  Vmag_max = max(sqrt(U.^2 + V.^2));
  if Vmag_max > cfg.logging.velocity_explosion_threshold
    critical_issues{end + 1} = sprintf('CRITICAL velocity explosion: max|V|=%.3e > %.3e', Vmag_max, cfg.logging.velocity_explosion_threshold);
  end

  % 4. Velocity Gradient Analysis
  if cfg.logging.velocity_gradients
    try
      % Compute velocity gradients (simplified check)
      dU = U - Uprev;
      dV = V - Vprev;
      vel_change_rate = max(sqrt(dU.^2 + dV.^2)) / dt;
      if vel_change_rate > cfg.logging.velocity_explosion_threshold / dt
        warnings{end + 1} = sprintf('High velocity change rate: %.3e (rapid acceleration detected)', vel_change_rate);
      end
    catch
      % Skip if gradient computation fails
    end
  end

  % 5. Energy Conservation Analysis
  if cfg.logging.energy_conservation
    ke_current = 0.5 * sum(U.^2 + V.^2);
    ke_prev = 0.5 * sum(Uprev.^2 + Vprev.^2);
    energy_change_rate = abs(ke_current - ke_prev) / (dt * max(ke_prev, eps));
    if energy_change_rate > 10  % Arbitrary threshold for rapid energy change
      warnings{end + 1} = sprintf('Rapid energy change: dKE/dt=%.3e (potential numerical instability)', energy_change_rate);
    end
  end

  % 6. Boundary Analysis
  if cfg.logging.boundary_analysis
    % Check for extreme values near boundaries (simplified)
    n_interior = length(U);
    if n_interior < length(xy1)  % Has boundary nodes
      U_boundary = U((n_interior + 1):end);
      V_boundary = V((n_interior + 1):end);
      if any(abs(U_boundary) > cfg.logging.velocity_explosion_threshold / 10) || any(abs(V_boundary) > cfg.logging.velocity_explosion_threshold / 10)
        warnings{end + 1} = sprintf('High boundary velocities detected (BC enforcement issues?)');
      end
    end
  end

  % 7. Mesh Quality Indicators
  if cfg.logging.mesh_quality_check && j <= 5  % Only check early in simulation
    % Simple mesh quality check based on velocity field smoothness
    if length(U) > 100  % Only for reasonable mesh sizes
      U_std = std(U);
      V_std = std(V);
      U_mean = mean(abs(U));
      V_mean = mean(abs(V));
      if U_std > 5 * U_mean || V_std > 5 * V_mean
        warnings{end + 1} = sprintf('High velocity field variation (mesh quality issues?)');
      end
    end
  end

  % Output warnings and critical issues
  if ~isempty(critical_issues)
    fprintf('  *** CRITICAL STABILITY ISSUES ***\n');
    for i = 1:length(critical_issues)
      fprintf('  [CRITICAL] %s\n', critical_issues{i});
    end
  end

  if ~isempty(warnings)
    fprintf('  *** STABILITY WARNINGS ***\n');
    for i = 1:length(warnings)
      fprintf('  [WARNING] %s\n', warnings{i});
    end
  end

  % Provide recommendations
  if ~isempty(critical_issues) || ~isempty(warnings)
    provide_stability_recommendations(critical_issues, warnings, CFL_min, div_max, Vmag_max, cfg);
  end
end

function provide_stability_recommendations(critical_issues, warnings, CFL_min, div_max, Vmag_max, cfg)
  %PROVIDE_STABILITY_RECOMMENDATIONS Give actionable advice for stability issues

  fprintf('  *** STABILITY RECOMMENDATIONS ***\n');

  % CFL-based recommendations
  if isfinite(CFL_min) && CFL_min > cfg.logging.cfl_warning_threshold
    new_dt = cfg.simulation.time_step * cfg.logging.cfl_warning_threshold / CFL_min;
    fprintf('  [RECOMMEND] Reduce time step: dt = %.3e (current: %.3e)\n', new_dt, cfg.simulation.time_step);
  end

  % Mass conservation recommendations
  if div_max > cfg.logging.divergence_warning_threshold
    fprintf('  [RECOMMEND] Check pressure solver: increase iterations or reduce tolerance\n');
    fprintf('  [RECOMMEND] Verify boundary conditions: ensure proper no-slip/outlet BCs\n');
    fprintf('  [RECOMMEND] Check mesh quality: refine near obstacles (reduce refine_a1, refine_b1)\n');
  end

  % Velocity explosion recommendations
  if Vmag_max > cfg.logging.velocity_explosion_threshold
    fprintf('  [RECOMMEND] Reduce time step significantly: dt < %.3e\n', cfg.simulation.time_step * 0.1);
    fprintf('  [RECOMMEND] Check initial conditions: ensure smooth velocity initialization\n');
    fprintf('  [RECOMMEND] Verify geometry: check for sharp corners or mesh degeneracies\n');
  end

  % General recommendations
  fprintf('  [RECOMMEND] Monitor: Enable cfg.logging.trace_substeps to identify which substep fails\n');
  fprintf('  [RECOMMEND] Mesh: Try finer mesh (reduce cfg.mesh.dist) or better refinement\n');
  fprintf('  [RECOMMEND] Algorithm: Consider lower Reynolds number for testing\n');
end

function analyze_instability_cause(j, W, p0, bad_idx, loc_str, pt, dt, Dx, Dy, D0_12_x, D0_12_y, xy1, h_min, h_med, cfg) %#ok<DEFNU>
  %ANALYZE_INSTABILITY_CAUSE Detailed analysis of what caused the instability

  fprintf('\n=== INSTABILITY ROOT CAUSE ANALYSIS ===\n');

  % Extract current and previous solutions
  U_cur = W(1:end / 2, j + 1);
  V_cur = W(end / 2 + 1:end, j + 1);
  if j > 1
    U_prev = W(1:end / 2, j);
    V_prev = W(end / 2 + 1:end, j);
  else
    U_prev = zeros(size(U_cur));
    V_prev = zeros(size(V_cur));
  end

  % 1. Location-specific analysis
  fprintf('1. LOCATION ANALYSIS:\n');
  fprintf('   Failed at: %s node at (%.3g, %.3g)\n', loc_str, pt(1), pt(2));

  if strcmp(loc_str, 'obstacle')
    fprintf('   -> OBSTACLE BOUNDARY: Likely boundary condition enforcement issue\n');
    fprintf('   -> Check: No-slip BC implementation, mesh quality near obstacle\n');
  elseif strcmp(loc_str, 'inlet')
    fprintf('   -> INLET BOUNDARY: Inflow condition issue\n');
    fprintf('   -> Check: Inlet velocity profile, boundary condition implementation\n');
  elseif strcmp(loc_str, 'outlet')
    fprintf('   -> OUTLET BOUNDARY: Outflow condition issue\n');
    fprintf('   -> Check: Pressure BC, convective outlet conditions\n');
  elseif strcmp(loc_str, 'wall')
    fprintf('   -> WALL BOUNDARY: Wall boundary condition issue\n');
    fprintf('   -> Check: No-slip enforcement, mesh quality at walls\n');
  else
    fprintf('   -> INTERIOR: Flow field instability\n');
    fprintf('   -> Check: CFL condition, mesh quality, pressure-velocity coupling\n');
  end

  % 2. Temporal analysis
  fprintf('\n2. TEMPORAL ANALYSIS:\n');
  if j <= 3
    fprintf('   Failed in startup phase (step %d)\n', j);
    fprintf('   -> Likely: Initial condition issues, startup scheme problems\n');
    fprintf('   -> Check: Initial velocity field, startup time step\n');
  elseif j <= 20
    fprintf('   Failed in early simulation (step %d)\n', j);
    fprintf('   -> Likely: Transient adjustment issues, mesh/BC problems\n');
    fprintf('   -> Check: Mesh quality, boundary condition implementation\n');
  else
    fprintf('   Failed in developed flow (step %d)\n', j);
    fprintf('   -> Likely: Accumulated numerical errors, CFL violation\n');
    fprintf('   -> Check: Time step size, numerical scheme stability\n');
  end

  % 3. Field magnitude analysis
  fprintf('\n3. FIELD MAGNITUDE ANALYSIS:\n');
  U_finite = U_cur(isfinite(U_cur));
  V_finite = V_cur(isfinite(V_cur));

  if ~isempty(U_finite) && ~isempty(V_finite)
    U_max_finite = max(abs(U_finite));
    V_max_finite = max(abs(V_finite));
    fprintf('   Max finite |U|: %.3e, Max finite |V|: %.3e\n', U_max_finite, V_max_finite);

    if U_max_finite > 100 || V_max_finite > 100
      fprintf('   -> VELOCITY EXPLOSION detected before NaN\n');
      fprintf('   -> Cause: Time step too large, mesh too coarse, or BC issues\n');
    end
  end

  % 4. CFL analysis at failure
  fprintf('\n4. CFL CONDITION AT FAILURE:\n');
  if ~isempty(U_finite) && ~isempty(V_finite) && isfinite(h_min) && h_min > 0
    Vmag_max = max(sqrt(U_finite.^2 + V_finite.^2));
    CFL_at_failure = Vmag_max * dt / h_min;
    fprintf('   CFL number at failure: %.3f\n', CFL_at_failure);

    if CFL_at_failure > 1.0
      fprintf('   -> CRITICAL CFL VIOLATION: CFL > 1.0\n');
      fprintf('   -> Solution: Reduce time step by factor %.1f\n', CFL_at_failure);
    elseif CFL_at_failure > 0.5
      fprintf('   -> HIGH CFL: Approaching stability limit\n');
      fprintf('   -> Solution: Reduce time step by factor %.1f\n', CFL_at_failure * 2);
    else
      fprintf('   -> CFL appears reasonable, other causes likely\n');
    end
  end

  % 5. Mass conservation at failure
  fprintf('\n5. MASS CONSERVATION AT FAILURE:\n');
  if ~isempty(U_finite) && ~isempty(V_finite)
    try
      div = D0_12_x * U_finite + D0_12_y * V_finite;
      div_max = max(abs(div));
      div_rms = sqrt(mean(div.^2));
      fprintf('   Max |divergence|: %.3e, RMS divergence: %.3e\n', div_max, div_rms);

      if div_max > 1e-1
        fprintf('   -> CRITICAL MASS CONSERVATION FAILURE\n');
        fprintf('   -> Cause: Pressure solver issues, poor mesh quality\n');
      elseif div_max > 1e-3
        fprintf('   -> MASS CONSERVATION DEGRADED\n');
        fprintf('   -> Cause: Numerical errors accumulating\n');
      end
    catch
      fprintf('   -> Cannot compute divergence (too many NaNs)\n');
    end
  end

  % 6. Pressure field analysis
  fprintf('\n6. PRESSURE FIELD ANALYSIS:\n');
  if ~isempty(p0) && any(isfinite(p0))
    p_finite = p0(isfinite(p0));
    p_max = max(abs(p_finite));
    fprintf('   Max |pressure|: %.3e\n', p_max);

    if p_max > 1000
      fprintf('   -> PRESSURE EXPLOSION detected\n');
      fprintf('   -> Cause: Pressure solver instability, BC issues\n');
    end
  end

  % 7. Specific recommendations
  fprintf('\n7. SPECIFIC RECOMMENDATIONS FOR THIS FAILURE:\n');

  % Time step recommendation
  if isfinite(h_min) && h_min > 0
    safe_dt = 0.1 * h_min;  % Very conservative
    fprintf('   [CRITICAL] Reduce time step to: dt = %.3e (current: %.3e)\n', safe_dt, dt);
  end

  % Mesh recommendation
  if strcmp(loc_str, 'obstacle') || strcmp(loc_str, 'wall')
    fprintf('   [CRITICAL] Refine mesh near boundaries: reduce refine_a1 to %.3f\n', cfg.mesh.refine_a1 * 0.5);
  end

  % Geometry-specific recommendations
  if strcmp(cfg.geometry.type, 'airfoil')
    fprintf('   [AIRFOIL] Try: smaller angle of attack, finer mesh, smaller time step\n');
  elseif strcmp(cfg.geometry.type, 'multi')
    fprintf('   [MULTI-OBSTACLE] Try: increase obstacle separation, reduce Reynolds number\n');
  end

  % Algorithm recommendations
  fprintf('   [ALGORITHM] Enable all debug flags: set debug_master = true in config.m\n');
  fprintf('   [ALGORITHM] Try lower Reynolds number: Re = 50 or Re = 20\n');
  fprintf('   [ALGORITHM] Check mesh generation: look for degenerate triangles\n');

  fprintf('\n=== END INSTABILITY ANALYSIS ===\n\n');
end

function [loc_str, pt] = classify_node_idx(idx, L_W, boundary_y, boundary_out, boundary_in, boundary_obs, xy1) %#ok<DEFNU>
  %CLASSIFY_NODE_IDX Determine location type of a node index
  N_y = size(boundary_y, 1);
  N_out = size(boundary_out, 1);
  N_inlet = size(boundary_in, 1);
  N_obs = size(boundary_obs, 1);

  if idx <= L_W
    loc_str = 'interior';
  else
    k = idx - L_W;
    if k <= N_y
      loc_str = 'wall';
    elseif k <= N_y + N_out
      loc_str = 'outlet';
    elseif k <= N_y + N_out + N_inlet
      loc_str = 'inlet';
    else
      loc_str = 'obstacle';
    end
  end
  pt = xy1(idx, :);
end

function save_debug_snapshot(cfg, j, xy1, W, p, dt_val, Dx, Dy, D0_12_x, D0_12_y, G) %#ok<DEFNU>
  %SAVE_DEBUG_SNAPSHOT Save simulation state for debugging
  try
    ddir = cfg.logging.snapshot_dir;
    if exist(ddir, 'dir') ~= 7
      mkdir(ddir);
    end
    fname = fullfile(ddir, sprintf('snapshot_step_%05d.mat', j));
    meta.cfg = cfg;
    meta.step = j;
    meta.time_step = dt_val;
    meta.timestamp = datestr(now);
    U_prev = W(1:end / 2, max(j, 1));
    V_prev = W(end / 2 + 1:end, max(j, 1));
    U_cur  = W(1:end / 2, j + 1);
    V_cur  = W(end / 2 + 1:end, j + 1);
    div_cur = D0_12_x * U_cur + D0_12_y * V_cur;
    save(fname, 'xy1', 'W', 'p', 'Dx', 'Dy', 'D0_12_x', 'D0_12_y', 'G', ...
         'U_prev', 'V_prev', 'U_cur', 'V_cur', 'div_cur', 'meta', '-v7.3');
    fprintf('[snapshot] saved %s\n', fname);
  catch ME
    warning('%s', ['snapshot failed: ' ME.message]);
  end
end
