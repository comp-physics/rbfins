function d = dairfoil(p, x_center, y_center, naca_digits, chord, angle_deg)
  %DAIRFOIL Signed distance function for NACA 4-digit airfoil
  %
  % This function computes the signed distance from points to a NACA 4-digit
  % airfoil boundary. The airfoil is defined by the NACA 4-digit series
  % parameterization with specified chord length, center position, and angle of attack.
  %
  % INPUTS:
  %   p           - Points to evaluate [N x 2] array of (x,y) coordinates
  %   x_center    - X-coordinate of airfoil center (leading edge position)
  %   y_center    - Y-coordinate of airfoil center (leading edge position)
  %   naca_digits - NACA 4-digit parameters [m, p, t1, t2] where:
  %                 m = maximum camber as percentage of chord (first digit)
  %                 p = position of maximum camber as percentage of chord (second digit)
  %                 t1, t2 = thickness as percentage of chord (last two digits)
  %   chord       - Chord length of the airfoil
  %   angle_deg   - Angle of attack in degrees
  %
  % OUTPUTS:
  %   d           - Signed distance to airfoil boundary [N x 1]
  %                 Negative inside, positive outside
  %
  % ALGORITHM:
  %   Uses NACA 4-digit series equations for mean camber line and thickness distribution.
  %   Distance is computed by finding the minimum distance to the airfoil surface.
  %
  % REFERENCE:
  %   Abbott, I. H., & Von Doenhoff, A. E. (1959). Theory of Wing Sections.
  %   Dover Publications.

  % Extract NACA parameters
  m = naca_digits(1) / 100;        % Maximum camber (as fraction of chord)
  p_camber = naca_digits(2) / 10;  % Position of maximum camber (as fraction of chord)
  t = (naca_digits(3) * 10 + naca_digits(4)) / 100; % Maximum thickness (as fraction of chord)

  % Convert angle to radians
  angle_rad = angle_deg * pi / 180;

  % Translate points to airfoil coordinate system (leading edge at origin)
  p_translated = p - [x_center, y_center];

  % Rotate points by negative angle of attack to align with airfoil coordinates
  cos_a = cos(-angle_rad);
  sin_a = sin(-angle_rad);
  rotation_matrix = [cos_a, -sin_a; sin_a, cos_a];
  p_rotated = (rotation_matrix * p_translated')';

  % Normalize by chord length
  x_norm = p_rotated(:, 1) / chord;
  y_norm = p_rotated(:, 2) / chord;

  % Initialize distance array
  d = zeros(size(p, 1), 1);

  % Process each point
  for i = 1:size(p, 1)
    x = x_norm(i);
    y = y_norm(i);

    % Check if point is within airfoil x-bounds
    if x < 0 || x > 1
      % Point is outside chord bounds - compute distance to nearest endpoint
      if x < 0
        % Distance to leading edge (0, 0)
        d(i) = sqrt(x^2 + y^2);
      else
        % Distance to trailing edge (1, 0)
        d(i) = sqrt((x - 1)^2 + y^2);
      end
    else
      % Point is within chord bounds - compute distance to airfoil surface
      d(i) = compute_airfoil_distance(x, y, m, p_camber, t);
    end
  end

  % Scale back to physical coordinates
  d = d * chord;
end

function dist = compute_airfoil_distance(x, y, m, p_camber, t)
  %COMPUTE_AIRFOIL_DISTANCE Compute distance from point to NACA airfoil surface
  %
  % Uses the NACA 4-digit series equations to compute the mean camber line
  % and thickness distribution, then finds the minimum distance to the surface.

  % Compute mean camber line and its derivative
  if x <= p_camber && p_camber > 0
    % Forward portion of camber line
    yc = (m / p_camber^2) * (2 * p_camber * x - x^2);
    dyc_dx = (m / p_camber^2) * (2 * p_camber - 2 * x);
  elseif x > p_camber && p_camber > 0
    % Aft portion of camber line
    yc = (m / (1 - p_camber)^2) * ((1 - 2 * p_camber) + 2 * p_camber * x - x^2);
    dyc_dx = (m / (1 - p_camber)^2) * (2 * p_camber - 2 * x);
  else
    % Symmetric airfoil (m = 0 or p_camber = 0)
    yc = 0;
    dyc_dx = 0;
  end

  % Compute thickness distribution using NACA equation
  yt = (t / 0.2) * (0.2969 * sqrt(x) - 0.1260 * x - 0.3516 * x^2 + 0.2843 * x^3 - 0.1015 * x^4);

  % Compute angle of camber line
  theta = atan(dyc_dx);

  % Upper and lower surface coordinates
  xu = x - yt * sin(theta);
  yu = yc + yt * cos(theta);
  xl = x + yt * sin(theta);
  yl = yc - yt * cos(theta);

  % Compute distances to upper and lower surfaces
  dist_upper = sqrt((x - xu)^2 + (y - yu)^2);
  dist_lower = sqrt((x - xl)^2 + (y - yl)^2);

  % Minimum distance to surface
  dist = min(dist_upper, dist_lower);

  % Determine sign (inside/outside) using simple y-coordinate comparison
  % Point is inside if it's between upper and lower surfaces at this x-location
  if y >= yl && y <= yu
    % Point is between upper and lower surfaces - inside airfoil
    % If point is outside (y > yu or y < yl), dist remains positive
    dist = -dist;
  end
end
