function [] = LISA_animation()
    % LISA orbitography animation
    %
    % Author & support nicolas.douillet (at) free.fr, 2007-2020.
    
    % Compute LISA orbits
    [Phi, E, Sp1path, Sp2path, Sp3path, C_total, Spsup, Spinf, R] = simu_LISA_orbits();

    % Computational parameters
    step = 5;

    % Display parameters
    time_lapse = 0.1;
    title_text = 'Orbitography of space interferometer project LISA, following Earth in the solar system';
    filename = 'LISA_orbitography.gif';
    az = -75;
    el = 20;

    % Display settings
    h = figure;
    set(h, 'Position', get(0, 'ScreenSize'));
    set(gcf, 'Color', [0 0 0]);
    axis tight manual;
    orbit_bounds = false;
    LISA_relative_circular_orbit = false;

    % Plot LISA orbits
    for k = 1:step:length(Phi)
        clf;
        hold on;
        plot3(E(1, :), E(2, :), E(3, :), 'Color', [0 0 1], 'LineWidth', 3);
        plot3(Sp1path(1, :), Sp1path(2, :), Sp1path(3, :), 'Color', [1 1 0], 'LineWidth', 2);
        plot3(Sp2path(1, :), Sp2path(2, :), Sp2path(3, :), 'Color', [1 0 1], 'LineWidth', 2);
        plot3(Sp3path(1, :), Sp3path(2, :), Sp3path(3, :), 'Color', [0 1 1], 'LineWidth', 2);

        if orbit_bounds
            plot3(Spinf(1, :), Spinf(2, :), Spinf(3, :), 'Color', [0 1 0]);
            plot3(Spsup(1, :), Spsup(2, :), Spsup(3, :), 'Color', [0 1 0]);
        end

        plot3(Sp1path(1, k), Sp1path(2, k), Sp1path(3, k), 'o', 'Color', [1 1 0], 'LineWidth', 3);
        plot3(Sp2path(1, k), Sp2path(2, k), Sp2path(3, k), 'o', 'Color', [1 0 1], 'LineWidth', 3);
        plot3(Sp3path(1, k), Sp3path(2, k), Sp3path(3, k), 'o', 'Color', [0 1 1], 'LineWidth', 3);

        Spcft = [Sp1path(1, k) Sp2path(1, k) Sp3path(1, k) Sp1path(1, k);
                 Sp1path(2, k) Sp2path(2, k) Sp3path(2, k) Sp1path(2, k);
                 Sp1path(3, k) Sp2path(3, k) Sp3path(3, k) Sp1path(3, k)];
    
        plot3(Spcft(1, :), Spcft(2, :), Spcft(3, :), 'Color', [1 0 0], 'LineWidth', 3);

        if LISA_relative_circular_orbit
            plot3(C_total(1, :, k), C_total(2, :, k), C_total(3, :, k), 'Color', [1 1 1]);
        end

        view(az, el);
        title(title_text, 'Color', [1 1 1]);
        set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1], 'ZColor', [1 1 1]);
        drawnow;

        % Capture frame for GIF
        frame = getframe(h);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        if k == 1
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', time_lapse);
        else
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', time_lapse);
        end

        hold off;
    end

    %% F_+ and F_x calculation
    % Polar
    theta = 0:0.05:pi;
    % Azimuthal
    phi = 0:0.05:(2 * pi);

    [A, D] = meshgrid(phi, theta);
    X = sin(D) .* cos(A);
    Y = sin(D) .* sin(A);
    Z = cos(D);

    % Generate function values
    fPlus = zeros(length(theta), length(phi));
    fCross = zeros(length(theta), length(phi));
    for lp1 = 1:length(phi)
        for lp2 = 1:length(theta)
            [fPlus(lp2, lp1), fCross(lp2, lp1)] = detframefpfc(theta(lp2), phi(lp1));
        end
    end

    % Plot fPlus
    figure;
    surf(X, Y, Z, abs(fPlus));
    shading interp;
    axis equal;
    colorbar;
    title('F_+');

    % Plot fCross
    figure;
    surf(X, Y, Z, abs(fCross));
    shading interp;
    axis equal;
    colorbar;
    title('F_x');
end
