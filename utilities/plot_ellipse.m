% plot ellipse x'*E*x <= r^2 around point cxy

% [Input]
%   r:        radius 
%   cxy:      center of circle
%   E:        2-by-2 symmetry matrix describing x'*E*x <= r^2
%   options:  color of line, 

% [Ref]
%   https://stackoverflow.com/questions/15915625/plotting-an-ellipse-in-matlab-given-in-matrix-form
    
function handle = plot_ellipse(r, cxy, E, line_color, line_style)
    switch nargin 
        case {1,2,3}
            line_color = [1 0 0];    % red
        
        case 4
            line_style = '-';

    end
    
    % find the Cholesky decomposition of E
    % x'*E*x = x'*(R'*R)*x = (Rx)'*(Rx) = z'*z
    R = chol(E);
    
    % sample on a circle
    N = 1e3;
    th = linspace(0, 2*pi, N);
    z = r*[cos(th); sin(th)];
    
    % transform the z back to original cordinate
    % R*x = z
    ellipse = R\z;
    
%     handle = patch(ellipse(1,:)+cxy(1),ellipse(2,:)+cxy(2),'red', ...
%         'FaceColor', color, 'FaceAlpha', 0.3, ...
%         'EdgeColor','none');
    handle = plot(ellipse(1,:)+cxy(1), ellipse(2,:)+cxy(2), ...
                    'color', line_color,'LineStyle',line_style);
end