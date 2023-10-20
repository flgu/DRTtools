function varargout = plot_DRT(x, y, options)

    arguments
        x
        y
        options.Color = "";
        options.LineWidth double = 2;
        options.LineStyle string = "-";
        options.DisplayName string = "DRT Spectrum";
    end % args

    fig = figure();
    ax = axes(fig);
    
    hold(ax, "on")
    plt = plot(ax, x, y);
    plt.LineWidth = options.LineWidth;
    plt.LineStyle = options.LineStyle;
    plt.DisplayName = options.DisplayName;
    
    hold(ax, "off")
    
    if strcmp(options.Color, "") == 0
        plt.Color = options.Color;
    end % if

    ax.XScale = 'log';
    
    grid(ax, "on")
    grid(ax, "minor")
    
    ax.XLabel.String = "f [Hz]";
    ax.YLabel.String = "\gamma(ln(f)) [\Omega]";

    % add optional outputs
    varargout{1} = ax;
    varargout{2} = fig;
end % fun def