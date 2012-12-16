classdef getdist_defaults < handle
    properties
        plot_size_inch = 6;
        lab_fontsize=12;
        axes_fontsize=12;
        likes=true;
        colormap='jet';
        lineM = {'-k','-r','-b','-m','-g','-c','-y','--k','--r','--b'};
        lineL = {':k',':r',':b',':m',':g',':c',':y','-.k','-.r','-.b'};
        colstr='krbmgcykrb';
        lw1=1;
        lw2=1;
        printex='pdf';
        print='-dpdf';
        last;
    end
    methods
        function obj = getdist_defaults(plot_data) % constructor
        if nargin > 0 && ~isempty(plot_data)
            path(plot_data,path)
        else
            p=getenv('getdist_plot_data');
            if p(length(p))=='/'
             p=p(1:length(p)-1);
            end
            if ~isempty(p) && isempty(strfind(path,p))
             path(p,path)
            end
        end
        end
    end
end



