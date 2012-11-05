classdef getdist_defaults < handle
    properties
        plot_size_inch = 6;
        lab_fontsize=12;
        axes_fontsize=12;
        likes=true;
        colormap='jet';
        lineM = {'-k','-r','-b','-m','-g','-c','-y','--k','--r','--b'};
        lineL = {':k',':r',':b',':m',':g',':c',':y','-.k','-.r','-.b'};
        lw1=1;
        lw2=1;
        printex='pdf';
        print='-dpdf';
        last;
    end
    methods
        function obj = getdist_defaults(plot_data) % constructor
        if nargin > 0 && path ~=''
            path(plot_data,path)
        else
            p=getenv('getdist_plot_data');
            if  size(p,1)>0
             path(p,path)
            end
        end
        end
    end
end



