function [samples, names]=readsampsf(filescell,namesfile,nburn,nthin)
%READSAMPSF Read in samples files.
%
%   [samples,NAMES] = READSAMPS(FILESCELL) reads in the samples in 
%   file  into array samples. It also creates 
%   an array NAMES which corresponds to the parameter names in the default
%   cosmomc implementation.
%

% copied from READSAMPS.M slb4apr04 so just give filename

%%if (~exist('ncols')) ncols=21; end % for pre 2004 cosmomc
%if (~exist('ncols')) ncols=22; end
if (~exist('nburn')) nburn=0; end
if (~exist('nthin')) nthin=1; end

nfiles=length(filescell);
fprintf(1,'\n');

% set up default answers in case there is a problem
samples=[];
names=[];

for ifile=1:nfiles
    filename=filescell{ifile};
    fid=fopen(filename);
    if (fid>1) 
        
        fprintf(1,'\nSucessfully opened file %i of %i: %s\n',ifile,nfiles,filename)
        
        % find out how many columns there are then restart
        tmp=fgetl(fid);
        tmp1=sscanf(tmp,'%g ');
        ncols=length(tmp1);
        fclose(fid);
        fid=fopen(filename);
        
        samplestmp=fscanf(fid,'%g');
        fclose(fid);
        nsamp=length(samplestmp)/ncols;

        % check it is an integer. If not then this is a problem.
        if ~(round(nsamp)==nsamp)
            fprintf(1,'!!! ERRORR !!!\n')
            fprintf(1,'Number of file entries is not a multiple of %i\n',...
                ncols);
            fprintf(1,'Number of file entries = %i\n',length(samplestmp));
            fprintf(1,'Ignoring incomplete rows for now\n')
%            fprintf(1,'Check ncols is correct!\n')
            nsamp=floor(nsamp);
            return
        end

        samplestmp=samplestmp(1:(nsamp*ncols));
        samplestmp=reshape(samplestmp,[ncols,nsamp]);
        % remove burn-in
        is=(nburn+1):nthin:nsamp;
        samplestmp=samplestmp(:,is);
        fprintf(1,'Read %i samples. ',nsamp)
        fprintf(1,'Removed the first %i for burn-in.\n',nburn)
        if (nthin>1) fprintf(1,'Thinned by a factor %i\n',nthin); end
        fprintf(1,'%i samples remaining\n',length(is));
        if (length(is)==0)
            fprintf(1,'\n!!! You have %i samples in file %s, and %i samples are being burned.\n',length(is),filename,nburn)
            fprintf(1,'Either get more samples or change No. to Burn (see text box on the load samples gui)');
            fprintf(1,' if you would like to use samples from %s\n',filename)
        end
        
        % append this to the full set of samples
        if (ifile==1)
            samples=samplestmp;
        else
            sz_samples=size(samples)
            sz_samplestmp=size(samplestmp)
            isampled=find(min(samples,[],2)-max(samples,[],2));
            isampledtmp=find(min(samplestmp,[],2)-max(samplestmp,[],2));
            if (length(isampled)==length(isampledtmp))
                if (all(isampled==isampledtmp)&(sz_samples(1)==sz_samplestmp(1)))
                    samples=[samples samplestmp];
                else
                    samples=[];
                    fprintf(1,'ERROR! Samples files are not compatible\n')
                    fprintf(1,'Different parameters sampled in each\n')
                    return
                end
            else
                samples=[];
                fprintf(1,'ERROR! Samples files are not compatible\n')
                fprintf(1,'Different numbers of parameters are sampled in each\n')
                return
            end
        end
        
    else
        fprintf(1,['Couldnt open file:',filename]);
    end
end
[ncolstmp totsamps]=size(samples);
fprintf(1,'Kept %i samples in total\n\n',totsamps);

if (exist('namesfile'))
    fid=fopen(namesfile);
    if (fid>0)
        fclose(fid);
    names=textread(namesfile,'%s','delimiter','\n') ;
    else
        fprintf(1,'ERROR: %s not found\n',namesfile);
        return;
    end
else
    for i=1:ncols
        if (i<10)
            nameschar(i,:)=['Column  ',num2str(i)];
        elseif (i>=100)
            fprintf('Need to modify readsampsf.m for ncols>=100\n')
        else
            nameschar(i,:)=['Column ',num2str(i)]; 
        end
    end
    names=cellstr(nameschar);
end
