function zlevels=confndlevelz(z,clevels,zlevelstep)
% z is probability (ie. exp(-loglike))
% zlevels are the conv levels in exp(-loglike)
% based on confndlevel

if (~exist('zlevelstep')) zlevelstep=500; end

z=reshape(z,prod(size(z)),[]); % irrelevant that is nd array

znorm=max(z);
z=z/znorm;

%% now get the values of -loglike for clevels (eg. 0.68) confidence

zlevels=zeros(size(clevels));

% integrate over whole range
norm=sum(z);
zmax=1;

% step through levels in z, evaluating confidence in each 
zlevel=zmax;
conf=0;
% I assume contour levels are increasing
for i=1:length(clevels)
    fprintf(1,['Finding contour level ',num2str(i),' ...\n'])
    while (conf<clevels(i)) % so it stops when conf is just greater than clevels(i)
        zlevel=zlevel-zmax/zlevelstep;
        %% make a vector to integerate which is the same as z, but zero where z.lt.lim
        tempmat=z;
        tempmat(find(z<zlevel))=0.;
        %% integrate it, normalise, compare with confidence level reqd
        conf=(sum(sum(tempmat)))/norm;
    end;
    
    % the rest is only important when things go wrong with the above
    while (conf==1.0)
        if (i>=2) zlevel=zlevels(i-1); end
        if (i==1) zlevel=zmax; end
        zlevelstep=zlevelstep*2;
        fprintf(1,['Increasing zlevelstep to ',num2str(zlevelstep),' in conf2dlevel to improve accuracy\n'])
        conf=0;
        while (conf<clevels(i))
            zlevel=zlevel-zmax/zlevelstep;
            %% make a vector to integerate which is the same as z, but zero where z.lt.lim
            tempmat=z;
            tempmat(find(z<zlevel))=0.;
            %% integrate it, normalise, compare with confidence level reqd
            conf=(sum(sum(tempmat)))/norm;
        end;
    end;

    fprintf(1,['Contour level ',num2str(i),' is actually the ',num2str(conf),' contour\n'])
    zlevels(i)=zlevel;
end;

if (any(zlevels<0)) 
    (['!!!! ERROR in conf2dlevel !!!! a zlevel is less than zero'])
end

%zlevel
%contour(x,y,z,[zlevel,zlevel])
%contour(xcoarse,ycoarse,zcoarse,[zlevel,zlevel])

zlevels=zlevels*znorm;





