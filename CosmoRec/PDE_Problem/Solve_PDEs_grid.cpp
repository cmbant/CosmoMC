//====================================================================================================================
// multiple resonances
//====================================================================================================================
void init_PDE_xarr(int npts, double *xarr, double xmin, double xmax, double xres, double xD, int mess=0)
{
    //===================================================================================   
    xD*=xres;
    vector<double> Dxarr(npts);
    
    //===================================================================================
    // lowest region
    //===================================================================================
    double core_width=10.0;
    int ncore=500;
    double enhance_fac=2.0;
    //
    int nlow_high=(npts-ncore);
    
    if(nlow_high<0){ cout << " init_PDE_xarr:: Need more points per resonance. Add about " << nlow_high << " more points" << endl; exit(0); }
    
    if(mess==1) cout << " xres= " << xres << endl;
    
    //===================================================================================
    double dec_up=log10(xmax-xres)-log10(core_width*xD);
    double dec_tot=log10(xres-xmin)-log10(core_width*xD) + enhance_fac*dec_up;
    double w_up=enhance_fac*dec_up/dec_tot;
    //===================================================================================   
    int nup=npts*w_up-1-ncore, nlow=npts-nup-ncore;
    
    //===================================================================================
    init_xarr(core_width*xD, xmax-xres, &xarr[nlow+ncore], nup, 1, mess);
    for(int k=0; k<nup; k++) xarr[nlow+ncore+k]+=xres;
    //===================================================================================
    init_xarr(xres-core_width*xD, xres+core_width*xD, &xarr[nlow], ncore+1, 0, mess);
    //===================================================================================
    init_xarr(core_width*xD, xres-xmin, &Dxarr[0], nlow+1, 1, mess);
    for(int k=0; k<nlow+1; k++) xarr[k]=xres-Dxarr[nlow-k];
    
    for(int k=1; k<npts; k++) if(xarr[k-1]>xarr[k]){ cerr << " init_PDE_xarr:: grid non-monotonic " << k << " " << xarr[k] << endl; exit(0); }
    
    return;
}

void init_PDE_xarr_cores(int npts, int npts_res, vector<double> &xarr, double xmin, double xmax, vector<double> &resonances)
{
    if(show_messages>=1) cout << "\n init_PDE_xarr_cores:: setting up grid. " << endl;
    
    int mess=(show_messages>=2 ? 1 : 0);
    int nres=resonances.size();
    int np_2s=2000; 
    double xcrit=0.8;
    
    if(nres*npts_res+np_2s>npts){ cerr << " init_PDE_xarr_cores:: Need more points per resonance. Add about " << nres*npts_res+np_2s-npts << " more points" << endl; exit(0); }
    
    vector<double> dx(nres);
    dx[0]=xmax-xmin;
    //===========================================
    // define regions around resonances
    //===========================================
    if(nres>1) dx[0]=(resonances[0]+resonances[1])*0.5-xmin;
    for(int i=1; i<nres-1; i++) dx[i]=(resonances[i+1]-resonances[i-1])*0.5;
    if(nres>1) dx[nres-1]=xmax-(resonances[nres-2]+resonances[nres-1])*0.5;
    
    //===========================================
    // below Ly-a (2s-1s part)
    //===========================================
    double xl=xmin;
    int istart=0;
    //
    if(xmin<=xcrit)
    {
        //===========================================
        // normal log-scale for x<xcrit
        //===========================================
        init_xarr(xl, xcrit, &xarr[istart], np_2s, 1, mess);    
        istart+=np_2s-1;
        xl=xcrit;
        
        dx[0]+=xmin-xl;
    }
    
    if(nres>1) 
    {
        //===========================================
        // Ly-a
        //===========================================
        init_PDE_xarr(npts_res, &xarr[istart], xl, xl+dx[0], 1.0, 2.0e-5, mess);
        istart+=npts_res-1;
        
        //===========================================
        // Ly-n
        //===========================================
        for(int i=1; i<nres-1; i++) 
        {
            xl+=dx[i-1];
            init_PDE_xarr(npts_res, &xarr[istart], xl, xl+dx[i], resonances[i], 2.0e-5, mess);
            
            istart+=npts_res-1;
        }
        
        xl+=dx[nres-2];
        init_PDE_xarr(npts-istart, &xarr[istart], xl, xl+dx[nres-1], resonances[nres-1], 2.0e-5, mess);
    }
    else 
    {
        //===========================================
        // Ly-a
        //===========================================
        init_PDE_xarr(npts-istart, &xarr[istart], xl, xl+dx[0], 1.0, 2.0e-5, mess);
    }
    
    for(int k=1; k<npts; k++) if(xarr[k-1]>xarr[k]){ cerr << " init_PDE_xarr_cores::grid non-monotonic " << k << " " << xarr[k] << endl; exit(0); }
    
    if(show_messages>=1) cout << " init_PDE_xarr_cores:: done. " << endl;

    return;
}   



