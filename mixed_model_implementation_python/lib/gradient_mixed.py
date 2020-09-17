from . import mixmem_aux as aux
import time
import numpy as np

def gradient_descent(n, r0 = None, c0=None, N=1000,sig_digs=4, N_stable=600,
                     wr=1E-4, wc=1E-4,
                     up=1.2, down=0.8,
                     verbose_time = False, progressBar=False,
                     exp = lambda x: (1+x**2)**.5+x, nblocks = None):
    ''' 
    Optimize a mixed membership model using gradient descent.

    Parameters:

    -  n: Real number of trips between stations, indexed as [from, to, time_bin]
    - r0: Block strength initial guess, indexed as [block_1,block_2,time]
    - c0: Station block initial guess, indexed as [station, block]
    - wc: Initial guess for the gradient weight for c0
    - rc: Initial guess for the gradient weight for r0
    - Number of iterations and precision based stopping criterion:
        - sig_digs = # of significant digits to test for no change (default 4)
        - N_stable = # of iterations for the base-10 floating point expansion 
          unnormalized log-likelihood to be unchanged to consider the optimization 
          to be complete. (default 600)
          - N = the max number of iterations the algorithm will run, whether
          the stopping criterion is met in that many steps or not. (default 1000)
          - If sig_digs=None OR N_stable=None the algorithm will run for N steps.
            This is the same as ommitting a stopping criterion other than max steps
        

    - verbose: Whether the method should be verbose
    - progressBar: Whether a progressbar should be displayed
    
    - exp:  vectorized function to substitute the exponential in the logarithmic 
            gradient descent. Use at your own risk. The minimal requeriments on
            the function are:

            * f is positive
            * f(0) = 1
            * f is monotonous increasing, and strictly monotonous near 0.

            Numerical experiments indicate x -> (x^2+1)^1/2+x is a good function,
            but twice the inverse logit looks reasonable and way more robust. 
    - nblocks: Number of blocks. Only necessary if either r0 or c0 are not provided.
    
    - Returns the updated values of the parameters:
    
        (rfit,cfit, wrfit, wcfit)
        This allows for the following code:
        
        x = gradient_descent(100, r,c,wr,wc, **params)
        while (stuff(x)):
            x = gradient_descent(100,*x, **params)
            
        so that you can code your own stopping conditions
    '''
    #return first d significant digits
    round_to_d = lambda x, d: round(x, -int(np.floor(np.log10(np.abs(x)))) + (d - 1))
    
    
    
    
    #Set up the progress bar with ipywidgets
    if progressBar:
        import ipywidgets as widgets
        from IPython.display import display
        f = widgets.FloatProgress(min=0, max=N)
        display(f)
        
    #random initialisation
    if r0 is None:
        assert (nblocks is not None)
        ntimes = n.shape[2]
        r0 = np.random.lognormal(size=(nblocks,nblocks,ntimes))
    if c0 is None:
        assert (nblocks is not None)
        nstations = n.shape[0]
        c0 = np.random.lognormal(size=(nstations,nblocks))
        
    #Compute the original likelihood
    old_l=aux.likelihood(c0,r0,n);
    
    #Factors with which we will update the weights
    #so that the speed of the Gradient Descent is Optimal
    up = up
    down = down

    start_time = time.time()
    if verbose_time!=False: #print first line for the table
        print("N: \t C_RelChange \t R_RelChange \t log_like \t d_log_like \t time(s)")
        
    stable_counter=0
    for i in range(N):
        
        
        #Single gradient descent for the R
        #=================================
        
        rs = aux.r_score(c0,r0,n) #compute r gradient (from the mixmem_aux library)
        #Compute two new points with different log-gradient weights
        r1 = r0*exp(down*wr*rs);        r2 = r0*exp(up*wr*rs)
        r1 = np.maximum(r1, 1E-20);       r2 = np.maximum(r2, 1E-20);    
        l1 = aux.likelihood(c0,r1,n);   l2 = aux.likelihood(c0,r2,n)
        
        #choose (and keep) the best weight (if possible)
        #Return if we a wall of the only possibilities being nan or infinity
        if (l2>l1)&(not(np.isinf(l2)|np.isnan(l2))): r0=r2;  wr*=up;    l=l2;
        elif (l2<=l1)&(not(np.isinf(l1)|np.isnan(l1))): r0=r1;  wr*=down;  l=l1;
        else:  
            print(l1)
            return normalize(r0,c0)+(wr,wc);
        
        
       #########################
    
        #Copy of the code above for the C
        #================================
        
        cs = aux.c_score(c0,r0,n)
                 
        c1 = c0*exp(down*wc*cs);          c2 = c0*exp(up*wc*cs)
        c1 = np.maximum(c1, 1E-20);       c2 = np.maximum(c2, 1E-20); 
        l1 = aux.likelihood(c1,r0,n);     l2 = aux.likelihood(c2,r0,n)
        
        #choose (and keep) the best weight (if possible)
        #Return if we a wall of the only possibilities being nan or infinity
        if (l2>l1)&(not(np.isinf(l2)|np.isnan(l2))):  c0=c2;  wc*=up;  l=l2;
        elif (l2<=l1)&(not(np.isinf(l1)|np.isnan(l1))): c0=c1;  wc*=down;  l=l1;
        else:
            print((l1,l2))
            return normalize(r0,c0)+(wr,wc); 
        
        #update counter if floating point up to sig_digs precision does not change   
        if (sig_digs is not None)& (N_stable is not None):
            if round_to_d(l,sig_digs)==round_to_d(old_l,sig_digs):
                stable_counter+=1
            else:
                stable_counter=0
            
        #######################
        # If the user wants information to be printed every "verbose_time" steps, do so.
        if verbose_time!=False and i%verbose_time==0:
            deltar = np.mean(np.abs(exp(up*wr*rs)-1)*r0)/np.mean(r0)
            deltac = np.mean(np.abs(exp(up*wc*cs)-1)*c0)/np.mean(c0)
            print(("{}"+(4*"\t {:.4e}")+"\t {:.2f}")
                  .format(i,deltac, deltar,l, l-old_l, time.time()-start_time))
        
        #Update likelihood
        old_l = l;         
        if progressBar:
            f.value=i+1;
            
        #stop if l with sig_digs significant digits does not change for N_stable
        #steps
        if (sig_digs is not None)& (N_stable is not None):
            if stable_counter==N_stable:
                end_time = time.time()
                print(str((end_time-start_time)/60.0)+"Minutes")
                print("Stopping criterion met")
                return normalize(r0,c0)+(wr,wc);
    end_time = time.time()
    print(str((end_time-start_time)/60.0)+"Minutes")
    print("Stopping criterion not met")
    #end the progress bar
    if progressBar:
        f.close()
    
    #The set of minima is not unique: we can divide c by a constant
    #and multiply r by a constant and obtain the same model. We must
    #choose an (arbitrary) normalisation. 
    
    return normalize(r0,c0)+(wr,wc);


# Auxiliary normalizing function. Normalized so that the
# station blocks sum to 1. 
def normalize(r,c):
    cnorm = np.sum(c,axis=0)
    rnorm = np.outer(cnorm,cnorm)
    cn = c/cnorm[np.newaxis,:]
    rn = r*rnorm[:,:,np.newaxis]
    return rn,cn

# Function for initializing the parameters so that they have magnitudes close to the magnitude of the data. 
def realistic_initial(n,N,numblocks,T=24):
    r0 = np.exp(np.random.normal(size=(numblocks,numblocks,T)))#**2+1
    c0 = np.exp(np.random.normal(size=(N,numblocks)))#**2+1
    c0=c0/np.sum(c0,0,keepdims=True)
    r0=r0*float(np.sum(n))/(float(numblocks**2)*T)
    return (r0,c0)
