
def Matrix_plot(data, names, bg_colors = None,params = {}, ylim = None, x_values = None, xlim = None):

    '''
     @param data: matrix of functions to plot
    
     @param names: Data names, array
    
     @param params: Dictionary, does nothing right now,
                    in the future it should pass parameters
                    to plot and fill_between
    
     @param ylim: does what expected. Should be unified with params
    
     @param x_values: x values, should be the same for all data right now
                     (one could improve that)
    
     @param xlim: x limits to be used, also should go to params
    
     @returns: the created figure and axis
    
     @side-effects: creates the plot
     '''
    
    import matplotlib.pyplot as plt
    import numpy as np

    
    data_t = data.swapaxes(0,1)

    n_plots = data.shape[0]
    figure, axis = plt.subplots(n_plots, n_plots,figsize =(15,15), sharex='col', sharey='row')

    if x_values is None:
        x_values = range(data.shape[2])

    #Tight X-lim if it doesn't exist
    if xlim is None:
        xlim = (np.min(x_values), np.max(x_values))

    if ylim is None:
        ylim = np.max(np.abs(data))
        ylim = (-ylim,ylim)

    #Set transparent bg_color if no bg_color exists
    if bg_colors is None:
        bg_colors = [[(1,1,1,0)]*len(data)]*len(data)


    for row in zip(data,data_t,axis,bg_colors):
        for curve, curve_t, a, bgcolor in zip(*row):

            a.fill_between(x_values,0,   curve, color='k', 
                          linewidth=0,edgecolor='green')
            a.fill_between(x_values,-curve_t,0, color='k', 
                          linewidth=0,edgecolor='green')
            
            a.plot(x_values, [0]*len(x_values), linewidth=.5, c='k')


            a.set_xlim(xlim)
            a.set_ylim(ylim)

            a.set_facecolor(bgcolor)

            a.tick_params(axis='x', colors=(0,0,0,0))
            a.tick_params(axis='y', colors=(0,0,0,0))

            for x in ('top', 'bottom','left','right'):
                a.spines[x].set_visible(False)

    for a,n in zip(axis[-1],names):
        a.set_xlabel(n.replace('(', '\n('), rotation = 45, ha = 'right')
    for a,n in zip(axis,names):
        a = a[0]
        a.tick_params(axis='y', which='major', pad=-15)
        a.set_ylabel(n.replace('(', '\n('), rotation = 45, va = 'bottom', ha = 'right')

    figure.subplots_adjust(hspace=0,wspace = 0)

    return (figure,axis)