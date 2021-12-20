#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

import matplotlib, numpy as np
import matplotlib.pyplot as plt
def get_meshgrid_function_data(function, plot_limits, num_pts_1d, qoi = 0):
    """
    Generate data from a function in the format needed for plotting.
    Samples are generated between spcified lower and upper bounds
    and the function is evaluated on the tensor product of the 1d samples.

    Parameters
    ----------
    function : callable function
        The function must accept an np.ndarray of size (2, num_pts_1d**2)
        and return a np.ndarray of size (num_pts_1d,num_qoi)

    plot_limits : np.ndarray or list
        The lower and upper bounds used to sample the function
        [lb1,ub1,lb1,lb2]

    num_pts_1d : integer
        The number of samples in each dimension. The function is evaluated
        on the tensor product of the 1d samples

    qoi : integer
        function returns a np.ndarray of size (num_pts_1d,num_qoi) qoi
        specifies which column of the array to access.

    Returns
    -------
    X : np.ndarray of size (num_pts_1d,num_pts_1d)
        The 1st coordinate of the samples

    Y : np.ndarray of size (num_pts_1d,num_pts_1d)
        The 2nd coordinate of the samples

    Z : np.ndarray of size (num_pts_1d,num_pts_1d)
        The function values at each sample
    """
    x = np.linspace( plot_limits[0], plot_limits[1], num_pts_1d )
    y = np.linspace( plot_limits[2], plot_limits[3], num_pts_1d )
    X, Y = np.meshgrid( x, y )
    pts = np.vstack( ( X.reshape( ( 1, X.shape[0]*X.shape[1] ) ),
                          Y.reshape( ( 1, Y.shape[0]*Y.shape[1]) ) ) )

    Z = function( pts )
    if ( Z.ndim == 2 ):
        Z = Z[:,qoi]
    Z = np.reshape( Z, ( X.shape[0], X.shape[1]) )
    return X,Y,Z

def plot_contours(X, Y, Z, num_contour_levels=10, ax=None, offset=0, 
                  cmap=matplotlib.cm.jet, zorder=None):
    """
    Plot the contours of a two-dimensional function using matplotlib.contourf.

    Parameters
    ----------
    num_contour_levels : boolean (default=10)
       Plot the contours of the function to plot beneath the 3d surface.
       
    offset : boolean (default=0)
       Plot the contours offset from z=0 plane

    X,Y,Z: see documentation of Return of function get_meshgrid_function_data
    """
    # Plot contours of functions underneath surface
    if num_contour_levels>0:
        cset = ax.contourf(
            X, Y, Z, zdir='z', offset=offset,
            levels=np.linspace(Z.min(),Z.max(),num_contour_levels),
            cmap=cmap, zorder=zorder)

    return ax


def plot_surface(X, Y, Z, ax=None, samples=None, limit_state=None,
                 num_contour_levels=0, plot_axes=True, cmap=matplotlib.cm.jet,
                 axis_labels=None, angle=None, alpha=1., zorder=None):
    """
    Plot the three-dimensional surface of a two-dimensional function using
    matplotlib.plot_surface

    Parameters
    ----------
    samples : np.ndarray (num_vars+1 x num_samples)
        The x,y,z coordinates of samples to plot

    limit_state : callable_function
        Function must take samples np.ndarray (num_vars x num_samples)
        and values np.ndarray (num_samples), and return the indices
        a set of indices which correspond to samples within the limit state
        and the valuethat should be assigned to each of those samples , e.g.
        e.g. I, default_value = limit_state(samples,values)

    plot_axes : boolean (default=True)
       Plot the 3 coordinate axes and ticks

    angle : integer (default=None)
       The angle at which the plot is viewed. If None matplotlib will
       use its default angle

    axis_labels : list of strings size (3)
       Labels of the x, y, and z axes

    alpha : double
       Transperancy of plot

    zorder : the priority of the plot on the axes. A zorder=2 places plot
    above all plots with zorder=1

    X,Y,Z: see documentation of Return of function get_meshgrid_function_data

    """
    if ax is None:
        fig = plt.figure()
        try:
            ax = fig.gca(projection='3d')
        except ValueError :
            # Add following import to avoid error Unknown projection '3d'
            # when using ax = fig.gca(projection='3d')
            from mpl_toolkits.mplot3d import Axes3D
            ax = Axes3D(fig)

    # Define transperancy of plot
    if samples is not None:
        ax.scatter3D(samples[0,:], samples[1,:], samples[2,:], marker='o',
                     s=100, color='k')

    if limit_state is not None:
        pts = np.vstack((X.reshape((1, X.shape[0]*X.shape[1])),
                         Y.reshape((1, Y.shape[0]*Y.shape[1]))))
        vals = Z.flatten()
        I, default_value = limit_state(pts, vals)
        vals[I] = default_value
        Z_tmp = vals.reshape((X.shape[0],X.shape[1]))
    else:
        Z_tmp = Z

    # Plot surface
    ax.plot_surface(X, Y, Z_tmp, cmap=cmap, antialiased=False,
                    cstride=1, rstride=1, linewidth=0., alpha=alpha,
                    zorder=zorder)

    if not plot_axes:
        ax.set_xticks([]) 
        ax.set_yticks([]) 
        ax.set_zticks([])
        # Get rid of the panes
        ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

        # Get rid of the spines
        ax.w_xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
        ax.w_yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
        ax.w_zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))

    if ( axis_labels is not None ):
        ax.set_xlabel(axis_labels[0])
        ax.set_ylabel(axis_labels[1])
        ax.set_zlabel(axis_labels[2])

    # Set the angle of the plot (how it is viewed)
    if ( angle is not None ):
        ax.view_init(ax.elev, angle)

    return ax

