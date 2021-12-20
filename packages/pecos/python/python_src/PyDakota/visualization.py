#  _______________________________________________________________________
#
#  PECOS: Parallel Environment for Creation Of Stochastics
#  Copyright (c) 2011, Sandia National Laboratories.
#  This software is distributed under the GNU Lesser General Public License.
#  For more information, see the README file in the top Pecos directory.
#  _______________________________________________________________________

import matplotlib.pyplot as plt, numpy
def cross_validated_degree_bar_plot(cv_score_per_degree, poly_degrees):
    num_rhs = cv_score_per_degree.shape[1]
    f,axs=plt.subplots(1,num_rhs,sharey=True,figsize=(num_rhs*8, 6))
    axs = axs.ravel()
    for i in range(num_rhs):
        data = cv_score_per_degree[:,i]
        rects = axs[i].bar(poly_degrees, data, width=0.7,log=True)
        best_cv_score_index=numpy.argmin(data)
        rect=rects[best_cv_score_index]
        rect.set_color('r')
        axs[i].text(rect.get_x() + rect.get_width()/2., 1.05*rect.get_height(),
                    '%d' % poly_degrees[best_cv_score_index],
                    ha='center', va='bottom')
        axs[i].set_title('QoI %d'%(i+1))
        axs[i].set_xlabel('degree')
        axs[i].set_ylabel('CV error')
    # share y uses axis of left most plot which may not have the greatest range
    # so set range here
    plt.ylim(10**numpy.floor(numpy.log10(cv_score_per_degree.min())),
             10**numpy.ceil(numpy.log10(cv_score_per_degree.max())))
