import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

sigfigs = lambda value, n_sigfigs: np.round(value, -int(np.ceil(np.sign(value) * np.log10(abs(value)))) + n_sigfigs)

class Metrics:

    @classmethod
    def pearson_cor_coef(cls, X, Y):

        return np.corrcoef(X, Y)[0, 1]

    @classmethod
    def spearman_cor_coef(cls, X, Y):

        return np.corrcoef(X.argsort().argsort(), Y.argsort().argsort())[0, 1]

    @classmethod
    def concordance_cor_coef(cls, X, Y):

        _X = X / np.sqrt(np.sum(X ** 2))
        _Y = Y / np.sqrt(np.sum(Y ** 2))

        return 2. * cls.pearson_cor_coef(_X, _Y) * np.std(_X) * np.std(_Y) / (np.std(_X) ** 2 + np.std(_Y) ** 2 + (np.mean(_X) - np.mean(_Y)) ** 2)

    @classmethod
    def root_mean_square_error(cls, X, Y):

        return np.sqrt(np.mean((X / np.sqrt(np.sum(X ** 2)) - Y / np.sqrt(np.sum(Y ** 2))) ** 2))

def calculate_correlation_make_plots(beta, betags, cell_types, name, precision = 3, noPlot = False,
                                     labelEachSample = False, sampleLabels = None, extension = '.png', dpi = 300):

    figWidth = 8.
    figHeight = 10.

    fig = plt.figure(figsize=(figWidth, figHeight))
    plt.suptitle(name, fontsize=15)

    grid = matplotlib.gridspec.GridSpec(5 if len(cell_types) > 10 else 4,
                                        3 if len(cell_types) > 10 else 2,
                                        hspace=0.5, wspace=0.8,
                                        bottom=0.5 / figHeight, top=1. - 0.8 / figHeight,
                                        left=0.5 / figWidth, right=1. - 1. / figWidth)

    sigfigs = lambda value, n_sigfigs: np.round(value, -int(np.ceil(np.sign(value) * np.log10(abs(value)))) + n_sigfigs) if (value == value) and (value != 0.) else 0.


    Pcoeffs, Scoeffs = [], []

    for j in range(len(cell_types)):
        wh = np.where(~np.isnan(betags[:,j]))[0]
        Pcoeffs.append(sigfigs(Metrics.pearson_cor_coef(beta[:,j][wh], betags[:,j][wh]), precision))
        Scoeffs.append(sigfigs(Metrics.spearman_cor_coef(beta[:,j][wh], betags[:,j][wh]), precision))

    df_result = pd.DataFrame(data=np.vstack([Pcoeffs, Scoeffs]),
                             columns = cell_types,
                             index = ['Pearson', 'Spearman'])

    if noPlot:
        return df_result

    for j in range(len(cell_types)):
        ax1 = plt.subplot(grid[j], label="1")
        ax1.plot(beta[:,j], betags[:,j], 'o', markersize=1.5, mec='k', mew=0.5, alpha=0.6)

        if labelEachSample and (not sampleLabels is None):
            for i in range(len(beta[:,j])):
                ax1.text(beta[i,j], betags[i,j], sampleLabels, va='top', fontsize=2)

        if np.isnan(betags[:,j]).all():

            ax1.set_xticks([])
            ax1.set_xticklabels([])
            ax1.set_yticks([])
            ax1.set_yticklabels([])

            df_result.iloc[:,j] = np.nan

            continue

        cell_types[j], Pcoeffs[j]
        cell_types[j], Scoeffs[j]


        beta_sub = beta[:,j]
        betags_sub = betags[:,j]

        ax1.set_xlim([np.min(beta_sub[~np.isnan(beta_sub)]) if len(beta_sub) > 0 else 0, np.max(beta_sub[~np.isnan(beta_sub)]) if len(beta_sub) > 0 else 1.0])
        ax1.set_ylim([np.min(betags_sub[~np.isnan(betags_sub)]) if len(betags_sub) > 0 else 0, np.max(betags_sub[~np.isnan(betags_sub)]) if len(betags_sub) > 0 else 1.0])
        ax1.tick_params(labelsize=4)
        ax1.tick_params(axis='x', labelrotation=90, colors='black')

        ax2 = fig.add_subplot(grid[j], label="2", frame_on=False)

        ax2.plot(beta[:,j].argsort().argsort(), betags[:,j].argsort().argsort(), 'o', markersize=1.5, mec='orange', mew=0.5, alpha=0.6)
        ax2.set_xlim([0.0, sigfigs(np.max(beta[:,j].argsort().argsort()), 1)])
        ax2.set_ylim([0.0, np.trunc(1.25 * np.max(betags[:,j].T.argsort().argsort()) * 10.) / 10.])
        ax2.tick_params(labelsize=4, colors='orange')
        ticks = np.round(np.linspace(0, sigfigs(np.max(beta[:,j].argsort().argsort()), 1), num=5, endpoint=True), 5).tolist()
        ticks += [np.round(ticks[-1] + ticks[1] - ticks[0], 5)]
        ax2.set_xticks(ticks)
        ax2.set_xticklabels(ticks)
        ax2.xaxis.tick_top()
        ax2.yaxis.tick_right()
        ax2.xaxis.set_label_position('top') 
        ax2.yaxis.set_label_position('right') 
        ax2.tick_params(axis='x', labelrotation=90)

        ax2.text(1.125, 1., 'Pearson\n%s' % (Pcoeffs[j]), transform=ax1.transAxes, va='top', fontsize=7, color='k', weight='semibold')
        ax2.text(1.125, 0.8, 'Spearman\n%s' % (Scoeffs[j]), transform=ax1.transAxes, va='top', fontsize=7, color='orange', weight='semibold')

        ax2.text(0.5, 0.95, ' '.join(cell_types[j].split(' ')[:3]), transform=ax1.transAxes, va='top', ha='center', fontsize=10, color='b')

    fig.savefig(name + extension, dpi=dpi)

    plt.close(fig)

    return df_result



