import numpy as np
import pandas as pd

import scipy.stats

import matplotlib.pyplot as plt

def rapid_linear_reg(x, y, plot=False):
    ox = np.array(x)
    oy = np.array(y)

    not_nan = ~(np.isnan(ox)|np.isnan(oy))
    x, y = ox[not_nan], oy[not_nan]

    reg = scipy.stats.linregress(x,y)
    print("The linear model is: Y = {:.5} * X + {:.5}".format(reg.slope, reg.intercept))

    pearsonr = scipy.stats.pearsonr(x,y)
    spearmanr = scipy.stats.spearmanr(x,y)
    
    print(pearsonr)
    print(spearmanr)
    py = reg.slope * ox + reg.intercept

    if plot:
        plt.scatter(ox,oy)
        plt.plot(ox, py, c = 'red')
        plt.show()

    return py, pearsonr, spearmanr


def rapid_process_result(loss, r_df, path=True, plot=False):
    if path:
        loss = np.load(loss)
        r_df = pd.read_csv(r_df, index_col=0)
    print('epochs: {0}'.format(loss.shape[1]))

    r = r_df['r'].to_numpy().reshape(-1,1)
    p = r_df.to_numpy()[:,1:]
    ae = np.abs(p-r)

    df = pd.DataFrame(index=r_df.index, columns=['r','p','ae'], dtype=float)
    df['r'] = r_df['r']

    for i in range(loss.shape[0]):
        val_idx = list(range(ae.shape[0]))[i::loss.shape[0]]
        r_i = r[val_idx, :]
        p_i = p[val_idx, :]
        ae_i = np.abs(p_i-r_i)

        Min_loss_i = [np.where(loss[i,:] == np.min(loss[i,:]))[0][0]]
        Min_val_i = [np.where(np.mean(ae_i**2, axis=0) == np.min(np.mean(ae_i**2, axis=0)))[0][0]]
        
        df.loc[val_idx, 'p'] = p_i[:,Min_loss_i]
    
    final_p = df['p'].to_numpy().reshape(-1,1)
    final_ae = abs(final_p - r)

    print('Median Absolute Error:', np.median(final_ae))
    print('Mean Absolute Error:', np.mean(final_ae))

    Min_loss = np.argmin(np.mean(loss, axis=0))
    Min_mae = np.argmin(np.median(ae, axis=0))
    print('Min_loss:', Min_loss, np.mean(loss, axis=0)[Min_loss])

    if plot:
        plt.figure(dpi=100,figsize = (12,6))
        plt.plot(np.mean(loss, axis=0))
        plt.plot(np.mean(ae, axis=0))
        plt.plot(np.median(ae, axis=0))
        plt.scatter(Min_loss, np.mean(loss, axis=0)[Min_loss], c='r', zorder = 2)
        plt.scatter(np.argmin(np.mean(ae, axis=0)), np.min(np.mean(ae, axis=0)), c='r', zorder = 5)
        plt.scatter(Min_mae, np.median(ae, axis=0)[Min_mae], c='r', zorder = 2)

        plt.axhline(y=np.median(final_ae), linewidth=1,color='grey')
        plt.ylim(0, 10)
        plt.show()

    return r, final_p, final_ae