#By downloading, installing, or otherwise accessing or using the Software , you (“Recipient”) agree to receive and use the above-
#identified SOFTWARE subject to the following terms, obligations and restrictions. If you do not agree to all of the following terms, 
#obligations and restrictions you are not permitted to download, install,
#execute, access, or use the SOFTWARE:
#1.	Originators of the SOFTWARE.  Provider is willing to license its rights in the SOFTWARE (“Provider’s Rights”) to academic researchers to use free of charge solely for academic, non-commercial research purposes subject to the terms and conditions outlined herein. The SOFTWARE was created at the University of Wisconsin ("UW") by Dakai Tamada. Please note Provider's Rights may include, but are not limited to, certain patents or patent applications owned by the Wisconsin Alumni Research Foundation (“WARF”). 
#2.	Limited License.  Provider hereby grants to Recipient a non-commercial, non-transferable, royalty-free, non-exclusive license, without the right to sublicense, under Provider’s Rights to  download, install, access, execute and use the SOFTWARE solely for academic, non-commercial research purposes. SOFTWARE may not be used, directly or indirectly, to perform services for a fee or for the production or manufacture of products for sale to third parties. The foregoing license does not include any license to third party intellectual property that may be contained in the SOFTWARE; obtaining a license to such rights is Recipient’s responsibility. 
#3.	Restrictions on SOFTWARE use and distribution.  Recipient shall not take, authorize, or permit any of the following actions with the SOFTWARE: (1) modify, translate or otherwise create any derivative works; or (2) publicly display (e.g., Internet) or publicly perform (e.g., present at a press conference); or (3) sell, lease, rent or lend; or (4) use it for any commercial purposes whatsoever. Recipient must fully reproduce and not obscure, alter or remove any of the Provider’s proprietary notices that appear on the SOFTWARE, including copyright notices or additional license terms included with any the third party software contained in the SOFTWARE. Recipient may not provide any third party with access to the SOFTWARE or use the SOFTWARE on a timeshare or service bureau basis. Recipient represents that it is compliance with all applicable export control provisions and is not prohibited from receiving the SOFTWARE. 
#4.	Reservation of rights.  Provider retains all rights and title in the SOFTWARE, including without limitation all intellectual property rights (e.g., patent, copyright and trade secret rights) that may now or in the future exist in the SOFTWARE, regardless of form or medium. Provider retains ownership and all of Its rights in the SOFTWARE, including all of its intellectual property rights (e.g., patent, copyright and trade secret rights) that may now or in the future cover the SOFTWARE or any uses of the SOFTWARE, regardless of form or medium; title remains with Provider and the SOFTWARE is merely being loaned to Recipient for the specific purposes and under the specific restrictions stated herein. Nothing in this Agreement grants Recipient any additional rights to the SOFTWARE, any right to obtain any updates or new releases of the SOFTWARE, any commercial license for the SOFTWARE, or any other intellectual property owned or licensed by Provider. Provider has no obligation to provide any support, updates, or bug fixes.
#5.	Disclaimer of Warranty. PROVIDER IS PROVIDING THE SOFTWARE TO RECIPIENT ON AN “AS IS” BASIS. PROVIDER MAKES NO REPRESENTATIONS OR WARRANTIES CONCERNING THE SOFTWARE OR ANY OUTCOME THAT MAY BE OBTAINED BY USING THE SOFTWARE, AND EXPRESSLY DISCLAIMS ALL SUCH WARRANTIES, INCLUDING WITHOUT LIMITATION ANY EXPRESS OR IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT OF INTELLECTUAL PROPERTY RIGHTS. PROVIDER MAKES NO REMEDY THAT THE SOFTWARE WILL OPERATE ERROR FREE OR UNINTERRUPTED.
#6.	Limitation of Liability; Indemnity.  TO THE FULLEST EXTENT PERMITTED BY LAW, IN NO EVENT SHALL PROVIDER BE LIABLE TO RECIPIENT FOR ANY LOST PROFITS OR ANY DIRECT, INDIRECT, EXEMPLARY, PUNITIVE, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING FROM THE SOFTWARE OR ITS USE. FURTHERMORE, IN NO EVENT WILL PROVIDER’S LIABILITY TO RECIPIENT EXCEED $100. PROVIDER HAS NO LIABILITY FOR ANY DECISION, ACT OR OMISSION MADE BY RECIPIENT AS A RESULT OF USE OF THE SOFTWARE. To the extent permitted by applicable law, Recipient agrees to indemnify, defend and hold harmless Provider, UW, and the SOFTWARE authors against all claims and expenses, including legal expenses and reasonable attorneys fees, arising from Recipient’s use of the SOFTWARE.
#7.	No use of names/trademarks.  Recipient shall not use Provider’s name, or the name of any author of the SOFTWARE or that of UW, in any manner without the prior written approval of the entity or person whose name is being used.
#8.	Termination.  Without prejudice to any other rights, Provider may terminate this Agreement if Recipient fails to comply with the terms of this Agreement for any reason. Upon termination for any reason, Recipient must immediately destroy all copies of the SOFTWARE in Recipient’s possession, custody, or control.	

import scipy.io
from os import listdir
from os.path import isfile, join
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import norm
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from sklearn.linear_model import LinearRegression
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.font_manager as fm



# Figure t1 correction:
def figure_t1_correction(rc_parameters):
    # path to .mat file with results
    base_dir = 'T1_correction.mat'

    # read in matlab .mat file
    results = scipy.io.loadmat(base_dir, squeeze_me = True)

    # load results into dataframe
    df = pd.DataFrame()
    r2_pb = np.array(results['R2_pb'])
    r2_fit = np.array(results['R2_fit'])
    r1_all = np.array(results['R1_test'])
    pb_label = np.full(shape=len(r2_pb), fill_value='pb')
    fit_label = np.full(shape=len(r2_fit), fill_value='fit')
    df['R2 (s' + r'$^{-1}$' + ')'] = np.concatenate((r2_pb, r2_fit))
    df['R1 (s' + r'$^{-1}$' + ')'] = np.concatenate((r1_all, r1_all))
    df['Method'] = np.concatenate((pb_label, fit_label))

    # formatting
    sns.set_context("poster")
    sns.set(style="ticks", font_scale=3.5, rc=rc_parameters)

    # plot
    g = sns.lineplot(data = df, x = 'R1 (s' + r'$^{-1}$' + ')', y = 'R2 (s' + r'$^{-1}$' + ')', hue = 'Method', style = 'Method', linewidth = 20, palette = dict(pb='#9932cc', fit='#F57047'))
    plt.legend(loc = 'right', labels=[r'$R2_{PB}(R1)$', r'$R2_{fit}(R1)$'], edgecolor='inherit')
    # plt.title('T2 = 12 ms       T2 = 7 ms')
    # g.set(xlim=(0, 10))
    g.set(ylim=(100, 500))
    sns.despine()
    plt.subplots_adjust(left=0.15, bottom=0.15, right=0.7, top=0.88)
    plt.show()
    


# parameters needed to make figures have dark or light background and certain linewidth
def get_rc_parameters(dark_or_light, line_width):
    if dark_or_light == 'light':
        rc_parameters = {"lines.linewidth": line_width}
    else:
        rc_parameters = {
            "lines.color": "white",
            "patch.edgecolor": "white",
            "text.color": "white",
            "axes.facecolor": "black",
            "axes.edgecolor": "lightgray",
            "axes.labelcolor": "white",
            "xtick.color": "white",
            "ytick.color": "white", 
            "grid.color": "lightgray",
            "figure.facecolor": "black",
            "figure.edgecolor": "black",
            "savefig.facecolor": "black",
            "savefig.edgecolor": "black",
            "lines.linewidth": line_width}
    return rc_parameters



if __name__ == "__main__":

    # make all figures have dark or light background
    dark_or_light = 'light'
    rc_parameters = get_rc_parameters(dark_or_light, 7)

    figure_t1_correction(rc_parameters)
