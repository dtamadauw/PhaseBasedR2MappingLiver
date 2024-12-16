%TITLE:MRI simulation code and phantom data sets
%Matlab script to generate figures using Sobol's GRE simulation
%Author: Daiki Tamada
%Affiliation: Department of Radiology, University of Wisconsin-Madison
%Date: 10/28/2024
%Email: dtamada@wisc.edu

%By downloading, installing, or otherwise accessing or using the Software , you (“Recipient”) agree to receive and use the above-
%identified SOFTWARE subject to the following terms, obligations and restrictions. If you do not agree to all of the following terms, 
%obligations and restrictions you are not permitted to download, install,
%execute, access, or use the SOFTWARE:

%1.	Originators of the SOFTWARE.  Provider is willing to license its rights in the SOFTWARE (“Provider’s Rights”) to academic researchers to use free of charge solely for academic, non-commercial research purposes subject to the terms and conditions outlined herein. The SOFTWARE was created at the University of Wisconsin ("UW") by Dakai Tamada. Please note Provider's Rights may include, but are not limited to, certain patents or patent applications owned by the Wisconsin Alumni Research Foundation (“WARF”). 
%2.	Limited License.  Provider hereby grants to Recipient a non-commercial, non-transferable, royalty-free, non-exclusive license, without the right to sublicense, under Provider’s Rights to  download, install, access, execute and use the SOFTWARE solely for academic, non-commercial research purposes. SOFTWARE may not be used, directly or indirectly, to perform services for a fee or for the production or manufacture of products for sale to third parties. The foregoing license does not include any license to third party intellectual property that may be contained in the SOFTWARE; obtaining a license to such rights is Recipient’s responsibility. 
%3.	Restrictions on SOFTWARE use and distribution.  Recipient shall not take, authorize, or permit any of the following actions with the SOFTWARE: (1) modify, translate or otherwise create any derivative works; or (2) publicly display (e.g., Internet) or publicly perform (e.g., present at a press conference); or (3) sell, lease, rent or lend; or (4) use it for any commercial purposes whatsoever. Recipient must fully reproduce and not obscure, alter or remove any of the Provider’s proprietary notices that appear on the SOFTWARE, including copyright notices or additional license terms included with any the third party software contained in the SOFTWARE. Recipient may not provide any third party with access to the SOFTWARE or use the SOFTWARE on a timeshare or service bureau basis. Recipient represents that it is compliance with all applicable export control provisions and is not prohibited from receiving the SOFTWARE. 
%4.	Reservation of rights.  Provider retains all rights and title in the SOFTWARE, including without limitation all intellectual property rights (e.g., patent, copyright and trade secret rights) that may now or in the future exist in the SOFTWARE, regardless of form or medium. Provider retains ownership and all of Its rights in the SOFTWARE, including all of its intellectual property rights (e.g., patent, copyright and trade secret rights) that may now or in the future cover the SOFTWARE or any uses of the SOFTWARE, regardless of form or medium; title remains with Provider and the SOFTWARE is merely being loaned to Recipient for the specific purposes and under the specific restrictions stated herein. Nothing in this Agreement grants Recipient any additional rights to the SOFTWARE, any right to obtain any updates or new releases of the SOFTWARE, any commercial license for the SOFTWARE, or any other intellectual property owned or licensed by Provider. Provider has no obligation to provide any support, updates, or bug fixes.
%5.	Disclaimer of Warranty. PROVIDER IS PROVIDING THE SOFTWARE TO RECIPIENT ON AN “AS IS” BASIS. PROVIDER MAKES NO REPRESENTATIONS OR WARRANTIES CONCERNING THE SOFTWARE OR ANY OUTCOME THAT MAY BE OBTAINED BY USING THE SOFTWARE, AND EXPRESSLY DISCLAIMS ALL SUCH WARRANTIES, INCLUDING WITHOUT LIMITATION ANY EXPRESS OR IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT OF INTELLECTUAL PROPERTY RIGHTS. PROVIDER MAKES NO REMEDY THAT THE SOFTWARE WILL OPERATE ERROR FREE OR UNINTERRUPTED.
%6.	Limitation of Liability; Indemnity.  TO THE FULLEST EXTENT PERMITTED BY LAW, IN NO EVENT SHALL PROVIDER BE LIABLE TO RECIPIENT FOR ANY LOST PROFITS OR ANY DIRECT, INDIRECT, EXEMPLARY, PUNITIVE, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING FROM THE SOFTWARE OR ITS USE. FURTHERMORE, IN NO EVENT WILL PROVIDER’S LIABILITY TO RECIPIENT EXCEED $100. PROVIDER HAS NO LIABILITY FOR ANY DECISION, ACT OR OMISSION MADE BY RECIPIENT AS A RESULT OF USE OF THE SOFTWARE. To the extent permitted by applicable law, Recipient agrees to indemnify, defend and hold harmless Provider, UW, and the SOFTWARE authors against all claims and expenses, including legal expenses and reasonable attorneys fees, arising from Recipient’s use of the SOFTWARE.
%7.	No use of names/trademarks.  Recipient shall not use Provider’s name, or the name of any author of the SOFTWARE or that of UW, in any manner without the prior written approval of the entity or person whose name is being used.
%8.	Termination.  Without prejudice to any other rights, Provider may terminate this Agreement if Recipient fails to comply with the terms of this Agreement for any reason. Upon termination for any reason, Recipient must immediately destroy all copies of the SOFTWARE in Recipient’s possession, custody, or control.	


function [LUTs] = build_LUTs_vdphi_T1B1(dphi, alpha, TRin, varargin)
% dphi: RF phase increment in radian
% alpha: Flip angle in radian
% TR: Repetition Time in millisecond
% varargin (optional): if phantom or not 


% Range of T2 for the lookup table
R2s = 10:600;
T2s = 1./R2s;
betas = (1.0);

if length(varargin) > 0
    isPhantom = varargin{1};
else
    isPhantom = false;
end
if isPhantom
    T1s = (25:25:1000)*1e-3;
    T1s_int = 10:25:1000;
    R2s_re = 1:600;
    T1s_re = (10:1000);
else
    T1s = (100:100:1200)*1e-3;
    T1s_int = 100:100:1200;
    R2s_re = 1:300;
    T1s_re = (100:1200);
end

if length(varargin) > 1
    coeff = varargin{2};
    R1s = coeff(1)*R2s + coeff(2);
    T1s = [1000]*1e-3;
else
    T1s = [1000]*1e-3;
end

TR = TRin*1e-3;

LUTs_eta = zeros(length(betas), length(T1s), length(T2s));
LUTs_theta = zeros(length(betas), length(T1s), length(T2s));
LUTs_index = zeros(length(betas), length(T1s), length(T2s),2);



C = dphi/2;
beta_ind = 1;
for beta = betas
    T1_ind = 1;
    for T1 = T1s
        T2_ind = 1;
        for T2 = T2s
            
            if length(varargin) > 1
                %R1 corrected
                [y1] = analytical_SPGR(TR, 1/R1s(T2_ind), T2, beta*alpha, C);
                [y2] = analytical_SPGR(TR, 1/R1s(T2_ind), T2, beta*alpha, -1*C);
                LUTs_eta(beta_ind, T1_ind, T2_ind) = T2;
                LUTs_theta(beta_ind, T1_ind, T2_ind) = abs(angle(y1*conj(y2)))/2;
                LUTs_index(beta_ind, T1_ind, T2_ind,1) = beta;
                LUTs_index(beta_ind, T1_ind, T2_ind,2) = T2;
                LUTs_index(beta_ind, T1_ind, T2_ind,3) = abs(angle(y1*conj(y2)))/2;
            else
                %W/O R1 correction
                [y1] = analytical_SPGR(TR, T1, T2, beta*alpha, C);
                [y2] = analytical_SPGR(TR, T1, T2, beta*alpha, -1*C);
                LUTs_eta(beta_ind, T1_ind, T2_ind) = T2;
                LUTs_theta(beta_ind, T1_ind, T2_ind) = abs(angle(y1*conj(y2)))/2;
                LUTs_index(beta_ind, T1_ind, T2_ind,1) = beta;
                LUTs_index(beta_ind, T1_ind, T2_ind,2) = T2;
                LUTs_index(beta_ind, T1_ind, T2_ind,3) = abs(angle(y1*conj(y2)))/2;
            end

            
            T2_ind = T2_ind + 1;
        end
        
        
        T1_ind = T1_ind+1;
    end

    beta_ind = beta_ind + 1;
end


if length(T1s) == 1

    LUTs_eta = reshape(LUTs_eta,length(betas)*length(T1s)*length(T2s),1);
    LUTs_theta = reshape(LUTs_theta,length(betas)*length(T1s)*length(T2s),1);
    LUTs_index = reshape(LUTs_index,length(betas)*length(T1s)*length(T2s),3);

    LUTs_eta_re = [];
    pred_eta = [];
    mdl = [];
    R2s_re = [];
    T1s_re = [];
    
else
    
    [Y,X] = meshgrid(R2s,T1s);
    [Yq,Xq] = meshgrid(R2s_re,T1s_re);


    LUTs_eta_re = interp2(Y,X,squeeze(LUTs_eta),Yq,Xq);

    LUTs_eta = reshape(LUTs_eta,length(betas)*length(T1s)*length(T2s),1);
    LUTs_theta = reshape(LUTs_theta,length(betas)*length(T1s)*length(T2s),1);
    LUTs_index = reshape(LUTs_index,length(betas)*length(T1s)*length(T2s),3);
    
    LUTs.fitmodel_name = 'Linear function T2 = (beta, T1, phase)';
    mdl = fitlm(LUTs_index,LUTs_eta,'poly346');
    
    pred_eta = predict(mdl,LUTs_index);
    pred_eta= reshape(pred_eta, length(betas),length(T1s),length(T2s));
    LUTs_eta = reshape(LUTs_eta, length(betas),length(T1s),length(T2s));
end

LUTs.mdl = mdl;
LUTs.index = LUTs_index;
LUTs.eta = LUTs_eta;
LUTs.theta = LUTs_theta;
LUTs.T2s = T2s;
LUTs.T1s = T1s;
LUTs.betas = betas;
LUTs.pred_eta = pred_eta;
LUTs.era_re = LUTs_eta_re;
LUTs.T1s_re = T1s_re;
LUTs.R2s_re = R2s_re;
