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


addpath('../tools');
%Stanisz, G. J., Odrobina, E. E., Pun, J., Escaravage, M., Graham, S. J., Bronskill, M. J., & Henkelman, R. M. (2005). T1, T2 relaxation and magnetization transfer in tissue at 3T. Magnetic Resonance in Medicine: An Official Journal of the International Society for Magnetic Resonance in Medicine, 54(3), 507-512.
%T1/T2 of the liver at3T: 812/42ms
%T1 of the liver at 1.5T: 576/54ms

%Hamilton, G., Middleton, M. S., Hooker, J. C., Haufe, W. M., Forbang, N. I., Allison, M. A., ... & Sirlin, C. B. (2015). In vivo breath‐hold 1H MRS simultaneous estimation of liver proton density fat fraction, and T1 and T2 of water and fat, with a multi‐TR, multi‐TE sequence. Journal of Magnetic Resonance Imaging, 42(6), 1538-1543.
%T1/T2 of the liver fat at3T: 312/53.4ms
%Hamilton, G., Middleton, M. S., Bydder, M., Yokoo, T., Schwimmer, J. B., Kono, Y., ... & Sirlin, C. B. (2009). Effect of PRESS and STEAM sequences on magnetic resonance spectroscopic liver fat quantification. Journal of Magnetic Resonance Imaging: An Official Journal of the International Society for Magnetic Resonance in Medicine, 30(1), 145-152.
%Drew Z, Haouimi A, Murphy A, et al. T1 values (1.5 T). Reference article, Radiopaedia.org (Accessed on 24 Jan 2024) https://doi.org/10.53347/rID-60740
%T1 of the liver fat at 1.5T: 260/64ms
TE = 0.9e-3;
%fat_ppm = -3.4;
field_strength = 1.5;
%freq = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60]*42.58*1.5*2*pi* TE;
relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];
fat_ppm = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
fat_Hz = fat_ppm*42.58*field_strength;
fat_phase = 2*pi* fat_Hz * TE;
%fat_Hz = fat_ppm*42.58*field_strength;
%fat_phase = 2*pi* fat_Hz * TE;
TR_1p5T = 3.2e-3;
T1w_1p5T = 576e-3;
T2w_1p5T = 54e-3;
T1f_1p5T = 260e-3;
T2f_1p5T = 64e-3;

dphi = 1.5*(pi/180.0);
alpha = 18.0*(pi/180.0);
C = dphi/2.0;

[LUTs] = build_LUTs_vdphi(dphi, alpha, TR_1p5T*1000, 10.0, T1w_1p5T*1000);


%T2s = (1:0.1:50)*1e-3;
R2s = (20:1:210);
FFs = [0.0:0.005:0.32];
phase_ftg = zeros(length(FFs),length(R2s));
T2_gt = zeros(length(FFs),length(R2s));
for jj = 1:length(FFs)
    FF = FFs(jj);
    for ii=1:length(R2s)
        T2wf = 1/R2s(ii);
        R1 = (R2s(ii) + 189)/154;
        T1w_1p5T = 1./R1;
                
        
        [f_wat, f_1, epsilon_eta, beta] = analytical_SPGR_W_diffusion(TR_1p5T, T1w_1p5T, T2wf, alpha, C, 2000e-9, 100.0, TR_1p5T);
        [f_fat, f_1, epsilon_eta, beta] = analytical_SPGR_W_diffusion(TR_1p5T, T1f_1p5T, T2f_1p5T, alpha, C, 2000e-9, 100.0, TR_1p5T);
        [f_wat2, f_1, epsilon_eta, beta] = analytical_SPGR_W_diffusion(TR_1p5T, T1w_1p5T, T2wf, alpha, -C, 2000e-9, 100.0, TR_1p5T);
        [f_fat2, f_1, epsilon_eta, beta] = analytical_SPGR_W_diffusion(TR_1p5T, T1f_1p5T, T2f_1p5T, alpha, -C, 2000e-9, 100.0, TR_1p5T);
    
        f_fat = sum(f_fat.*exp(complex(0,fat_phase)).*relAmps);
        f_fat2 = sum(f_fat2.*exp(complex(0,fat_phase)).*relAmps);
        f_sum = (1.0-FF)*f_wat+FF*f_fat; f_sum2 = (1.0-FF)*f_wat2+FF*f_fat2;

        phase_ftg(jj,ii) = angle(f_sum.*conj(f_sum2))/2.0;
        T2_gt(jj,ii) = T2wf;
        
        %f_sum = (1.0-FF)*f_wat+FF.*f_fat.*sum(relAmps.*exp(complex(0,freq)));
        %f_sum2 = (1.0-FF)*f_wat2+FF*f_fat2*sum(relAmps.*exp(complex(0,freq)));
        %phase_ftg(jj,ii) = angle(f_sum.*conj(f_sum2))/2.0;
        %T2_gt(jj,ii) = T2wf;
    end
end


T2map = T2map_phase(phase_ftg, LUTs);
R2map = 1./T2map;

%%

R2_var = abs(R2map-1./T2_gt) ./ (1./T2_gt);
R2_var = convn(R2_var, ones(5,5)/25,'same');


figure;
imagesc(R2s, FFs*100, 100.0*R2_var, [0 50]);title('Bias', 'FontWeight', 'normal');xlabel('R2 (s^{-1})');ylabel('PDFF (%)');
ax = gca;ax.YDir = 'normal';
set(gca, 'fontname', 'Arial', 'FontSize',16,'FontWeight','normal','LineWidth',1);

c3=contourc(R2s(1:(end-5)), FFs(1:(end-5))*100, 100.0*R2_var(1:(end-5),1:(end-5)), [-10 10]); %indt = (c3(2,:) < R2s_crsq) & (c3(1,:) < 1) & (c3(2,:) < R2s_max); c3 = c3(:,indt);
c0=contourc(R2s(1:(end-5)), FFs(1:(end-5))*100, 100.0*R2_var(1:(end-5),1:(end-5)), [-20 20]); %indt = (c0(2,:) > R2s_crsq) & (c0(1,:) < 1) & (c0(2,:) < R2s_max) ;c0 = c0(:,indt);

hold on;
[~,indt] = sort(c3(2,:)); c3_re = c3(:,indt);
[t,indt] = sort(c0(2,:)); c0_re = c0(:,indt);
c_re = [c3_re];
plot(c0_re(1,1:end), c0_re(2,1:end), '--w','LineWidth',2); 
ylim([0 30])
xlim([20 200])


%%
%3T

TE = 0.86e-3;
field_strength = 3.0;
fat_ppm = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
fat_Hz = fat_ppm*42.58*field_strength;
fat_phase = 2*pi* fat_Hz * TE;
TR_3T = 3.1e-3;
T1w_3T = 812e-3;
T2w_3T = 54e-3;
T1f_3T = 312e-3;
T2f_3T = 53.4e-3;
dphi = 1.5*(pi/180.0);
alpha = 18.0*(pi/180.0);
C = dphi/2.0;
relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];

[LUTs] = build_LUTs_vdphi(dphi, alpha, TR_3T*1000, 10.0, T1w_3T*1000);

R2s = (20:1:210);
FFs = [0.0:0.005:0.31];
phase_ftg = zeros(length(FFs),length(R2s));
B1_ftg = ones(length(FFs),length(R2s));
T2_gt = zeros(length(FFs),length(R2s));
for jj = 1:length(FFs)
    FF = FFs(jj);
    for ii=1:length(R2s)
        T2wf = 1/R2s(ii);
        R1_curr = 0.0038.*R2s(ii) + 0.89;
        T1w_3T = 1./R1_curr;
        [f_wat] = analytical_SPGR(TR_3T, T1w_3T, T2wf, alpha, C);
        [f_fat] = analytical_SPGR(TR_3T, T1f_3T, T2f_3T, alpha, C);
        [f_wat2] = analytical_SPGR(TR_3T, T1w_3T, T2wf, alpha, -C);
        [f_fat2] = analytical_SPGR(TR_3T, T1f_3T, T2f_3T, alpha, -C);
    
        f_fat = sum(f_fat.*exp(complex(0,fat_phase)).*relAmps);
        f_fat2 = sum(f_fat2.*exp(complex(0,fat_phase)).*relAmps);
        f_sum = (1.0-FF)*f_wat+FF*f_fat; f_sum2 = (1.0-FF)*f_wat2+FF*f_fat2;
        phase_ftg(jj,ii) = angle(f_sum.*conj(f_sum2))/2.0;
        T2_gt(jj,ii) = T2wf;
    end
end


T2map = T2map_phase(phase_ftg, LUTs);
R2map = 1./T2map;

R2ss = (20:1:200);
FFss = [0.0:0.001:0.3];
[X,Y] = meshgrid(R2s,FFs);
[Xq,Yq] = meshgrid(R2ss,FFss);
T2_gt_ss = interp2(X, Y, T2_gt, Xq, Yq);
R2map_ss = interp2(X, Y, R2map, Xq, Yq);

%%

R2_var = abs(R2map-1./T2_gt) ./ (1./T2_gt);
R2_var = convn(R2_var, ones(5,5)/25,'same');

figure;
imagesc(R2s, FFs*100, (100.0*R2_var), [0 50]);title('Bias', 'FontWeight', 'normal');xlabel('R2 (s^{-1})');ylabel('PDFF (%)');
ax = gca;ax.YDir = 'normal';
set(gca, 'fontname', 'Arial', 'FontSize',16,'FontWeight','normal','LineWidth',1);

c3=contourc(R2s(1:(end-5)), FFs(1:(end-5))*100, 100.0*R2_var((1:(end-5)),(1:(end-5))), [-20 20]); %indt = (c3(2,:) < R2s_crsq) & (c3(1,:) < 1) & (c3(2,:) < R2s_max); c3 = c3(:,indt);
hold on;
[~,indt] = sort(c3(2,:)); c3_re = c3(:,indt);
[t,indt] = sort(c0(2,:)); c0_re = c0(:,indt);
c_re = [c3_re];
c_re(:,305:433) = [];

[~,c_re_grad] = gradient(c_re); c_re_grad_abs = sqrt(sum(c_re_grad.^2,1));
plot(c_re(1,1:end), c_re(2,1:end), '--w','LineWidth',2); %plot(c0_re(1,1:end), c0_re(2,1:end), '--w','LineWidth',2);
ylim([0 30])
xlim([20 200])

