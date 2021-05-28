function flg = egformat(ishot,data,a)
%eg登録データフォーマット この関数の概要をここに記述
%   詳細説明をここに記述

     filename = ['dytrans_ts_9@',num2str(ishot),'.dat'];
     A1    = ['# [Parameters]'];
     A2  = ['# Name = ''dytrans_ts_9'''];
     A3    = ['# ShotNo = ',num2str(ishot)];

     formatOut = 'mm/dd/yyyy HH:MM';

     A4  = ['# Date = ''',datestr(now, formatOut),''''];
     A5  = ['#'];
     A6  = ['# DimNo = 2'];
     A7  = ['# DimName = ''Time'', ''reff'''];
     A8  = ['# DimSize = ',num2str(a(1)-1),', ',num2str(a(2))];
     A9  = ['# DimUnit = ''s'', ''m'''];
     A10 = ['#'];
     A11 = ['# ValNo = 39'];
     A12 = ['# ValName = ''reff/a99'', ''S(dV/dreff)'', ''ne'', ''Te'', ''Ti'', ''Qe_NBI'', ''Qe_ECH'', ''Qe_dnTdt'', ''Qei_ep'', ''Qe_total'', ''Qi_NBI'', ''Qi_dnTdt'', ''Qie_ep'', ''Qi_total'', ''Qe_NBI/S'', ''Qe_ECH/S'', ''Qe_dnTdt/S'', ''Qei_ep/S'', ''Qe_total/S'', ''Qi_NBI/S'', ''Qi_dnTdt/S'', ''Qie_ep/S'', ''Qi_total/S'', ''Qe_total/S/ne'', ''Qi_total/S/ne'', ''dTe/dreff'', ''dTi/dreff'', ''dTe/dreff/Te'', ''dTi/dreff/Ti'', ''chi_e'', ''chi_i'', ''chi_e_Qei_ep_eq_zero'', ''chi_i_Qie_ep_eq_zero'', ''chi_eff'', ''Wpe'', ''Wpi'', ''Wp_kinetic'', ''tauE_kinetic'', ''tauE_kinetic_Te_eq_Ti'''
];
     A13 = ['# ValUnit = '' '', ''m2'', ''e19m-3'', ''keV'', ''keV'', ''MW'', ''MW'', ''MW'', ''MW'', ''MW'', ''MW'', ''MW'', ''MW'', ''MW'', ''MWm-2'', ''MWm-2'', ''MWm-2'', ''MWm-2'', ''MWm-2'', ''MWm-2'', ''MWm-2'', ''MWm-2'', ''MWm-2'', ''keVms-1'', ''keVms-1'', ''keVm-1'', ''keVm-1'', ''m-1'', ''m-1'', ''m2s-1'', ''m2s-1'', ''m2s-1'', ''m2s-1'', ''m2s-1'', ''kJ'', ''kJ'', ''kJ'', ''ms'', ''ms'''];
     A14 = ['# [Comments]'];
     A15 = ['# Contact: motojima.gen@nifs.ac.jp or tel.2142.'];
     A16 = ['# Ver2.2']; %versionを書く
     A17 = ['# [data]'];

     Aset = [A1 newline A2 newline A3 newline A4 newline A5 newline A6 newline A7 newline A8 newline A9 newline A10 newline A11 newline A12 newline A13 newline A14 newline A15 newline A16 newline A17];
            
     try
            flg=1;
            %cd '/Users/gmotojima/Documents/Research/MATLAB/flxgradanalysis'
            fid = fopen(filename,'wt'); %
            fprintf(fid,'%s\n',Aset);
            fprintf(fid,'%.5f,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n',data(:,:)');
            fclose(fid); % 
     catch
            flg=0;
     end
end

