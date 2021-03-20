function [] = dytrans_ts_9(shotn)

clearvars -except shotn

% shotn = 169623; %Fedelico
disp('ver1.3');
disp(shotn);
%=======================================
%fit3d_sd_autoanaを読み込む
    signame1 = 'fit3d_sd_autoana';
    filename1 = strcat('ncm/',signame1,'@',num2str(shotn),'.dat');
%    url      = strcat('http://egftp1.lhd.nifs.ac.jp/data/',signame1,'/',num2str(round(shotn-500,-3)),'/',num2str(shotn),'/000001/',signame1,'@',num2str(shotn),'.dat.zip');
        try
            %ncmFiles = unzip(url,'ncm');           
            flg_fit3d   = igetfile(shotn,'fit3d_sd_autoana'); %for local PC
% %             command     = ['igetfile -s ',num2str(shotn),' -d ',signame1,' -o ncm/',signame1,'@',num2str(shotn),'.dat']; %for egcalc
% %             system(command); %for egcalc
% %             flg_fit3d = 1;   %for egcalc
            
            header      = importdata(filename1);

            tmp         = extractAfter(header.textdata(8),'=');
            timen_auto  = str2num(tmp{1}); 
            rown_auto   = header.data;
            rhon_auto   = rown_auto(2);
            nbn_auto    = rown_auto(1);
            fit3d  = readmatrix(filename1);         %for local PC
% %             fit3d  = csvread(filename1, 19, 0); %for egcalc
            fit3dn = size(fit3d);               
            
             delete(filename1);
        catch
           warning(strcat('something wrong! no file? ',signame1,'@',num2str(shotn)));
            flg_fit3d = 0;
            return;
        end   

%=======================================

data_fit3d = zeros(timen_auto,rhon_auto,fit3dn(2),nbn_auto); %Time, rho, val, nbiNo
intsn = 0;
intfn = 0;
time_auto = zeros(1,timen_auto);
        for j = 1:timen_auto
            time_auto(j) = fit3d((j-1).*rhon_auto.*nbn_auto+1,1);
            for i = 1:nbn_auto
                
                intsn = intfn + 1;
                intfn = intfn + rhon_auto;
                
                data_fit3d(j,1:rhon_auto,:,i) = fit3d(intsn:intfn,:);

            end
        end
        
%autoanaの読み込み終了

%=======================================
%tsmap_smoothを読み込む
    signame1 = 'tsmap_smooth';
    filename1 = strcat('ncm/',signame1,'@',num2str(shotn),'.dat');
%    url      = strcat('http://egftp1.lhd.nifs.ac.jp/data/',signame1,'/',num2str(round(shotn-500,-3)),'/',num2str(shotn),'/000001/',signame1,'@',num2str(shotn),'.dat.zip');
        try
            %ncmFiles = unzip(url,'ncm');
            %filename1 = 'ncm/tsmap_smooth@168641.dat';
            flg_tsmap = igetfile(shotn,'tsmap_smooth');%for local PC
% %             command     = ['igetfile -s ',num2str(shotn),' -d ',signame1,' -o ncm/',signame1,'@',num2str(shotn),'.dat']; %for egcalc
% %             system(command); %for egcalc
% %             flg_tsmap = 1;   %for egcalc
            
            header     = importdata(filename1);
            
            tmp    = extractAfter(header.textdata(8),'=');
            timen_tsmap   = str2num(tmp{1}); 
            rown_tsmap   = header.data;
            rhon_tsmap   = rown_tsmap(1);
            tsmap_s      = readmatrix(filename1);         %for local PC
% %             tsmap_s      = csvread(filename1, 28, 0); %for egcalc
            tsmap_sn     = size(tsmap_s);             
            
             delete(filename1);

            %=======================================
            data_tsmap = zeros(timen_tsmap,rhon_tsmap,tsmap_sn(2)); %Time, rho, val, nbiNo
            intsn = 0;
            intfn = 0;
            time_tsmap_s = zeros(1,timen_tsmap);
                    for i = 1:timen_tsmap
                            time_tsmap(i) = tsmap_s((i-1).*rhon_tsmap+1,1);

                            intsn = intfn + 1;
                            intfn = intfn + rhon_tsmap;

                            data_tsmap(i,:,:) = tsmap_s(intsn:intfn,:);

                    end
        catch
            warning(strcat('something wrong! no file? ',signame1,'@',num2str(shotn)));
            flg_tsmap_s = 0;
            return;
        end   
%tsmapの読み込み終了
%=======================================
%=======================================
%lhdgaussを読み込む
    signame1 = 'LHDGAUSS_DEPROF';
    filename1 = strcat('ncm/',signame1,'@',num2str(shotn),'.dat');
%    url      = strcat('http://egftp1.lhd.nifs.ac.jp/data/',signame1,'/',num2str(round(shotn-500,-3)),'/',num2str(shotn),'/000001/',lower(signame1),'@',num2str(shotn),'.dat.zip');
        try
            %ncmFiles = unzip(url,'ncm');
            %filename1 = 'ncm/lhdgauss_deprof@168641.dat';
            flg_gauss = igetfile(shotn,'LHDGAUSS_DEPROF'); %for local PC
% %             command     = ['igetfile -s ',num2str(shotn),' -d ',signame1,' -o ncm/',signame1,'@',num2str(shotn),'.dat']; %for egcalc
% %             system(command); %for egcalc
% %             flg_gauss = 1;   %for egcalc
            
            header = importdata(filename1);

            tmp    = extractAfter(header.textdata(8),'=');
            timen_gauss   = str2num(tmp{1}); 
            rown_gauss   = header.data;
            rhon_gauss   = rown_gauss(1);
            gauss  = readmatrix(filename1,'NumHeaderLines',18);  %for local PC
% %             gauss  = csvread(filename1, 18, 0);     %for egcalc
            gauss_n = size(gauss);                  
            
             delete(filename1);

        data_gauss = zeros(timen_gauss,rhon_gauss,gauss_n(2)); %Time, rho, val, nbiNo
        intsn = 0;
        intfn = 0;
        time_gauss = zeros(1,timen_gauss);
                for i = 1:timen_gauss
                        time_gauss(i) = gauss((i-1).*rhon_gauss+1,1);

                        intsn = intfn + 1;
                        intfn = intfn + rhon_gauss;

                        data_gauss(i,:,:) = gauss(intsn:intfn,:);

                end
        catch
            warning(strcat('something wrong! no file? ',signame1,'@',num2str(shotn)));
            flg_gauss = 0;
            return;
        end   

%lhdgaussを読み込終了
%=======================================

%=======================================
%cx9を読み込む
    signame1 = 'cxsmap9_poly6';
    filename1 = strcat('ncm/',signame1,'@',num2str(shotn),'.dat');
%    url      = strcat('http://egftp1.lhd.nifs.ac.jp/data/',signame1,'/',num2str(round(shotn-500,-3)),'/',num2str(shotn),'/000001/',signame1,'@',num2str(shotn),'.dat.zip');
        try
%            ncmFiles = unzip(url,'ncm');
            flg_cx9 = igetfile(shotn,'cxsmap9_poly6'); %for local PC
% %             command  = ['igetfile -s ',num2str(shotn),' -d ',signame1,' -o ncm/',signame1,'@',num2str(shotn),'.dat']; %for egcalc
% %             system(command); %for egcalc
% %             flg_cx9 = 1;     %for egcalc
            
            %filename1 = 'ncm/lhdcx9_deprof@168641.dat';
            header = importdata(filename1);

            tmp    = extractAfter(header.textdata(8),'=');
            timen_cx9    = str2num(tmp{1}); 
            rown_cx9     = header.data;
            rhon_cx9_0     = rown_cx9(1);
            cx9_0        = readmatrix(filename1,'NumHeaderLines',21); %ヘッダー21行 %for local PC
% %             cx9_0        = csvread(filename1,21, 0); %for egcalc
            cx9_n_0      = size(cx9_0);
            
            %=======================================
            %array 5,7抽出
            cx9_sele   = (cx9_0(:,3) == 5 | cx9_0(:,3) == 7);
            cx9        = cx9_0(cx9_sele,:);
            cx9_n      = size(cx9);
            rhon_cx9   = cx9_n(1)./timen_cx9;
            %=======================================
             delete(filename1);


            data_cx9 = zeros(timen_cx9,rhon_cx9,cx9_n(2)); %Time, rho, val
            intsn = 0;
            intfn = 0;
            time_cx9 = zeros(1,timen_cx9);
                    for j = 1:timen_cx9
                            time_cx9(j) = cx9((j-1).*rhon_cx9+1,1);

                            intsn = intfn + 1;
                            intfn = intfn + rhon_cx9;

                            data_cx9(j,:,:) = cx9(intsn:intfn,:);

                    end
            %=======================================
            %dniTi/dtを求めるために有効時間スライスを抽出
            cx9_ex = cx9;
            cx9_ex(sum(cx9(:,4:10),2) == 0,:) = []; %必要なTiデータだけを抽出
            cx9_n_ex = size(cx9_ex);
            timen_cx9_ex = cx9_n_ex(1)./rhon_cx9;

            data_cx9_ex = zeros(timen_cx9_ex,rhon_cx9,cx9_n(2)); %Time, rho, val
            intsn = 0;
            intfn = 0;
            time_cx9_ex = zeros(1,timen_cx9_ex);
                    for j = 1:timen_cx9_ex
                            time_cx9_ex(j) = cx9_ex((j-1).*rhon_cx9+1,1);

                            intsn = intfn + 1;
                            intfn = intfn + rhon_cx9;

                            data_cx9_ex(j,:,:) = cx9_ex(intsn:intfn,:);

                    end
                    data_cx9_ex0 = data_cx9_ex; %map-Ti前の有効Tiデータを抽出
            %=======================================
            
            %=======================================
            %マッピングと書き換え
            data_cx9_ex_map      = mapTi(data_cx9_ex, time_cx9_ex, time_tsmap, timen_tsmap); %Tiをtsmapにマッピング
      %data_cx9_ex_map_sele = mapTi_select(data_cx9_ex_map,timen_cx9_ex,time_cx9_ex); %mapTi_selectは使用しない。
            data_cx9        = data_cx9_ex_map; %data_cx9をマッピングに書き換える。
            data_cx9_ex     = data_cx9_ex_map; %data_cx9_exをマッピングに書き換える。data_cx9とdata_cx9_exを区別なくす
            %マッピングと書き換え終了
            %=======================================
            
        catch
            warning(strcat('something wrong! no file? ',signame1,'@',num2str(shotn)));
            flg_cx9 = 0;
            return;
        end   
        
%cx9を読み込終了

%tを選択
t1 =3.7003;

aaQ_out = zeros(timen_tsmap,rhon_tsmap,41); 
%'reff/a99', 'ne', 'Te', 'Ti', 'Qe_NBI', 'Qe_ECH', 'Qe_dnTdt', 'Qei_ep', 'Qe_total', 'Qe_NBI/S', 'Qe_ECH/S', 'Qe_dnTdt/S', 'Qei_ep/S', 'Qe_total/S', 'Qi_NBI', 'Qi_dnTdt', 'Qie_ep', 'Qi_total', 'Qi_NBI/S', 'Qi_dnTdt/S', 'Qie_ep/S', 'Qi_total/S', 'dVdreff', 'dTe/dreff', 'dTe/dreff/Te', 'dTi/dreff', 'dTi/dreff/Ti', 'chi_e', 'chi_e_noep', 'chi_i', 'chi_i_noep', 'chi_i_noep_Te=Ti', 'chi_eff', 'chi_eff_noep_Te=Ti'  

for i=1:timen_tsmap-1 %-1にしているのは、Qedeltaを求める際、次のフレームのtsmapを評価する必要があるため。
t1 = time_tsmap(1,i);
%=======================================
%NBIとneTe
        [~,idx] = min(abs(data_tsmap(:,1,1)-t1));
        d_tsmap = data_tsmap(idx,:,:);
        d_tsmap = reshape(d_tsmap,rhon_tsmap,tsmap_sn(2));
        %plot(d_tsmap(:,5),d_tsmap(:,4));

        rho_tsmap     = d_tsmap(:,5);    %reff/a99:5, reff/a999:6
        reff_tsmap    = d_tsmap(:,2);
        Te_tsmap      = d_tsmap(:,3);
        ne_tsmap      = d_tsmap(:,4);
        dVdreff_tsmap = d_tsmap(:,8);

         %grad_Te     = gradient(Te_tsmap);
         %grad_reff   = gradient(reff_tsmap);        
         %dTedreff    = grad_Te./grad_reff;   %重要
         [flg_gradTe,dTedreff] = gradT(rhon_tsmap,reff_tsmap,Te_tsmap,3);%重要 数値は前後何点取るか
         
         L_Te        = dTedreff./Te_tsmap;

                [val_nb,t1_nb] = min(abs(data_fit3d(:,1,1,1)-t1));
    
    if val_nb > 1/30
        rho_flxgrad(i,:)          = rho_tsmap;
        dTedreff_flxgrad(i,:)     = -dTedreff;      
        Qe_flxne_flxgrad(i,:)     = 0;
        
        Qe_NBI        = zeros(rhon_tsmap,1);
        intQe_Pe      = zeros(rhon_tsmap,1);
        intQe_delta   = zeros(rhon_tsmap,1);
        Qi_NBI        = zeros(rhon_tsmap,1);
        intQi_Pi      = zeros(rhon_tsmap,1);
        intQi_delta   = zeros(rhon_tsmap,1);

        Qe_flx_set    = zeros(rhon_tsmap,1);
        chi_i_TeTi    = zeros(rhon_tsmap,1);
        Qi_calc_flxni = zeros(rhon_tsmap,1);
        intvol_neTe   = zeros(rhon_tsmap,1);
        
    else
    
        %t1_nb = find(data_fit3d(:,1,1,1)==t1);

        d_fit3_nb1 = data_fit3d(t1_nb,:,:,1);
        d_fit3_nb2 = data_fit3d(t1_nb,:,:,2);
        d_fit3_nb3 = data_fit3d(t1_nb,:,:,3);
        d_fit3_nb4 = data_fit3d(t1_nb,:,:,4);
        d_fit3_nb5 = data_fit3d(t1_nb,:,:,5);

        d_fit3_nb1 = reshape(d_fit3_nb1,rhon_auto,fit3dn(2));
        d_fit3_nb2 = reshape(d_fit3_nb2,rhon_auto,fit3dn(2));
        d_fit3_nb3 = reshape(d_fit3_nb3,rhon_auto,fit3dn(2));
        d_fit3_nb4 = reshape(d_fit3_nb4,rhon_auto,fit3dn(2));
        d_fit3_nb5 = reshape(d_fit3_nb5,rhon_auto,fit3dn(2));

        %===========================
        %tsmapのマッピングに合わせる場合
             
         rho_auto      = d_fit3_nb1(:,3); %rho_autoana

         d_fit3_nb1_sp = interp1(rho_auto,d_fit3_nb1,rho_tsmap);
         d_fit3_nb2_sp = interp1(rho_auto,d_fit3_nb2,rho_tsmap);
         d_fit3_nb3_sp = interp1(rho_auto,d_fit3_nb3,rho_tsmap);
         d_fit3_nb4_sp = interp1(rho_auto,d_fit3_nb4,rho_tsmap);
         d_fit3_nb5_sp = interp1(rho_auto,d_fit3_nb5,rho_tsmap);

         
         %電子系
         Qe1 = d_fit3_nb1_sp(:,9);
         Qe2 = d_fit3_nb2_sp(:,9);
         Qe3 = d_fit3_nb3_sp(:,9);
         Qe4 = d_fit3_nb4_sp(:,9);
         Qe5 = d_fit3_nb5_sp(:,9);

         Qe_NBI = (Qe1+Qe2+Qe3+Qe4+Qe5).*1e-3.*6.242e18; %6.242e+18 eV/J

         Qe_flx1 = Qe1./dVdreff_tsmap;
         Qe_flx2 = Qe2./dVdreff_tsmap;
         Qe_flx3 = Qe3./dVdreff_tsmap;
         Qe_flx4 = Qe4./dVdreff_tsmap;
         Qe_flx5 = Qe5./dVdreff_tsmap;

         Qe_flx_set = (Qe_flx1+Qe_flx2+Qe_flx3+Qe_flx4+Qe_flx5).*1e-3.*6.242e18;
         Qe_flxne = Qe_flx_set./ne_tsmap./1e19; %重要２
        %===========================

        %======================================
        %次のフレームのtsmap
        d_tsmap_nxt = data_tsmap(idx+1,:,:);  
        d_tsmap_nxt = reshape(d_tsmap_nxt,rhon_tsmap,tsmap_sn(2)); 

        reff_tsmap_nxt    = d_tsmap_nxt(:,2);
        Te_tsmap_nxt      = d_tsmap_nxt(:,3);
        ne_tsmap_nxt      = d_tsmap_nxt(:,4);
        dVdreff_tsmap_nxt = d_tsmap_nxt(:,8);

        Te_tsmap_nxt_sp           = interp1(reff_tsmap_nxt,Te_tsmap_nxt,reff_tsmap);
        ne_tsmap_nxt_sp           = interp1(reff_tsmap_nxt,ne_tsmap_nxt,reff_tsmap);
        dVdreff_tsmap_nxt_sp      = interp1(reff_tsmap_nxt,dVdreff_tsmap_nxt,reff_tsmap);
        
        neTe = ne_tsmap.*Te_tsmap;
        neTe_nxt = ne_tsmap_nxt_sp.*Te_tsmap_nxt_sp;

        dneTedt       = (neTe_nxt-neTe)./(d_tsmap_nxt(1,1)-d_tsmap(1,1));
        
        dvol_neTe = 3./2.*neTe.*dVdreff_tsmap.*gradient(reff_tsmap).*1e19./6.242e+18; %dvolume keV -> kJ %6.242e+18 eV/J
        intvol_neTe = cumtrapz(dvol_neTe); %volume integral kJ
        %======================================

        %Peのみの寄与
        Pe1  = d_fit3_nb1_sp(:,6); %MW m-3
        Pe2  = d_fit3_nb2_sp(:,6); %MW m-3
        Pe3  = d_fit3_nb3_sp(:,6); %MW m-3
        Pe4  = d_fit3_nb4_sp(:,6); %MW m-3
        Pe5  = d_fit3_nb5_sp(:,6); %MW m-3

        totalPe  = (Pe1+Pe2+Pe3+Pe4+Pe5).*1e6.*1e-3.*6.242e18; %keV*m^-3s^-1 MW->W, eV->keV
        dQe_Pe   = totalPe.*dVdreff_tsmap.*gradient(reff_tsmap); %keVs^-1

        dQe_Pe1                 = dQe_Pe;
        dQe_Pe1(isnan(dQe_Pe1)) = 0;
        intQe_Pe                = cumtrapz(dQe_Pe1); %keVs^-1 出力1

        %neTeのみの寄与
        Qe_delta    =  3./2.*dneTedt.*1e19; %keV*m^-3s^-1
        dQe_delta   = Qe_delta.*dVdreff_tsmap.*gradient(reff_tsmap); %keVs^-1

        dQe_delta1 = dQe_delta;
        dQe_delta1(isnan(dQe_delta1)) = 0;
        intQe_delta = cumtrapz(dQe_delta1); %keVs^-1 出力2   
    
%NBIとneTe終了

%=======================================   
%Qi, Pi

        %===========================        
        %イオン系
         Qi1 = d_fit3_nb1_sp(:,10);
         Qi2 = d_fit3_nb2_sp(:,10);
         Qi3 = d_fit3_nb3_sp(:,10);
         Qi4 = d_fit3_nb4_sp(:,10);
         Qi5 = d_fit3_nb5_sp(:,10);

         Qi_NBI = (Qi1+Qi2+Qi3+Qi4+Qi5).*1e-3.*6.242e18; %6.242e+18 eV/J

         Qi_flx1 = Qi1./dVdreff_tsmap;
         Qi_flx2 = Qi2./dVdreff_tsmap;
         Qi_flx3 = Qi3./dVdreff_tsmap;
         Qi_flx4 = Qi4./dVdreff_tsmap;
         Qi_flx5 = Qi5./dVdreff_tsmap;

         Qi_flx_set = (Qi_flx1+Qi_flx2+Qi_flx3+Qi_flx4+Qi_flx5).*1e-3.*6.242e18;
         Qi_flxne = Qi_flx_set./ne_tsmap./1e19; %重要２
        %===========================

        %chi_i 計算 Te=Tiと仮定の場合
        Qi_calc       = Qi_NBI-intQe_delta;
        Qi_calc_flx   = Qi_calc./dVdreff_tsmap;
        Qi_calc_flxni = Qi_calc_flx./ne_tsmap./1e19;
        chi_i_TeTi    = Qi_calc_flxni./(-dTedreff);
        %chi_i 計算 Te=Tiと仮定の場合終                                  
        
        
        
        %Piのみの寄与
        Pi1  = d_fit3_nb1_sp(:,7); %MW m-3
        Pi2  = d_fit3_nb2_sp(:,7); %MW m-3
        Pi3  = d_fit3_nb3_sp(:,7); %MW m-3
        Pi4  = d_fit3_nb4_sp(:,7); %MW m-3
        Pi5  = d_fit3_nb5_sp(:,7); %MW m-3

        totalPi  = (Pi1+Pi2+Pi3+Pi4+Pi5).*1e6.*1e-3.*6.242e18; %keV*m^-3s^-1 MW->W, eV->keV
        dQi_Pi   = totalPi.*dVdreff_tsmap.*gradient(reff_tsmap); %keVs^-1

        dQi_Pi1                 = dQi_Pi;
        dQi_Pi1(isnan(dQi_Pi1)) = 0;
        intQi_Pi                = cumtrapz(dQi_Pi1); %keVs^-1 出力1      
        %===========================
        
%Qi, Pi終了
    end

%=======================================
%ECH
if flg_gauss == 0
   
    Qe_ech_sp = zeros(rhon_tsmap,1);
    
end


if flg_gauss == 1
    [val_ech,t1_ech] = min(abs(data_gauss(:,1,1)-t1));
        if val_ech > 1/30
                           Qe_ech_sp = zeros(rhon_tsmap,1);
        else

            d_gauss = data_gauss(t1_ech,:,:);
            d_gauss = reshape(d_gauss,rhon_gauss,gauss_n(2));

            Qe_ech   = d_gauss(:,4).*1e6.*1e-3.*6.242e18; %MW -> keVs^-1 
            reff_ech = d_gauss(:,2);

            Qe_ech_sp = interp1(reff_ech,Qe_ech,reff_tsmap); %出力3
        end        
end
%ECH終了
%=======================================

        %chi_e 計算 Te=Tiと仮定の場合
        Qe_calc       = Qe_NBI+Qe_ech_sp-intQe_delta;
        Qe_calc_flx   = Qe_calc./dVdreff_tsmap;
        Qe_calc_flxne = Qe_calc_flx./ne_tsmap./1e19;
        chi_e_noep    = Qe_calc_flxne./(-dTedreff);
        %chi_e 計算 Te=Tiと仮定の場合終

        chi_eff_noep  = (Qe_calc_flxne+Qi_calc_flxni)./(-2.*dTedreff); %chi_eff noequipの場合(Te=Ti)



%Ti
if flg_cx9 == 1
    [val_cx9,t1_cx9] = min(abs(data_cx9(:,1,1)-t1));
    %[val_cx9,t1_cx9] = min(abs(data_cx9_ex(:,1,1)-t1)); %抽出データを用いる

        if val_cx9 < 1/30 %もしTi mapを使う場合は、1/100以下でも良い。なんなら一致でもいいと思う。default 1/30 @2021-03-15

            d_cx9 = data_cx9(t1_cx9,:,:);
            d_cx9 = reshape(d_cx9,rhon_cx9,cx9_n(2));

            Ti_cx9 = [d_cx9(:,11), d_cx9(:,4), d_cx9(:,22), d_cx9(:,20)]; %reff,Ti, ni, dVdreff
            Ti_cx9_sort_tmp = sortrows(Ti_cx9);
            
            [C1,idex1]   = unique(Ti_cx9_sort_tmp(:,1)); %重複したデータを探す
            Ti_cx9_sort = Ti_cx9_sort_tmp(idex1,:);      %重複したデータを消す
            
try
            Ti_cx9_sp      = interp1(Ti_cx9_sort(:,1),Ti_cx9_sort(:,2),reff_tsmap);
            ni_cx9_sp      = interp1(Ti_cx9_sort(:,1),Ti_cx9_sort(:,3),reff_tsmap);
            dVdreff_cx9_sp = interp1(Ti_cx9_sort(:,1),Ti_cx9_sort(:,4),reff_tsmap);
catch
   warning(strcat('something wrong! interporation might be missing T=',num2str(data_cx9(t1_cx9,1,1)),'s(',num2str(t1_cx9),')'));
   return;
end
            Ti_cx9_sp(isnan(Ti_cx9_sp)) = 0;       
            if sum(Ti_cx9_sp) > 0
            %eqipartionを評価    
                            ln_lambda = 17;      %クーロン対数
                            me = 9.10938356e-31; %kg
                            mi = 1.6726231e-27; %kg

                            tau_e = 1.09e16.*Te_tsmap.^(3/2)./(ne_tsmap.*1e19.*ln_lambda); %s Tokamaks 14.6
                            nu_e  = 1./tau_e;
                            nu_ep = (2.*me)./mi.*nu_e; %HLOM+, CPP 2017

                            dPeq = nu_ep.*(Te_tsmap-Ti_cx9_sp).*ne_tsmap.*1e19.*dVdreff_tsmap.*gradient(reff_tsmap); %s^-1 keV

                                dPeq1 = dPeq;
                                dPeq1(isnan(dPeq1)) = 0;
                                intPeq = cumtrapz(dPeq1); %keVs^-1
                                intPeqflx = intPeq./dVdreff_tsmap;
                                
                                  aaQ_out(i,:,11) = -intPeq.*1e-3./6.242e18;    %keVs-1 -> MW
                                  aaQ_out(i,:,20) = -intPeqflx.*1e-3./6.242e18; %keVs-1m-2 -> MW/m2
            %eqipartionを評価終了
            
                %grad_Ti         = gradient(Ti_cx9_sp);
                %dTidreff        = grad_Ti./grad_reff; %重要
                [flg_gradTi,dTidreff] = gradT(rhon_tsmap,reff_tsmap,Ti_cx9_sp,3);%重要 数値は前後何点取るか

                L_Ti            = dTidreff./Ti_cx9_sp;
            
            %次のフレームを求める。
            %======================================
            %次のフレームのcx9

            [val_cx9_ex,t1_cx9_ex] = min(abs(data_cx9_ex(:,1,1)-t1)); %抽出データを用いる

s_data_cx9_ex = size(data_cx9_ex);            
if t1_cx9_ex < s_data_cx9_ex(1)
    
                d_cx9_ex_nxt = data_cx9_ex(t1_cx9_ex+1,:,:);  
                d_cx9_ex_nxt = reshape(d_cx9_ex_nxt,rhon_cx9,cx9_n(2)); 

        %         reff_cx9_ex_nxt    = d_cx9_ex_nxt(:,11);
        %         Ti_cx9_ex_nxt      = d_cx9_ex_nxt(:,4);
        %         ni_cx9_ex_nxt      = d_cx9_ex_nxt(:,22);
        %         dVdreff_cx9_ex_nxt = d_cx9_ex_nxt(:,20);
                Ti_cx9_nxt = [d_cx9_ex_nxt(:,11), d_cx9_ex_nxt(:,4), d_cx9_ex_nxt(:,22), d_cx9_ex_nxt(:,20)]; %reff,Ti, ni, dVdreff
                Ti_cx9_nxt_sort_tmp = sortrows(Ti_cx9_nxt);
                [C2,idex2]          = unique(Ti_cx9_nxt_sort_tmp(:,1)); %重複したデータを探す
                Ti_cx9_nxt_sort     = Ti_cx9_nxt_sort_tmp(idex2,:);     %重複したデータを消す
    try
                Ti_cx9_ex_nxt_sp      = interp1(Ti_cx9_nxt_sort(:,1),Ti_cx9_nxt_sort(:,2),reff_tsmap);
                ni_cx9_ex_nxt_sp      = interp1(Ti_cx9_nxt_sort(:,1),Ti_cx9_nxt_sort(:,3),reff_tsmap);
                dVdreff_cx9_ex_nxt_sp = interp1(Ti_cx9_nxt_sort(:,1),Ti_cx9_nxt_sort(:,4),reff_tsmap);

                niTi = ni_cx9_sp.*Ti_cx9_sp;
                niTi_ex_nxt = ni_cx9_ex_nxt_sp.*Ti_cx9_ex_nxt_sp;


                dniTidt       = (niTi_ex_nxt-niTi)./(d_cx9_ex_nxt(1,1)-d_cx9(1,1));
                
                dvol_niTi = 3./2.*niTi.*dVdreff_tsmap.*gradient(reff_tsmap).*1e19./6.242e+18; %dvolume keV -> kJ %6.242e+18 eV/J
                intvol_niTi = cumtrapz(dvol_niTi); %volume integral kJ
                
                intvol_nT = intvol_neTe + intvol_niTi;
                %======================================  

                %niTiのみの寄与
                Qi_delta    =  3./2.*dniTidt.*1e19; %keV*m^-3s^-1
                dQi_delta   = Qi_delta.*dVdreff_tsmap.*gradient(reff_tsmap); %keVs^-1

                dQi_delta1 = dQi_delta;
                dQi_delta1(isnan(dQi_delta1)) = 0;
                intQi_delta = cumtrapz(dQi_delta1); %keVs^-1 出力2
  
    catch
          warning(strcat('something wrong! interporation might be missing T=',num2str(data_cx9_ex(t1_cx9_ex,1,1)),'s(',num2str(t1_cx9_ex),')'));
          return;
    end    
end
                
            else
               
                intPeq            = zeros(rhon_tsmap,1);
                Ti_cx9_sp         = zeros(rhon_tsmap,1);
                Qe_total          = zeros(rhon_tsmap,1);
                Qe_total_flx      = zeros(rhon_tsmap,1);
                intQi_delta       = zeros(rhon_tsmap,1);
                Qi_total          = zeros(rhon_tsmap,1);
                dTidreff          = zeros(rhon_tsmap,1);
                L_Ti              = zeros(rhon_tsmap,1);
                chi_e             = zeros(rhon_tsmap,1);
                chi_i             = zeros(rhon_tsmap,1);
                chi_eff           = zeros(rhon_tsmap,1);
                intvol_nT         = zeros(rhon_tsmap,1);
                intvol_niTi       = zeros(rhon_tsmap,1);
                
            end
    
            %Qeトータルの計算
            if sum(intPeq) ~= 0
                 Qe_total = Qe_NBI + Qe_ech_sp - intQe_delta - intPeq;
                 Qi_total = Qi_NBI - intQi_delta + intPeq;
                 Qi_noep  = Qi_NBI - intQi_delta;
                 Qeff_total = Qe_total + Qi_total;
                 
                 Qe_total_flx     = Qe_total./dVdreff_tsmap;
                 Qi_total_flx     = Qi_total./dVdreff_tsmap;
                 Qi_noep_flx      = Qi_noep./dVdreff_tsmap;
                 Qeff_total_flx   = Qeff_total./dVdreff_tsmap;
                 
                 Qe_total_flxne   = Qe_total_flx./ne_tsmap./1e19;
                 Qi_total_flxni   = Qi_total_flx./ni_cx9_sp./1e19;
                 Qi_noep_flxni    = Qi_noep_flx./ni_cx9_sp./1e19;
                 Qeff_total_flxne = Qeff_total_flx./ne_tsmap./1e19;
                 
                 chi_e = Qe_total_flxne ./ (-dTedreff);
                 chi_i = Qi_total_flxni ./ (-dTidreff);
                 chi_i_noep = Qi_noep_flxni ./ (-dTidreff);
                 chi_eff = Qeff_total_flxne ./ (-dTedreff-dTidreff);
            else
                 Qi_total_flx      = zeros(rhon_tsmap,1);
                 Qe_total_flxne    = zeros(rhon_tsmap,1);
                 Qi_total_flxni    = zeros(rhon_tsmap,1);
                 chi_i_noep        = zeros(rhon_tsmap,1);
            end
                                        aaQ_out(i,:,7)  = Ti_cx9_sp;                                     
                                        aaQ_out(i,:,12) = Qe_total.*1e-3./6.242e18;    %keVs-1 -> MW; 
                                        aaQ_out(i,:,14) = intQi_delta.*1e-3./6.242e18; %keVs-1 -> MW;
                                        aaQ_out(i,:,15) = intPeq.*1e-3./6.242e18;      %keVs-1 -> MW;
                                        aaQ_out(i,:,16) = Qi_total.*1e-3./6.242e18;    %keVs-1 -> MW;
                                        aaQ_out(i,:,21) = Qe_total_flx.*1e-3./6.242e18;               %keVs-1m-2 -> MW/m2;
                                        aaQ_out(i,:,23) = intQi_delta./dVdreff_tsmap.*1e-3./6.242e18; %keVs-1m-2 -> MW/m2;
                                        aaQ_out(i,:,24) = intPeq./dVdreff_tsmap.*1e-3./6.242e18;      %keVs-1m-2 -> MW/m2;
                                        aaQ_out(i,:,25) = Qi_total_flx.*1e-3./6.242e18;               %keVs-1m-2 -> MW/m2;
                                        aaQ_out(i,:,26) = Qe_total_flxne;
                                        aaQ_out(i,:,27) = Qi_total_flxni;
                                        aaQ_out(i,:,29) = dTidreff;
                                        aaQ_out(i,:,31) = L_Ti;
                                        aaQ_out(i,:,32) = chi_e;
                                        aaQ_out(i,:,33) = chi_i;
                                        aaQ_out(i,:,35) = chi_i_noep;
                                        aaQ_out(i,:,36) = chi_eff;
                                        %aaQ_out(i,:,38) = chi_eff_noep;
                                        aaQ_out(i,:,37) = intvol_neTe;
                                        aaQ_out(i,:,38) = intvol_niTi;
                                        aaQ_out(i,:,39) = intvol_nT;
                                        aaQ_out(i,:,40) = intvol_nT./(Qe_total.*1e-3./6.242e18+Qi_total.*1e-3./6.242e18); %kJ/MW = ms
                                        aaQ_out(i,:,41) = 2.*intvol_neTe./(Qe_total.*1e-3./6.242e18+Qi_total.*1e-3./6.242e18); %kJ/MW = ms
                                        
        end                            
end                              
%Ti終了
%======================================= 


                                aaQ_out(i,:,1)  = time_tsmap(1,i);
                                aaQ_out(i,:,2)  = reff_tsmap;
                                aaQ_out(i,:,3)  = rho_tsmap;
                                aaQ_out(i,:,4)  = dVdreff_tsmap;
                                aaQ_out(i,:,5)  = ne_tsmap;
                                aaQ_out(i,:,6)  = Te_tsmap;
                                aaQ_out(i,:,8)  = Qe_NBI.*1e-3./6.242e18;       %keVs-1 -> MW;       %autoanaから                                
                                aaQ_out(i,:,9)  = Qe_ech_sp.*1e-3./6.242e18;    %keVs-1 -> MW;    %lhdgaussから
                                aaQ_out(i,:,10) = intQe_delta.*1e-3./6.242e18;  %keVs-1 -> MW;  %neTe分
                                aaQ_out(i,:,13) = Qi_NBI.*1e-3./6.242e18;       %keVs-1 -> MW;                                  
                                aaQ_out(i,:,17) = Qe_flx_set.*1e-3./6.242e18;                 %keVs-1m-2 -> MW/m2;
                                aaQ_out(i,:,18) = Qe_ech_sp./dVdreff_tsmap.*1e-3./6.242e18;   %keVs-1m-2 -> MW/m2;     
                                aaQ_out(i,:,19) = intQe_delta./dVdreff_tsmap.*1e-3./6.242e18; %keVs-1m-2 -> MW/m2;                        
                                aaQ_out(i,:,22) = Qi_NBI./dVdreff_tsmap.*1e-3./6.242e18;      %keVs-1m-2 -> MW/m2;
                                aaQ_out(i,:,28) = dTedreff;
                                aaQ_out(i,:,30) = L_Te;                              
                                aaQ_out(i,:,34) = chi_e_noep;
                                %aaQ_out(i,:,36) = chi_i_TeTi;
                                

end

aaQ_out(isnan(aaQ_out)) = 0;
aaQ_out(~isfinite(aaQ_out))=0;

aaQ_out_mask = mask_data(aaQ_out,timen_cx9_ex,time_cx9_ex); %有効Tiショット以外をマスク

aaQ_out_regi = zeros(rhon_tsmap.*(timen_tsmap-1),41);
for j=1:timen_tsmap-1 %dneTe/dtを評価していたため１コマ少ない

aaQ_out_regi((j-1).*rhon_tsmap+1:j.*rhon_tsmap,:) = aaQ_out_mask(j,:,:);    


end

%equipartionの量をrho>1で0にする
aaQ_out_regi1 = aaQ_out_regi;
aaQ_out_regi1(aaQ_out_regi1(:,3) >1,8:27)  = 0;
aaQ_out_regi1(aaQ_out_regi1(:,3) >1,32:41) = 0;

size_aaQ = size(aaQ_out);
flg_save = egformat(shotn, aaQ_out_regi1,size_aaQ);

if flg_cx9 == 0
    disp('cxsmap9_poly6 missing')
end

if flg_fit3d == 0
    disp('fit3d_sd_autoana missing')
end

if flg_tsmap == 0
    disp('tsmap_smooth missing')
end

if flg_gauss == 0
    disp('LHDGAUSS_DEPROF missing')
end

if flg_tsmap == 1 & flg_cx9 == 1 & flg_tsmap == 1 & flg_gauss == 1
    message = ['dytrans_ts_9@',num2str(shotn),', iregist succeeded'];
    %flg_iregist = iregist_data('dytrans_ts_9',shotn);
else
    message = ['dytrans_ts_9@',num2str(shotn),', iregist failed'];
    disp(message);
end

if flg_save == 1
    registfilename = ['dytrans_ts_9@',num2str(shotn),'.dat'];
     delete(registfilename); %for egcalc
    
end

end

% !git checkout master
% !git branch -a
% !git add *.m
% !git commit -m ''
% !git push origin master
