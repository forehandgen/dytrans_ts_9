clearvars -except shotn

%shotn = 166256; %Fedelico

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
            %array 9抽出
            cx9_sele   = (cx9_0(:,3) == 9); 
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