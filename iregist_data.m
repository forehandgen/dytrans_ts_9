function flg = iregist_data(signame,ishot)

    %ishot = 169623;
    filename = [signame,'@',num2str(ishot),'.dat'];

    %For eg regist and/or idelete
    try
                 flg = 1;         
                 command_del =  ['idelete -s ',num2str(ishot),' -m 1 -d dytrans_ts_9'];
                 status = system(command_del);        
                 command = ['iregist -s ' ,num2str(ishot), ' -m 1 -d dytrans_ts_9 -p ' ,filename , ' -u kc-motojima -w tN8P5_De'];
                 status = system(command);
    catch
        warning(strcat('something wrong! iregist failed? ',signame,'@',num2str(ishot)));
        flg = 0;
    end
end
