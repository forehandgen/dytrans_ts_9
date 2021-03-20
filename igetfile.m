function flg = igetfile(shotn,diagname)
%UNTITLED2 この関数の概要をここに記述
%   詳細説明をここに記述
command1 = 'export PATH="/Users/gmotojima/Documents/Research/igetfile/mac:$PATH";';
command2 = 'export CLIENT_INI2="/Users/gmotojima/Documents/Research/igetfile/client2.ini";';
command3 = 'export INDEXSERVERNAME="dasindex.lhd.nifs.ac.jp";igetfile -s ';
command4 = [' -d ',diagname,' -o '];
command5 = ['ncm/',diagname,'@'];
command6 = '.dat';
%diagname = 'tsmap_calib';
%shotn    = 162111;

allcommand = [command1,command2,command3,num2str(shotn),command4,command5,num2str(shotn),command6];
try
    system(allcommand);
    flg = 1;
catch
    flg = 0;
end
end

