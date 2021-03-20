function mmeshTi = mapTi(data_cx7_ex, time_cx7_ex, time_tsmap, timen_tsmap)
%eg登録データフォーマット この関数の概要をここに記述
%   詳細説明をここに記述

size_data_cx7_ex  = size(data_cx7_ex);
mmeshTi           = zeros(timen_tsmap,size_data_cx7_ex(2),size_data_cx7_ex(3));

size_time_cx7_ex = size(time_cx7_ex);
end_cx7_ex = size_time_cx7_ex(2);

[a,b] = min(abs(time_tsmap-time_cx7_ex(1)));
[c,d] = min(abs(time_tsmap-time_cx7_ex(end_cx7_ex)));

for j=1:size_data_cx7_ex(3)
    for k=1:size_data_cx7_ex(2)
       
        mmeshTi(:,k,j) = interp1(time_cx7_ex,data_cx7_ex(:,k,j),time_tsmap);
        
    end
end

if time_tsmap(b) < time_cx7_ex(1)
   b = b + 1;
end

if time_tsmap(d) > time_cx7_ex(end_cx7_ex)
   d = d - 1;
end


mmeshTi(b-1,:,:) = mmeshTi(b,:,:); %前後一フレームをおまけ
mmeshTi(d+1,:,:) = mmeshTi(d,:,:); %前後一フレームをおまけ

end
