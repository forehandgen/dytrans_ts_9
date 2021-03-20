
function aaQ_out1 = mask_data(aaQ_out,timen_cx7_ex,time_cx7_ex)

%aaQ_out = mapTi(data_cx7_ex, time_cx7_ex, time_tsmap, timen_tsmap); %Tiをtsmapにマッピング
s_data   = size(aaQ_out);
aaQ_out1 = aaQ_out; 

          aaQ_out1(:,:,7)     = 0;
          aaQ_out1(:,:,11:12) = 0;
          aaQ_out1(:,:,14:16) = 0;
          aaQ_out1(:,:,20:21) = 0;
          aaQ_out1(:,:,23:27) = 0;
          aaQ_out1(:,:,29)    = 0;
          aaQ_out1(:,:,31:33) = 0;
          aaQ_out1(:,:,35:41) = 0;
          


    for j=1:timen_cx7_ex 
    t1 = time_cx7_ex(1,j);

    [val_cx7,idx] = min(abs(aaQ_out(:,1,1)-t1));

        if val_cx7 > 1/30   %cx7のタイミングとTSのタイミングが30ms以内に収まるものを選ぶ
  
        else
%j,idx
            if idx ~= 1
                %map_sele(idx-1,:,:) = aaQ_out(idx-1,:,:); %一つ前のフレームをおまけ
            end

                aaQ_out1(idx,:,:)   = aaQ_out(idx,:,:);

            if idx ~= timen_cx7_ex
                %map_sele(idx+1,:,:) = aaQ_out(idx+1,:,:); %一つ後のフレームをおまけ
            end
        end
    end
    
end