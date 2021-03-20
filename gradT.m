function [flg,grad_T_linear] = gradT(rhon,reff,T,val_ave)

    % tiledlayout(3,4);
    % for k=1:2:20
try
        val_ave = 5; %前後何点を選ぶか
        grad_T_linear = zeros(rhon,1);

        for j=(1+val_ave):(rhon-val_ave)

            y = T(j-val_ave:j+val_ave);
            x = reff(j-val_ave:j+val_ave);
            p = polyfit(x,y,1);
            grad_T_linear(j)=p(1);

        end

    flg = 1;
        % nexttile;
        % plot(reff,grad_Te_linear,'-')
        % title("points="+num2str(k));
        % ylim=([0,-15]);
        % hold on;
        % plot(reff,dTedreff,'-');
        % hold off;


    %end
catch
    flg = 0;
end
end