for day = 1:4
    for time = 1:288
        price_noise = zeros(1,24);
        outputname = strcat('price_noise/price_noise_aug0',...
             num2str(day),'_time',num2str(time),'.mat')
        
        predict_time = 1;
        for predict_std = 1:9/23:10
            price_noise(1, predict_time)=...
                randn(1)*predict_std;
            predict_time = predict_time+1;
        end
        save(outputname, 'price_noise')
    end
end