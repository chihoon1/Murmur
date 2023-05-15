function h_hat = MomentMatchHurst_new(data, filt, L,a,b, ismean)

    J = log2(size(data,1));
    wddata = dwtr(data, J - L, filt);
    pairs = nchoosek( 1 : J - L, 2); 
    
    pairs = pairs(find( pairs(:,1) >= a & pairs(:,2 ) <=b ),:);
    
    H_hat = zeros(size(pairs, 1), 2);
    
    for p = 1: size(pairs,1)
        k1 = pairs(p,1); k2 = pairs(p,2); 

        l1 = J - k2; l2 = J - k1;
        k1_indx =  2^(l1) + 1 : 2^(l1 + 1); k2_indx =  2^(l2) + 1 : 2^(l2 + 1);
                
        help_k1 = wddata( k1_indx ); help_k2 = wddata( k2_indx );
                
        if  ismean == 1
            d_k1 = mean ( help_k1.^2);  
            d_k2 = mean ( help_k2.^2);
                    
        elseif  ismean == 0
            d_k1 = median ( help_k1.^2);  
            d_k2 = median ( help_k2.^2);
                    
        elseif  ismean == 2
            d_k1 = fastDcov(help_k1', help_k1');  
            d_k2 = fastDcov(help_k2', help_k2');
        end 
                
        d_k1 = log2(d_k1); d_k2 = log2(d_k2);

        A = d_k1 - d_k2; 
        B = psi( 2^( l1 - 1 ) ) - psi( 2^(l2 - 1) ) - log(2^(l1)) + log( 2^(l2) );
        B = log2( exp(1) ) * B;
        C = l1 - l2;
        
        h_p = ( 1/ ( 2*C ) ) * ( B - A  - C ); 
        h_p_var = ( 1/ ( 2^(2*C) ) ) * ( psi( 1, 2^( l1 - 2 ) ) - psi( 1, 2^(l2 - 2) ) );
                
        H_hat(p, :) =  [ h_p h_p_var ];
    end
      wgt = H_hat(:, 2)/sum(H_hat(:, 2));
            
      h_hat = H_hat(:,1)'* wgt;
end 
